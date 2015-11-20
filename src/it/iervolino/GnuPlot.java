package it.iervolino;
// import stuff
import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.layout.StripeLayout;
import com.panayotis.gnuplot.plot.AbstractPlot;
import com.panayotis.gnuplot.style.NamedPlotColor;
import com.panayotis.gnuplot.style.PlotStyle;
import com.panayotis.gnuplot.style.Style;

import com.panayotis.gnuplot.terminal.ImageTerminal;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;

/**
 * Continene i metodi per creare i plot dei segnali e dei relativi spettri
 * @author Daniele Iervolino
 * @version 1.0
 */
public class GnuPlot {

    private JavaPlot plotter;
    private RealMatrix matrix;
    private RealMatrix data;

    public GnuPlot(RealMatrix data, String type) throws IOException {

        if (type == null) {
            throw new IOException("Tipo plot non inserito");
        }
        
        this.plotter = new JavaPlot();
        this.data = data;
        matrix = null;

        ImageTerminal im = new ImageTerminal();
        im.set("size", "700,900");
        plotter.getAxis("y").set("grid");
        plotter.getAxis("x").set("grid", "noxtics");


        if (type.equals("Spettro")) {
            plotSpettro(data.getColumnDimension());

        } else {
            plot(data.getColumnDimension(), data.getRowDimension());
        }

        this.plotter.setTerminal(im);
        this.plotter.set("lmargin", "8");
        this.plotter.set("rmargin", "5");
        this.plotter.set("tmargin", "1");
        this.plotter.set("tmargin", "1");
        this.plotter.setKey(JavaPlot.Key.OFF);
        //this.plotter.setMultiTitle(type);

        StripeLayout lo = new StripeLayout();
        lo.setType(StripeLayout.EXPANDROWS);
        this.plotter.getPage().setLayout(lo);
        this.plotter.plot();

        BufferedImage bi = im.getImage();
        File outputfile = new File("/home/daniele/NetBeansProjects/ApdpWeb/web/resources/images/" + type + ".png");
        ImageIO.write(bi, "png", outputfile);

    }

    private void plot(int size, int numSig) {

        RealVector index = new ArrayRealVector(size);
        for (int i = 0; i < size; i++) {
            index.setEntry(i, i);
        }

        int indexSize = index.getDimension();

        RealMatrix mat = new Array2DRowRealMatrix(indexSize, numSig);
        mat.setColumnVector(0, index);
        for (int i = 1; i <= mat.getColumnDimension(); i++) {
            if (i != 1) {
                plotter.newGraph();
            }
            double ytics = (StatUtils.max(data.getRowVector(i - 1).getData()) - StatUtils.min(data.getRowVector(i - 1).getData()));
            ytics /= 7.0;
            mat.setColumnVector(1, data.getRowVector(i - 1));
            plotter.addPlot(mat.getData());
            PlotStyle stl = ((AbstractPlot) plotter.getPlots().get(0)).getPlotStyle();
            stl.setStyle(Style.LINES);
            stl.setLineType(NamedPlotColor.ROYALBLUE);
            this.plotter.getAxis("y").set("ytics", "" + ytics);
            this.plotter.getAxis("y").set("format y", "\"%.3g\"");
        }
    }

    private void plotSpettro(int size) {

        RealVector index = new ArrayRealVector();
        for (int i = 0; i < size; i++) {
            double value = (double) i * (62.5 / (double) (size));
            if (value <= 2) {

                index = index.append(value);
            } else {
                break;
            }
        }

        int indexSize = index.getDimension();
        matrix = new Array2DRowRealMatrix(indexSize, data.getRowDimension());
        matrix.setColumnVector(0, index);

        for (int i = 0; i < matrix.getColumnDimension(); i++) {
            RealVector vec = data.getRowVector(i).getSubVector(0, indexSize);
            matrix.setColumnVector(i, vec);
        }

        RealMatrix mat = new Array2DRowRealMatrix(indexSize, 2);
        mat.setColumnVector(0, index);
        for (int i = 1; i <= matrix.getColumnDimension(); i++) {
            if (i != 1) {
                plotter.newGraph();
            }
            mat.setColumnVector(1, matrix.getColumnVector(i - 1));
            plotter.addPlot(mat.getData());
            PlotStyle stl = ((AbstractPlot) plotter.getPlots().get(0)).getPlotStyle();
            stl.setStyle(Style.LINES);
            stl.setLineType(NamedPlotColor.GOLDENROD);
            this.plotter.getAxis("x").setLabel("Frequency (Hz)", "Monaco", 8);
//            double ytics = (StatUtils.max(data.getRowVector(i - 1).getData()) - StatUtils.min(data.getRowVector(i - 1).getData()));
//            ytics /= 7.0;
            this.plotter.getAxis("x").set("xtics", "0.25");

            this.plotter.getAxis("x").setBoundaries(0, 3);
            this.plotter.getAxis("y").setBoundaries(StatUtils.min(matrix.getColumnVector(i - 1).getData()),
                    StatUtils.max(matrix.getColumnVector(i - 1).getData()));
            //  this.plotter.getAxis("y").set("ytics", "" + ytics);

        }
    }
}
