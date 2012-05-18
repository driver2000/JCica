package it.iervolino;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Crea i plot dei dati
 * @param matrix matrice contente i dati da plottare
 * @param type Spettro se si vuole lo spettro, altrimenti genera un plot dei segnali separati
 */
public class Plot {

  private double[][] mat;
  private CombinedDomainXYPlot plot;
  private NumberAxis range;
  private XYPlot[] data;

  public Plot(double[][] mat, String type) {

    this.mat = mat;
    if (type.equals("Spettro")) {
      data = createDatasetSpettro();
      range = new NumberAxis("Frequenza (Hz)");
    } else {
      data = createDataset();
      range = new NumberAxis("Campioni");
    }

    plot = new CombinedDomainXYPlot(range);
    for (int i = 0; i < data.length; i++) {
      plot.add(data[i]);
    }
    plot.setOrientation(PlotOrientation.VERTICAL);
  }

  private XYPlot[] createDataset() {
    int cols = mat[0].length;
    int rows = mat.length;

    XYItemRenderer renderer = new StandardXYItemRenderer();
    Color myBlue = new Color(20, 125, 182);
    renderer.setSeriesPaint(0, myBlue);
    NumberAxis rangeAxis = new NumberAxis();
    XYPlot[] xyplot = new XYPlot[rows];

    for (int i = 0; i < rows; i++) {
      XYSeries series = new XYSeries(" ");
      for (int j = 0; j < cols; j++) {
        series.add(j, mat[i][j]);
      }
      xyplot[i] = new XYPlot(new XYSeriesCollection(series), null, rangeAxis, renderer);

    }
    return xyplot;
  }

  private XYPlot[] createDatasetSpettro() {
    int cols = mat[0].length;
    int rows = mat.length;
    double index = 0.0;
    XYItemRenderer renderer = new StandardXYItemRenderer();
    renderer.setSeriesPaint(0, Color.DARK_GRAY);

    XYPlot[] xyplot = new XYPlot[rows];
    for (int i = 0; i < rows; i++) {
      XYSeries series = new XYSeries(" ");
      NumberAxis rangeAxis = new NumberAxis();
      rangeAxis.setAutoRangeIncludesZero(false);
      rangeAxis.setVisible(false);


      for (int j = 0; j < cols; j++) {
        index = (double) j * (62.5 / (double) (cols));
        if (index <= 2.3) {
          series.add(index, mat[i][j]);

        } else {
          break;
        }
      }
      xyplot[i] = new XYPlot(new XYSeriesCollection(series), null, rangeAxis, renderer);
      xyplot[i].setRangeGridlinesVisible(false);
      xyplot[i].setRangeCrosshairVisible(false);
    }
    return xyplot;
  }

  public void plotToFile(String fileName) throws IOException {
    JFreeChart chart = new JFreeChart(null, null, plot, false);
    ChartUtilities.saveChartAsPNG(new File(fileName), chart, 700, 600);
  }
}
