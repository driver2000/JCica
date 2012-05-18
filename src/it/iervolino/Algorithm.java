/**
 * Contiene gli algoritmi
 */
package it.iervolino;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import it.iervolino.extern.HungarianAlgorithm;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.ComposableFunction;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.correlation.Covariance;

/**
 * Gestisce tutti i passi dell'algoritmo CICA
 * @author Daniele Iervolino
 * @version 1.0
 */
public class Algorithm {

  private static ArrayList<RealMatrix> whitening(RealMatrix workdb) {

    //Matrice di covarianza
    DoubleFactory2D factory = DoubleFactory2D.dense;
    Covariance c = new Covariance(workdb.transpose());
    DoubleMatrix2D cov = factory.make(c.getCovarianceMatrix().getData());

    //Autovalori e autovettori
    DenseDoubleAlgebra algebra = DenseDoubleAlgebra.DEFAULT;
    DenseDoubleEigenvalueDecomposition eig = algebra.eig(cov);
    DoubleMatrix2D Dx = eig.getD();

    RealMatrix Ex = new Array2DRowRealMatrix(eig.getV().toArray());

    //Whitening  
    DoubleMatrix2D invValue = algebra.inverse(Dx);
    for (int i = 0; i < invValue.rows(); i++) {
      invValue.set(i, i, Math.sqrt(invValue.get(i, i)));
    }

    RealMatrix invDx = new Array2DRowRealMatrix(invValue.toArray());
    RealMatrix Q = invDx.multiply(Ex.transpose());
    RealMatrix X = Q.multiply(workdb);

    ArrayList<RealMatrix> xq = new ArrayList<RealMatrix>(2);
    xq.add(0, X);
    xq.add(1, Q);
    return xq;
  }

  private static void amplitudeIndetermincy(ArrayList<RealMatrix> yTot, RealMatrix Q, RealMatrix W, int numSig, int index) {
    int size = yTot.get(0).getColumnDimension();

    RealMatrix invQ = new LUDecompositionImpl(Q).getSolver().getInverse();
    RealMatrix invW = new LUDecompositionImpl(W).getSolver().getInverse();
    RealMatrix mat = new Array2DRowRealMatrix(numSig, size);

    for (int i = 0; i < numSig; i++) {
      mat.setRowVector(i, yTot.get(index).getRowVector(i));
      RealVector means = meanCols(invQ.multiply(invW).multiply(mat));
      yTot.get(index).setRowVector(i, means);
      mat = mat.scalarMultiply(0.0);
    }
  }

  private static ArrayList<RealMatrix> createRicDP(ArrayList<RealMatrix> yTot, int numSig, int n) {

    //Costruzione DB di lavoro
    ArrayList<RealMatrix> db = new ArrayList<RealMatrix>();

    int m = yTot.get(0).getColumnDimension();
    RealMatrix workdb = new Array2DRowRealMatrix(n, m);

    RealMatrix workfft = null;
    for (int i = 0; i < numSig; i++) {
      for (int j = 0; j < n; j++) {
        workfft = yTot.get(j);
        workdb.setRowVector(j, workfft.getRowVector(i));

      }

      db.add(i, workdb);
      workdb = workdb.scalarMultiply(0.0);
    }

    return db;
  }

  private static double sumAbs(RealVector wOld, RealVector w) {
    double sum = 0.0;
    double a = 0.0;
    double b = 0.0;

    for (int i = 0; i < w.getDimension(); i++) {
      a = w.getEntry(i);
      b = wOld.getEntry(i);
      sum += Math.pow(Math.abs(a) - Math.abs(b), 2);
    }
    return Math.sqrt(sum);
  }

  private static RealVector decorrelation(RealMatrix W, RealVector w) {
    RealVector work = w.subtract(W.multiply(W.transpose()).operate(w));
    work = work.mapDivide(work.getNorm());
    return work;
  }

  private static RealVector tanh(RealMatrix X, RealVector w) {
    RealVector y = X.operate(w);
    for (int i = 0; i < y.getDimension(); i++) {
      double value = 1 - 2 / (Math.exp(2 * y.getEntry(i)) + 1);
      y.setEntry(i, value);
    }
    return y;
  }

  private static RealVector meanCols(RealMatrix mat) {
    RealVector means = new ArrayRealVector(mat.getColumnDimension());
    for (int j = 0; j < means.getDimension(); j++) {
      means.setEntry(j, StatUtils.mean(mat.getColumn(j)));
    }
    return means;
  }

  private static ArrayList<RealMatrix> hungDP(ArrayList<RealMatrix> ricDP) {

    int n = ricDP.get(0).getRowDimension();
    int col = ricDP.get(0).getColumnDimension();
    int m = ricDP.size();
    RealMatrix cost = new Array2DRowRealMatrix(m, m);

    RealMatrix xin = new Array2DRowRealMatrix(m, col);

    ArrayList<RealMatrix> xric = new ArrayList<RealMatrix>(m);

    for (int i = 0; i < m; i++) {
      xric.add(i, new Array2DRowRealMatrix(n, col));
      for (int j = 0; j < col; j++) {
        xin.setEntry(i, j, ricDP.get(i).getEntry(0, j));
        xric.get(i).setEntry(0, j, ricDP.get(i).getEntry(0, j));
      }
    }

    RealVector work = null;
    for (int i = 1; i < n; i++) {
      for (int j = 0; j < m; j++) {


        RealVector workJ = xin.getRowVector(j).ebeMultiply(xin.getRowVector(j));
        double sum = StatUtils.sum(workJ.getData());
        workJ = workJ.mapDivide(sum);

        for (int k = 0; k < m; k++) {

          RealMatrix ricDPkesima = ricDP.get(k);
          work = ricDPkesima.getRowVector(i);
          work = work.ebeMultiply(work);
          double sumK = StatUtils.sum(work.getData());
          work = work.mapDivide(sumK);

          double costK = crossEntropy(workJ, work);
          work = work.mapMultiplyToSelf(0.0);
          cost.setEntry(j, k, costK);
        }
      }

      int[] C = HungarianAlgorithm.hgAlgorithm(cost.getData(), "min");
      for (int k = 0; k < m; k++) {
        for (int j = 0; j < col; j++) {
          xric.get(k).setEntry(i, j, ricDP.get(C[k]).getEntry(i, j));
          xin.setEntry(k, j, ricDP.get(C[k]).getEntry(i, j));
        }
      }
    }

    return xric;
  }

  private static double crossEntropy(RealVector workJ, RealVector work) {

    double cost = 0.0;
    double tiny = Math.exp(-700);

    RealVector p = new ArrayRealVector(workJ.getDimension());
    p = workJ.mapAdd(tiny);
    p = p.ebeDivide(work.mapAdd(tiny));
    try {
      p = p.map(ComposableFunction.LOG);
    } catch (FunctionEvaluationException ex) {
      Logger.getLogger(Algorithm.class.getName()).log(Level.SEVERE, null, ex);
    }
    cost = StatUtils.sum(workJ.ebeMultiply(p).getData());
    return cost;
  }

  /**
   * Genera una matrice FFT per ogni segnale di input
   * @param signal il segnale
   * @param points i punti della FFT
   * @param window dimensione della finestra della FFT
   * @param overlap sovrapposizione delle finistre
   * @return la matrice FFT 
   */
  public static RealMatrix STFT(RealVector signal, int points, int window, int overlap) {

    int sizeColumn = signal.getDimension();

    //controllo finestra dispari
    if (window % 2 == 0) {
      window++;
    }

    //lunghezza dimezzata della finestra 
    int halfLen = (window - 1) / 2;

    //punto medio della finestra
    int midPoint = points / 2;
    int actminLen = Math.min(halfLen, midPoint);

    //finestra dimezzata
    RealVector halfWin = new ArrayRealVector(halfLen + 1);

    //inizializzazione finestra
    for (int i = 0; i < halfWin.getDimension(); i++) {
      double d = (double) i / halfLen;
      double value = 0.5 * (1 + Math.cos(Math.PI * d));
      halfWin.setEntry(i, value);
    }


    //particolare istanza della finestra per fft
    RealVector newWindow = new ArrayRealVector(points, 0.0);
    newWindow.setSubVector(midPoint, halfWin.getSubVector(0, actminLen));
    double[] vec = halfWin.getSubVector(0, actminLen).getData();

    for (int i = 0; i < vec.length; i++) {
      newWindow.setEntry(midPoint - i, vec[i]);
    }

    //matrice della STFT
    RealMatrix cMatrix = new Array2DRowRealMatrix(midPoint + 1, 1 + (sizeColumn - points) / overlap);

    int col = 0;
    double[] t = null;

    for (int i = 0; i < (sizeColumn - points); i += overlap) {
      RealVector u = newWindow.ebeMultiply(signal.getSubVector(i, points));

      RealVector work = new ArrayRealVector(u.getDimension() * 2, 0.0);
      work.setSubVector(0, u);
      //FFT
      DoubleFFT_1D fft = new DoubleFFT_1D(u.getDimension());
      t = work.getData();
      fft.realForwardFull(t);

      for (int j = 0; j <= midPoint; j++) {
        cMatrix.setEntry(j, col, t[j * 2]);
      }
      col++;
    }
    return cMatrix;
  }

  /**
   * Effettua la ISTFT dei dati di input, e restituisce i segnali separati
   * @param mat matrice FFT
   * @param points i punti della IFFT
   * @param window dimensione della finestra della IFFT
   * @param overlap sovrapposizione delle finistre
   * @return un singolo segnale
   */
  public static RealVector ISTFT(RealMatrix mat, int points, int window, int overlap) {

    int cols = mat.getColumnDimension();

    int xlen = points + (cols * overlap);


    //controllo finestra dispari
    if (window % 2 == 0) {
      window++;
    }

    //lunghezza dimezzata della finestra 
    int halfLen = (window - 1) / 2; // 256


    //punto medio della finestra
    int midPoint = points / 2; // 256
    int actminLen = Math.min(halfLen, midPoint); // 256

    //finestra dimezzata
    RealVector halfWin = new ArrayRealVector(halfLen + 1); //257

    //inizializzazione finestra
    for (int i = 0; i < halfWin.getDimension(); i++) {
      double d = (double) i / halfLen;
      double value = 0.5 * (1 + Math.cos(Math.PI * d));
      halfWin.setEntry(i, value);

    }

    //particolare istanza della finestra per fft
    RealVector newWindow = new ArrayRealVector(points, 0.0); // 512

    newWindow.setSubVector(midPoint, halfWin.getSubVector(0, actminLen));
    double[] vec = halfWin.getSubVector(0, actminLen).getData();


    for (int i = 0; i < vec.length; i++) {
      newWindow.setEntry(midPoint - i, vec[i]);
    }


    RealVector x = new ArrayRealVector(xlen, 0.0);
    double[] w = null;
    for (int i = 0; i < cols * overlap; i += overlap) {

      RealVector ft = mat.getColumnVector(i / overlap);

      RealVector ftPadd = new ArrayRealVector(ft.getDimension() + midPoint - 1);

      ftPadd.setSubVector(0, ft);

      for (int j = 0; j < midPoint - 1; j++) {
        ftPadd.setEntry(ft.getDimension() - 1 + j, ft.getEntry(midPoint - j));
      }

      //Nuova IFFT
      DoubleFFT_1D ifft = new DoubleFFT_1D(ftPadd.getDimension());
      RealVector work = new ArrayRealVector(ftPadd.getDimension() * 2, 0.0);
      work.setSubVector(0, ftPadd);
      w = work.getData();
      ifft.realInverseFull(w, true);

      RealVector px = new ArrayRealVector(w.length / 2);
      for (int j = 0; j < px.getDimension(); j++) {
        px.setEntry(j, w[j * 2]);
        // System.out.println(px.getEntry(j));
      }
      //px.*win
      RealVector mul = px.ebeMultiply(newWindow);

      //x((b+1):(b+ftsize)) = x((b+1):(b+ftsize)) + mul
      x.setSubVector(i, x.getSubVector(i, points).add(mul));

    }

    return x;
  }

  /**
   * Questo metodo costruisce un Database delle righe delle matrici contenute in  shortTimeFourier
   * @param shortTimeFourier contiene le matrici delle stft
   * @param numSig numero dei segnali
   * @param n rappresenta la lunghezza dei segnali
   * @return il database
   */
  public static ArrayList<RealMatrix> createDB(ArrayList<RealMatrix> shortTimeFourier, int numSig, int n) {
    //Costruzione DB di lavoro
    ArrayList<RealMatrix> db = new ArrayList<RealMatrix>();
    int m = shortTimeFourier.get(0).getColumnDimension();
    RealMatrix workdb = new Array2DRowRealMatrix(numSig, m);

    //Matrice di lavoro
    RealMatrix workfft = null;
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < numSig; i++) {
        workfft = shortTimeFourier.get(i);
        workdb.setRowVector(i, workfft.getRowVector(j));
      }
      db.add(j, workdb);
      workdb = workdb.scalarMultiply(0.0);
    }
    return db;
  }

  /**
   * Si occupa della decorrelazione dei segnali di input
   * @param db database contente le righe delle matrici della FFT
   * @param size dimensione delle righe delle matrici FFT
   * @param numSig numero segnali
   * @param cols lunghezza dei segnali
   * @return  matrici di segnali decorrelati
   */
  public static ArrayList<RealMatrix> computeDP(ArrayList<RealMatrix> db, int size, int numSig, int cols) {

    ArrayList<RealMatrix> yTot = new ArrayList<RealMatrix>(size);
    for (int index = 0; index < size; index++) {
      RealMatrix workdb = db.get(index);
      ArrayList<RealMatrix> xq = null;

      xq = whitening(workdb);



      RealMatrix X = xq.get(0);
      RealMatrix Q = xq.get(1);

      //Matrice di lavoro per operazione di decorrelazione
      RealMatrix W = new Array2DRowRealMatrix(numSig, numSig);

      int maxCounter = 1000;

      //Algoritmo per il calcolo delle componenti
      for (int i = 0; i < numSig; i++) {

        int count = 0;
        //Decorrelazione
        RealVector w = new ArrayRealVector(numSig);

        for (int j = 0; j < w.getDimension(); j++) {
          w.setEntry(i, Math.random() - 0.5);
        }

        w = decorrelation(W, w);

        RealVector wOld = new ArrayRealVector(numSig, 0.0);
        double tolerance = 0.00001;

        while (Math.min(sumAbs(wOld, w), maxCounter - count) > tolerance) {
          wOld = w;
          RealVector hypTan = tanh(X.transpose(), w);
          RealVector XhypTan = X.operate(hypTan);
          double sum = 0.0;

          for (int j = 0; j < hypTan.getDimension(); j++) {
            sum += 1 - Math.pow(hypTan.getEntry(j), 2);
          }

          w = XhypTan.subtract(w.mapMultiply(sum)).mapDivide(cols);
          w = decorrelation(W, w);
          count++;
        }
        W.setColumnVector(i, w);
      }
      yTot.add(index, W.transpose().multiply(X));
      amplitudeIndetermincy(yTot, Q, W.transpose(), numSig, index);
    }

    ArrayList<RealMatrix> ricDP = createRicDP(yTot, numSig, size);
    ricDP = hungDP(ricDP);

    return ricDP;
  }

  /**
   * Calcola lo spettro di potenza dei segnali separati
   * @param mat matrice dei segnali
   * @return lo spettro di potenza
   */
  public static RealMatrix spettro(RealMatrix mat) {

    int numSig = mat.getRowDimension();
    int cols = mat.getColumnDimension();

    double re = 0.0;
    double value = 0.0;
    double[] tc = null;

    DoubleFFT_1D fft = new DoubleFFT_1D(cols);
    RealMatrix spettro = new Array2DRowRealMatrix(numSig, cols / 2);
    RealVector vecTc = new ArrayRealVector(cols * 2, 0.0);
    for (int i = 0; i < numSig; i++) {
      vecTc.setSubVector(0, mat.getRow(i));
      tc = vecTc.getData();
      fft.realForwardFull(tc);
      for (int j = 0; j < tc.length / 4; j++) {
        re = tc[2 * j];
        value = re * re;
        spettro.setEntry(i, j, value);
      }
    }
    return spettro;
  }

  public static void print(RealMatrix asd) {
    System.out.println(asd.getRowDimension() + " x " + asd.getColumnDimension());

    for (int i = 0; i < asd.getRowDimension(); i++) {
      for (int j = 0; j < asd.getColumnDimension(); j++) {
        System.out.print(asd.getEntry(i, j) + " ");
      }
      System.out.println("");

    }
    System.out.println("");
  }
}
