package it.iervolino;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.StatUtils;

/**
 * Classe wrapper che permette di effettuare la separazione dei segnali
 * convoluti, effettuando i passi dell'algoritmo CICA
 * 
 * @author Daniele Iervolino
 * @version 1.0
 */
public class JCica {

  private int points;
  private int window;
  private int overlap;
  private int rows;
  private int columns;
  private RealMatrix mixMatrix;
  private RealMatrix demixMatrix;

  /**
   * Costruttore
   * @param points punti FFT
   * @param window lunghezza della finestra per le FFT
   * @param overlap sovrapposizione finestre
   * @param mixMatrix matrice contenete i segnali da separare
   */
  public JCica(int points, int window, int overlap, RealMatrix mixMatrix) {
    this.points = points;
    this.window = window;
    this.overlap = overlap;
    this.mixMatrix = mixMatrix;
    this.rows = mixMatrix.getRowDimension();
    this.columns = mixMatrix.getColumnDimension();

  }

  /**
   * Effettua la separazione dei segnali, richiamando i passi dell'algoritmo CICA
   * @see Algorithm
   * @return segnali separati
   */
  public RealMatrix demix() {

    // Normalizzazione 
    for (int i = 0; i < rows; i++) {
      mixMatrix.setRow(i, StatUtils.normalize(mixMatrix.getRow(i)));
    }

    //Short Time Fourier Transform
    ArrayList<RealMatrix> stft = new ArrayList<RealMatrix>();
    for (int i = 0; i < rows; i++) {
      stft.add(i, Algorithm.STFT(mixMatrix.getRowVector(i), points, window, overlap));
    }

    int sizeStft = stft.get(0).getRowDimension();

    ArrayList<RealMatrix> database = Algorithm.createDB(stft, rows, sizeStft);
    ArrayList<RealMatrix> xricDP = Algorithm.computeDP(database, sizeStft, rows, columns);

    columns = xricDP.get(0).getColumnDimension();
    columns = points + (columns * overlap);

    //Inverse Short Time Fourier Transform
    demixMatrix = new Array2DRowRealMatrix(rows, columns);
    for (int i = 0; i < rows; i++) {
      demixMatrix.setRowVector(i, Algorithm.ISTFT(xricDP.get(i), points, window, overlap));
    }

    return demixMatrix;
  }

  /**
   * Genera lo spettro dei sengnali separati
   * @see Algorithm
   * @return lo spettro dei sengali separati
   */
  public RealMatrix powerSpectrum() {

    return Algorithm.spettro(demixMatrix);

  }

  /**
   * Genera un file contenete i dati dei sengali separati
   * @param filename path del file 
   * @throws IOException 
   */
  public void saveToFile(String filename) throws IOException {
    File file = new File(filename);
    BufferedWriter writer = new BufferedWriter(new FileWriter(file));
    try {
      for (int i = 0; i < demixMatrix.getRowDimension(); i++) {
        double[] row = demixMatrix.getRow(i);
        for (int j = 0; j < row.length; j++) {
          writer.write(row[j] + ";");
        }
        writer.newLine();
      }

    } finally {
      writer.close();
    }
  }
}
