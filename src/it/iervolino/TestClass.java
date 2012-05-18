package it.iervolino;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

public class TestClass {

  public static void main(String[] args) {

    int points = 512;
    int window = 512;
    int overlap = 128;
    String fileName[] = new String[4];
//        fileName[0] = "SIN";
//        fileName[1] = "FUNNY";

    fileName[0] = "AMS2_E_18.in";
    fileName[1] = "ASB2_E_18.in";
    fileName[2] = "BGNG_E_18.in";
    fileName[3] = "TAGG_E_18.in";

    RealMatrix mix = new Array2DRowRealMatrix(fileName.length, 450000);
    String str = null;
    BufferedReader reader = null;
    try {
      for (int i = 0; i < mix.getRowDimension(); i++) {

        File doc = new File("/home/daniele/Database/" + fileName[i]);
        reader = new BufferedReader(new FileReader(doc));
        while ((str = reader.readLine()) != null) {

          StringTokenizer token = new StringTokenizer(str, ";");
          int j = 0;
          while (token.hasMoreTokens()) {
            mix.setEntry(i, j, Double.parseDouble(token.nextToken()));
            ++j;
          }

        }
      }
      reader.close();
    } catch (IOException ex) {
      Logger.getLogger(TestClass.class.getName()).log(Level.SEVERE, null, ex);
    }


    for (String string : fileName) {
      System.out.println(string);
    }

    JCica cica = new JCica(points, window, overlap, mix);
    RealMatrix demixSignals = cica.demix();

    RealMatrix spettro = cica.powerSpectrum();
    Plot plot = new Plot(demixSignals.getData(), "Segnali Separati");
    Plot spect =  new Plot(spettro.getData(), "Spettro");
   
    try {
//      plot.plotToFile(new File("/home/daniele/Segnali Separati.png"));
//      spect.plotToFile(new File("/home/daniele/Spettro.png"));
      cica.saveToFile("/home/daniele/Scrivania/Tesi/Daniele/demix.out");

    } catch (IOException ex) {
      Logger.getLogger(TestClass.class.getName()).log(Level.SEVERE, null, ex);
    }
    System.exit(0);
  }
}
