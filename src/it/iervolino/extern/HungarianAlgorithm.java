package it.iervolino.extern;

/**
 * Created on Apr 25, 2005
 * 
 * Munkres-Kuhn (Hungarian) Algorithm Clean Version: 0.11
 * 
 * 
 * This Java class implements the Hungarian algorithm [a.k.a Munkres' algorithm,
 * a.k.a. Kuhn algorithm, a.k.a. Assignment problem, a.k.a. Marriage problem,
 * a.k.a. Maximum Weighted Maximum Cardinality Bipartite Matching].
 *  <p>
 * [It can be used as a method call from within any main (or other function).]
 * It takes 2 arguments:
 * a. A 2-D array (could be rectangular or square).
 * b. A string ("min" or "max") specifying whether you want the min or max assignment.
 * [It returns an assignment matrix[array.length][2] that contains the row and col of
 * the elements (in the original inputted array) that make up the optimum assignment.]
 *  <p>
 * [This version contains only scarce comments. If you want to understand the 
 * inner workings of the algorithm, get the tutorial version of the algorithm
 * from the same website you got this one (http://www.spatial.maine.edu/~kostas/dev/soft/munkres.htm)]
 * 
 * Any comments, corrections, or additions would be much appreciated. 
 * Credit due to professor Bob Pilgrim for providing an online copy of the
 * pseudocode for this algorithm (http://216.249.163.93/bob.pilgrim/445/munkres.html)
 * 
 * Feel free to redistribute this source code, as long as this header--with
 * the exception of sections in brackets--remains as part of the file.
 * 
 * Requirements: JDK 1.5.0_01 or better.
 * [Created in Eclipse 3.1M6 (www.eclipse.org).]
 * 
 * 
 * @author Konstantinos A. Nedas                     
 * Department of Spatial Information Science & Engineering
 * University of Maine, Orono, ME 04469-5711, USA
 * kostas@spatial.maine.edu
 * http://www.spatial.maine.edu/~kostas   
 * @version 0.11
 *
 **/
public class HungarianAlgorithm {

  //Finds the largest element in a positive array.
  //works for arrays where all values are >= 0.
  public static double findLargest(double[][] array) {
    double largest = 0;
    for (int i = 0; i < array.length; i++) {
      for (int j = 0; j < array[i].length; j++) {
        if (array[i][j] > largest) {
          largest = array[i][j];
        }
      }
    }

    return largest;
  }

  //Transposes a double[][] array.
  public static double[][] transpose(double[][] array) {
    double[][] transposedArray = new double[array[0].length][array.length];
    for (int i = 0; i < transposedArray.length; i++) {
      for (int j = 0; j < transposedArray[i].length; j++) {
        transposedArray[i][j] = array[j][i];
      }
    }
    return transposedArray;
  }

  //Copies all elements of an array to a new array.
  public static double[][] copyOf(double[][] original) {
    double[][] copy = new double[original.length][original[0].length];
    for (int i = 0; i < original.length; i++) {
      //Need to do it this way, otherwise it copies only memory location
      System.arraycopy(original[i], 0, copy[i], 0, original[i].length);
    }

    return copy;
  }

  //**********************************//
  //METHODS OF THE HUNGARIAN ALGORITHM//
  //**********************************//
  public static int[] hgAlgorithm(double[][] array, String sumType) {
    //Create the cost matrix
    double[][] cost = copyOf(array);

    //Then array is weight array. Must change to cost.
    if (sumType.equalsIgnoreCase("max")) {
      double maxWeight = findLargest(cost);

      //Generate cost by subtracting.
      for (int i = 0; i < cost.length; i++) {
        for (int j = 0; j < cost[i].length; j++) {
          cost[i][j] = (maxWeight - cost[i][j]);
        }
      }
    }

    //Find largest cost matrix element (needed for step 6).
    double maxCost = findLargest(cost);

    //The mask array.
    int[][] mask = new int[cost.length][cost[0].length];

    //The row covering vector.
    int[] rowCover = new int[cost.length];

    //The column covering vector.
    int[] colCover = new int[cost[0].length];

    //Position of last zero from Step 4.
    int[] zero_RC = new int[2];
    int step = 1;
    boolean done = false;

    //main execution loop
    while (done == false) {
      switch (step) {
        case 1:
          step = hg_step1(step, cost);
          break;
        case 2:
          step = hg_step2(step, cost, mask, rowCover, colCover);
          break;
        case 3:
          step = hg_step3(step, mask, colCover);
          break;
        case 4:
          step = hg_step4(step, cost, mask, rowCover, colCover, zero_RC);
          break;
        case 5:
          step = hg_step5(step, mask, rowCover, colCover, zero_RC);
          break;
        case 6:
          step = hg_step6(step, cost, rowCover, colCover, maxCost);
          break;
        case 7:
          done = true;
          break;
      }
    }//end while

    //Create the returned array.
    int[][] assignment = new int[array.length][2];
    for (int i = 0; i < mask.length; i++) {
      for (int j = 0; j < mask[i].length; j++) {
        if (mask[i][j] == 1) {
          assignment[i][0] = i;
          assignment[i][1] = j;
        }
      }
    }

    int[] is = new int[assignment.length];
    for (int i = 0; i < assignment.length; i++) {
      is[i] = assignment[i][1];
    }
    return is;
  }

  public static int hg_step1(int step, double[][] cost) {

    //What STEP 1 does:
    //For each row of the cost matrix, find the smallest element
    //and subtract it from from every other element in its row. 

    double minval;

    for (int i = 0; i < cost.length; i++) {
      minval = cost[i][0];
      for (int j = 0; j < cost[i].length; j++) //1st inner loop finds min val in row.
      {
        if (minval > cost[i][j]) {
          minval = cost[i][j];
        }
      }
      for (int j = 0; j < cost[i].length; j++) //2nd inner loop subtracts it.
      {
        cost[i][j] = cost[i][j] - minval;
      }
    }

    step = 2;
    return step;
  }

  public static int hg_step2(int step, double[][] cost, int[][] mask, int[] rowCover, int[] colCover) {

    //What STEP 2 does:
    //Marks uncovered zeros as starred and covers their row and column.

    for (int i = 0; i < cost.length; i++) {
      for (int j = 0; j < cost[i].length; j++) {
        if ((cost[i][j] == 0) && (colCover[j] == 0) && (rowCover[i] == 0)) {
          mask[i][j] = 1;
          colCover[j] = 1;
          rowCover[i] = 1;
        }
      }
    }

    clearCovers(rowCover, colCover);	//Reset cover vectors.

    step = 3;
    return step;
  }

  public static int hg_step3(int step, int[][] mask, int[] colCover) {

    //What STEP 3 does:
    //Cover columns of starred zeros. Check if all columns are covered.

    for (int i = 0; i < mask.length; i++) //Cover columns of starred zeros.
    {
      for (int j = 0; j < mask[i].length; j++) {
        if (mask[i][j] == 1) {
          colCover[j] = 1;
        }
      }
    }

    //Check if all columns are covered.
    int count = 0;
    for (int j = 0; j < colCover.length; j++) {
      count = count + colCover[j];
    }

    //Should be cost.length but ok, because mask has same dimensions.	
    if (count >= mask.length) {
      step = 7;
    } else {
      step = 4;
    }

    return step;
  }

  public static int hg_step4(int step, double[][] cost, int[][] mask, int[] rowCover, int[] colCover, int[] zero_RC) {
    //What STEP 4 does:
    //Find an uncovered zero in cost and prime it (if none go to step 6). Check for star in same row:
    //if yes, cover the row and uncover the star's column. Repeat until no uncovered zeros are left
    //and go to step 6. If not, save location of primed zero and go to step 5.

    int[] row_col = new int[2];	//Holds row and col of uncovered zero.
    boolean done = false;
    while (done == false) {
      row_col = findUncoveredZero(row_col, cost, rowCover, colCover);
      if (row_col[0] == -1) {
        done = true;
        step = 6;
      } else {
        //Prime the found uncovered zero.
        mask[row_col[0]][row_col[1]] = 2;

        boolean starInRow = false;
        for (int j = 0; j < mask[row_col[0]].length; j++) {

          //If there is a star in the same row...
          if (mask[row_col[0]][j] == 1) {
            starInRow = true;

            //remember its column.
            row_col[1] = j;
          }
        }

        if (starInRow == true) {

          //Cover the star's row.
          rowCover[row_col[0]] = 1;

          //Uncover its column.
          colCover[row_col[1]] = 0;
        } else {

          //Save row of primed zero.
          zero_RC[0] = row_col[0];

          //Save column of primed zero.
          zero_RC[1] = row_col[1];
          done = true;
          step = 5;
        }
      }
    }

    return step;
  }

  //Aux 1 for hg_step4.
  public static int[] findUncoveredZero(int[] row_col, double[][] cost, int[] rowCover, int[] colCover) {

    //Just a check value. Not a real index.
    row_col[0] = -1;
    row_col[1] = 0;

    int i = 0;
    boolean done = false;
    while (done == false) {
      int j = 0;
      while (j < cost[i].length) {
        if (cost[i][j] == 0 && rowCover[i] == 0 && colCover[j] == 0) {
          row_col[0] = i;
          row_col[1] = j;
          done = true;
        }
        j = j + 1;
      }//end inner while

      i = i + 1;
      if (i >= cost.length) {
        done = true;
      }
    }//end outer while

    return row_col;
  }

  public static int hg_step5(int step, int[][] mask, int[] rowCover, int[] colCover, int[] zero_RC) {
    //What STEP 5 does:	
    //Construct series of alternating primes and stars. Start with prime from step 4.
    //Take star in the same column. Next take prime in the same row as the star. Finish
    //at a prime with no star in its column. Unstar all stars and star the primes of the
    //series. Erasy any other primes. Reset covers. Go to step 3.

    //Counts rows of the path matrix.
    int count = 0;

    //Path matrix (stores row and col).
    int[][] path = new int[(mask[0].length * mask.length)][2];

    //Row of last prime.
    path[count][0] = zero_RC[0];

    //Column of last prime.
    path[count][1] = zero_RC[1];

    boolean done = false;
    while (done == false) {
      int r = findStarInCol(mask, path[count][1]);
      if (r >= 0) {
        count = count + 1;

        //Row of starred zero.
        path[count][0] = r;

        //Column of starred zero.
        path[count][1] = path[count - 1][1];
      } else {
        done = true;
      }

      if (done == false) {
        int c = findPrimeInRow(mask, path[count][0]);
        count = count + 1;

        //Row of primed zero.
        path[count][0] = path[count - 1][0];

        //Col of primed zero.
        path[count][1] = c;
      }
    }//end while

    convertPath(mask, path, count);
    clearCovers(rowCover, colCover);
    erasePrimes(mask);

    step = 3;
    return step;

  }

  //Aux 1 for hg_step5.
  public static int findStarInCol(int[][] mask, int col) {

    //Again this is a check value.
    int r = -1;
    for (int i = 0; i < mask.length; i++) {
      if (mask[i][col] == 1) {
        r = i;
      }
    }

    return r;
  }

  //Aux 2 for hg_step5.
  public static int findPrimeInRow(int[][] mask, int row) {
    int c = -1;
    for (int j = 0; j < mask[row].length; j++) {
      if (mask[row][j] == 2) {
        c = j;
      }
    }

    return c;
  }

  //Aux 3 for hg_step5.
  public static void convertPath(int[][] mask, int[][] path, int count) {
    for (int i = 0; i <= count; i++) {
      if (mask[(path[i][0])][(path[i][1])] == 1) {
        mask[(path[i][0])][(path[i][1])] = 0;
      } else {
        mask[(path[i][0])][(path[i][1])] = 1;
      }
    }
  }

  //Aux 4 for hg_step5.
  public static void erasePrimes(int[][] mask) {
    for (int i = 0; i < mask.length; i++) {
      for (int j = 0; j < mask[i].length; j++) {
        if (mask[i][j] == 2) {
          mask[i][j] = 0;
        }
      }
    }
  }

  //Aux 5 for hg_step5 (and not only).
  public static void clearCovers(int[] rowCover, int[] colCover) {
    for (int i = 0; i < rowCover.length; i++) {
      rowCover[i] = 0;
    }
    for (int j = 0; j < colCover.length; j++) {
      colCover[j] = 0;
    }
  }

  public static int hg_step6(int step, double[][] cost, int[] rowCover, int[] colCover, double maxCost) {
    //What STEP 6 does:
    //Find smallest uncovered value in cost: a. Add it to every element of covered rows
    //b. Subtract it from every element of uncovered columns. Go to step 4.

    double minval = findSmallest(cost, rowCover, colCover, maxCost);

    for (int i = 0; i < rowCover.length; i++) {
      for (int j = 0; j < colCover.length; j++) {
        if (rowCover[i] == 1) {
          cost[i][j] = cost[i][j] + minval;
        }
        if (colCover[j] == 0) {
          cost[i][j] = cost[i][j] - minval;
        }
      }
    }

    step = 4;
    return step;
  }

  //Aux 1 for hg_step6.
  public static double findSmallest(double[][] cost, int[] rowCover, int[] colCover, double maxCost) {

    //There cannot be a larger cost than this.
    double minval = maxCost;

    //Now find the smallest uncovered value.
    for (int i = 0; i < cost.length; i++) {
      for (int j = 0; j < cost[i].length; j++) {
        if (rowCover[i] == 0 && colCover[j] == 0 && (minval > cost[i][j])) {
          minval = cost[i][j];
        }
      }
    }

    return minval;
  }
}