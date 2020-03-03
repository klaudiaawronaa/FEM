import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.sqrt;

public class Element {

    boolean area1 = false;
    boolean area2 = false;
    boolean area3 = false;
    boolean area4 = false;

    public int[] nodeID = new int[4];
    public static int ID = 1;
    double[][] bigJacobiMatrix;
    double[] bigDetJ;
    double[][] bigInverseJacobiMatrix;
    double[][] bigdNdX;
    double[][] bigdNdY;
    double[][] H;
    double[][] bigC;
    double[] P;
    double[][] H_BC;

    public Element(double nL, double nH, GlobalData globaldata) {
        if (ID % nH == 0) ID++;
        nodeID[0] = ID;
        nodeID[3] = ID + 1;
        nodeID[1] = ID + (int) (nH);
        nodeID[2] = ID + (int) (nH) + 1;

        ID++;
        prepareAreas(nL, nH);
        H_BC = new double[4][4];
        H_BC = prepareHBC(globaldata);
        P = new double[4];
        P = prepareP(globaldata);

    }

    public void prepareAreas(double nL, double nH) {
        for (int i = 0; i < nodeID.length; i++) {
            for (int j = 0; j < nodeID.length; j++) {
                if (nodeID[i] == (1 + nH * j)) area1 = true;
                if (nodeID[i] >= (nH * nL - nH + 1)) area2 = true;
                if (nodeID[i] == (nH * j)) area3 = true;
            }
            if (nodeID[i] <= nH) area4 = true;
        }

    }


    public void setBigJacobiMatrix(Jacobian[] jacobiMatrix) {
        List<Double> listOfValuesJakobians = new ArrayList<>();
        List<Double> listOfValuesReverseJakobians = new ArrayList<>();
        bigJacobiMatrix = new double[jacobiMatrix.length][4];
        bigInverseJacobiMatrix = new double[jacobiMatrix.length][4];
        bigDetJ = new double[jacobiMatrix.length];

        for (int i = 0; i < jacobiMatrix.length; i++) {
            for (int j = 0; j < jacobiMatrix[0].jacobiMatrix.length; j++) {
                for (int g = 0; g < jacobiMatrix[0].jacobiMatrix[0].length; g++) {
                    listOfValuesJakobians.add(jacobiMatrix[i].jacobiMatrix[j][g]);
                    listOfValuesReverseJakobians.add(jacobiMatrix[i].inverseJacobiMatrix[j][g]);

                }
            }
            bigDetJ[i] = jacobiMatrix[0].detJ;
        }

        int index = 0;
        for (int i = 0; i < bigJacobiMatrix.length; i++) {
            for (int j = 0; j < bigJacobiMatrix[0].length; j++) {
                bigJacobiMatrix[i][j] = listOfValuesJakobians.get(index);
                bigInverseJacobiMatrix[i][j] = listOfValuesReverseJakobians.get(index);
                index++;
            }
        }

        bigdNdX = new double[jacobiMatrix.length][4];
        ;
        for (int i = 0; i < bigdNdX.length; i++) {
            for (int j = 0; j < bigdNdX[0].length; j++) {
                bigdNdX[i][j] = jacobiMatrix[i].dNdX[j];
            }
        }

        bigdNdY = new double[jacobiMatrix.length][4];
        ;
        for (int i = 0; i < bigdNdY.length; i++) {
            for (int j = 0; j < bigdNdY[0].length; j++) {
                bigdNdY[i][j] = jacobiMatrix[i].dNdY[j];
            }
        }


        H = new double[4][4];
        double suma = 0;
        for (int i = 0; i < H.length; i++)
            for (int j = 0; j < H[0].length; j++) {
                for (int g = 0; g < jacobiMatrix.length; g++) {
                    suma = suma + jacobiMatrix[g].KsumdetJ[i][j];
                }
                H[i][j] = suma;
                suma = 0;
            }
    }


    public void setBigC(Jacobian[] jacobiMatrix) {
        bigC = new double[4][4];
        for (int i = 0; i < bigC.length; i++) {
            for (int j = 0; j < bigC[0].length; j++) {
                for (int g = 0; g < jacobiMatrix.length; g++) {
                    bigC[i][j] += jacobiMatrix[g].C[i][j];
                }
            }
        }
    }

    public double[][] prepareHBC(GlobalData global_data) {

        double ksi1;
        double ksi2;
        double eta1;
        double eta2;
        double[][] for1pc1 = new double[4][4];
        double[][] for2pc1 = new double[4][4];
        double[][] for1pc2 = new double[4][4];
        double[][] for2pc2 = new double[4][4];
        double[][] for1pc3 = new double[4][4];
        double[][] for2pc3 = new double[4][4];
        double[][] for1pc4 = new double[4][4];
        double[][] for2pc4 = new double[4][4];
        double[][] sum1 = new double[4][4];
        double[][] sum2 = new double[4][4];
        double[][] sum3 = new double[4][4];
        double[][] sum4 = new double[4][4];

        double[][] sum = new double[4][4];

        if (area1) {
            ksi1 = ((-1) / sqrt(3));
            eta1 = -1;
            for1pc1 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = (1 / sqrt(3));
            eta2 = -1;
            for2pc1 = fillTable(ksi2, eta2, global_data.alfa);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum1[i][j] = global_data.L / (global_data.nL - 1) / 2 * (for1pc1[i][j] + for2pc1[i][j]);
                }
            }
        }


        if (area2) {
            ksi1 = 1;
            eta1 = (-1) / sqrt(3);
            for1pc2 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = 1;
            eta2 = 1 / sqrt(3);
            for2pc2 = fillTable(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum2[i][j] = global_data.H / (global_data.nH - 1) / 2 * (for1pc2[i][j] + for2pc2[i][j]);
                }
            }
        }


        if (area3) {
            ksi1 = 1 / sqrt(3);
            eta1 = 1;
            for1pc3 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = (-1) / sqrt(3);
            eta2 = 1;
            for2pc3 = fillTable(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum3[i][j] = global_data.L / (global_data.nL - 1) / 2 * (for1pc3[i][j] + for2pc3[i][j]);
                }
            }
        }
        if (area4) {
            ksi1 = -1;
            eta1 = 1 / sqrt(3);
            for1pc4 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = -1;
            eta2 = -1 / sqrt(3);
            for2pc4 = fillTable(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum4[i][j] = global_data.H / (global_data.nH - 1) / 2 * (for1pc4[i][j] + for2pc4[i][j]);
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                sum[i][j] = sum1[i][j] + sum2[i][j] + sum3[i][j] + sum4[i][j];
            }
        }
        return sum;
    }


    public double[][] fillTable(double ksi, double eta, double alfa) {
        double[] helpShape = new double[4];
        double[][] filled = new double[4][4];
        helpShape[0] = 0.25 * (1 - ksi) * (1 - eta);
        helpShape[1] = 0.25 * (1 + ksi) * (1 - eta);
        helpShape[2] = 0.25 * (1 + ksi) * (1 + eta);
        helpShape[3] = 0.25 * (1 - ksi) * (1 + eta);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                filled[i][j] = helpShape[i] * helpShape[j] * alfa;
            }
        }
        return filled;
    }

    public double[] prepareP(GlobalData global_data) {

        double ksi1;
        double ksi2;
        double eta1;
        double eta2;

        double[] P1_1 = new double[4];
        double[] P2_1 = new double[4];
        double[] P3_1 = new double[4];
        double[] P4_1 = new double[4];
        double[] P1_2 = new double[4];
        double[] P2_2 = new double[4];
        double[] P3_2 = new double[4];
        double[] P4_2 = new double[4];
        double[] sum1 = new double[4];
        double[] sum2 = new double[4];
        double[] sum3 = new double[4];
        double[] sum4 = new double[4];

        double[] P = new double[4];

        if (area1) {
            ksi1 = ((-1) / sqrt(3));
            eta1 = -1;
            P1_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = (1 / sqrt(3));
            eta2 = -1;
            P1_2 = fillP(ksi2, eta2, global_data.alfa);

            for (int i = 0; i < 4; i++) {
                sum1[i] = global_data.L / (global_data.nL - 1) / 2 * (P1_1[i] + P1_2[i]);
            }
        }


        if (area2) {
            ksi1 = 1;
            eta1 = (-1) / sqrt(3);
            P2_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = 1;
            eta2 = 1 / sqrt(3);
            P2_2 = fillP(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                sum2[i] = global_data.L / (global_data.nL - 1) / 2 * (P2_1[i] + P2_2[i]);

            }
        }


        if (area3) {
            ksi1 = 1 / sqrt(3);
            eta1 = 1;
            P3_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = (-1) / sqrt(3);
            eta2 = 1;
            P3_2 = fillP(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                sum3[i] = global_data.L / (global_data.nL - 1) / 2 * (P3_1[i] + P3_2[i]);

            }
        }

        if (area4) {
            ksi1 = -1;
            eta1 = 1 / sqrt(3);
            P4_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = -1;
            eta2 = -1 / sqrt(3);
            P4_2 = fillP(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                sum4[i] = global_data.L / (global_data.nL - 1) / 2 * (P4_1[i] + P4_2[i]);

            }
        }

        for (int i = 0; i < 4; i++) P[i] = (sum1[i] + sum2[i] + sum3[i] + sum4[i]) * global_data.ambTemp;
        return P;
    }


    public double[] fillP(double ksi, double eta, double alfa) {
        double[] helpShape = new double[4];
        double[] filled = new double[4];
        helpShape[0] = 0.25 * (1 - ksi) * (1 - eta);
        helpShape[1] = 0.25 * (1 + ksi) * (1 - eta);
        helpShape[2] = 0.25 * (1 + ksi) * (1 + eta);
        helpShape[3] = 0.25 * (1 - ksi) * (1 + eta);
        for (int i = 0; i < 4; i++) {
            filled[i] = helpShape[i] * alfa;
        }
        return filled;
    }

    /*********************** printing *****************************/


    public void printBigJacobiMatrix() {
        System.out.println();
        System.out.println("Big Jacobi Matrix");
        for (int i = 0; i < bigJacobiMatrix.length; i++) {
            for (int j = 0; j < bigJacobiMatrix[0].length; j++) {
                System.out.print(bigJacobiMatrix[i][j] + " ");
                if (j == bigJacobiMatrix[0].length - 1) System.out.println();
            }
        }
    }


    public void printBigReverseJacobiMatrix() {
        System.out.println();
        System.out.println("Big inverse Jacobi Matrix");
        for (int i = 0; i < bigInverseJacobiMatrix.length; i++) {
            for (int j = 0; j < bigInverseJacobiMatrix[0].length; j++) {
                System.out.print(bigInverseJacobiMatrix[i][j] + " ");
                if (j == bigInverseJacobiMatrix[0].length - 1) System.out.println();
            }
        }
    }

    public void printBigDetJ() {
        System.out.println();
        System.out.println("Matrix of determinants");

        for (int i = 0; i < bigDetJ.length; i++)
            System.out.print(bigDetJ[i] + " ");

    }

    public void printBigdNdX() {
        System.out.println();
        System.out.println("Big dNdX");
        for (int i = 0; i < bigdNdX.length; i++) {
            for (int j = 0; j < bigdNdX[0].length; j++) {
                System.out.print(bigdNdX[i][j] + "   ");
                if (j == (bigdNdX[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printBigdNdY() {
        System.out.println();
        System.out.println("Big dNdY");
        for (int i = 0; i < bigdNdY.length; i++) {
            for (int j = 0; j < bigdNdY[0].length; j++) {
                System.out.print(bigdNdY[i][j] + "   ");
                if (j == (bigdNdY[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printH() {
        System.out.println();
        System.out.println("H matrix ");
        for (int i = 0; i < H.length; i++) {
            for (int j = 0; j < H[0].length; j++) {
                System.out.print(H[i][j] + "   ");
                if (j == (H[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printBigC() {
        System.out.println();
        System.out.println("C (big) matrix");
        for (int i = 0; i < bigC.length; i++) {
            for (int j = 0; j < bigC[0].length; j++) {
                System.out.print(bigC[i][j] + "   ");
                if (j == (bigC[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printHBC() {
        System.out.println();
        System.out.println("H_BC matrix");
        for (int i = 0; i < H_BC[0].length; i++) {
            for (int j = 0; j < H_BC.length; j++) {
                System.out.print(H_BC[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printP() {
        System.out.println();
        System.out.println("P vector");

        for (int i = 0; i < P.length; i++)
            System.out.print(P[i] + " ");
        System.out.println();

    }

}
