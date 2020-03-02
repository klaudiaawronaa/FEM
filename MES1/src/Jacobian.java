public class Jacobian {

    double dYdEta;
    double dYdKsi;
    double dXdKsi;
    double dXdEta;
    double[][] jacobiMatrix;
    double detJ;
    int numberElement;
    int integrationPoint;
    double[][] inverseJacobiMatrix;
    double dNdX[];
    double dNdY[];
    double dNdXdNdXT[][];
    double dNdYdNdYT[][];
    double KsumdetJ[][];
    double C[][];
    //double H_BC[][];
    //double[] P;

    public Jacobian(Element[] elements, Node[] nodes, elementUni elementuni,
                    globalData globaldata, int numberElement, int integrationPoint) {

        this.integrationPoint = integrationPoint;
        this.numberElement = numberElement;
        dYdEta = 0;
        dXdKsi = 0;
        dXdEta = 0;
        dYdKsi = 0;


        for (int i = 0; i < 4 /*(globaldata.PointQuantity * globaldata.PointQuantity)*/; i++) {
            dYdEta = dYdEta + nodes[elements[numberElement].nodeID[i]].getY() * elementuni.dNdEta[integrationPoint][i];
            dYdKsi = dYdKsi + nodes[elements[numberElement].nodeID[i]].getY() * elementuni.dNdKsi[integrationPoint][i];
            dXdKsi = dXdKsi + nodes[elements[numberElement].nodeID[i]].getX() * elementuni.dNdKsi[integrationPoint][i];
            dXdEta = dXdEta + nodes[elements[numberElement].nodeID[i]].getX() * elementuni.dNdEta[integrationPoint][i];
        }

        jacobiMatrix = new double[2][2];
        jacobiMatrix[0][0] = dXdKsi;
        jacobiMatrix[0][1] = dYdKsi;
        jacobiMatrix[1][0] = dXdEta;
        jacobiMatrix[1][1] = dYdEta;

        detJ = jacobiMatrix[0][0] * jacobiMatrix[1][1] - jacobiMatrix[0][1] * jacobiMatrix[1][0];

        //complement matrix multiplied by 1/det
        inverseJacobiMatrix = new double[2][2];
        inverseJacobiMatrix[0][0] = (1 / detJ) * dYdEta;
        inverseJacobiMatrix[0][1] = (1 / detJ) * (-1) * dYdKsi;
        inverseJacobiMatrix[1][0] = (1 / detJ) * (-1) * dXdEta;
        inverseJacobiMatrix[1][1] = (1 / detJ) * dXdKsi;

        //dNdX for integration point
        dNdX = new double[4];
        for (int i = 0; i < dNdX.length; i++) {
            dNdX[i] = (1 / detJ) * (dYdEta * elementuni.dNdKsi[integrationPoint][i] + (-1) * dXdEta
                    * elementuni.dNdEta[integrationPoint][i]);
        }
        //dNdY for integration point
        dNdY = new double[4];
        for (int i = 0; i < dNdY.length; i++) {
            dNdY[i] = (1 / detJ) * ((-1) * dXdEta * elementuni.dNdKsi[integrationPoint][i] + dXdKsi
                    * elementuni.dNdEta[integrationPoint][i]);
        }
        //dNdX * transposed dNdX for integration point
        dNdXdNdXT = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                dNdXdNdXT[i][j] = dNdX[i] * dNdX[j];
            }
        }
        //dNdY * transposed dNdY for integration point
        dNdYdNdYT = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                dNdYdNdYT[i][j] = dNdY[i] * dNdY[j];
            }
        }

        //K * ( dNdX * transposed dNdX + dNdY * transposed dNdY ) * det J
        KsumdetJ = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                KsumdetJ[i][j] = globaldata.K * (dNdXdNdXT[i][j] * globaldata.gauss[1][0]
                        + dNdYdNdYT[i][j] * globaldata.gauss[1][1]) * detJ;
            }
        }


        //C matrix
        C = new double[4][4];
        for (int i = 0; i < elementuni.shapeFunction[0].length; i++) {
            for (int j = 0; j < elementuni.shapeFunction.length; j++) {

                C[i][j] = detJ * globaldata.c * globaldata.ro
                        * elementuni.shapeFunction[integrationPoint][i] * elementuni.shapeFunction[integrationPoint][j];
            }
        }

    }

 //************************************* printing ************************************//
    public void printJacobiMatrix() {

        System.out.println();
        System.out.println("Element no.: " + numberElement);
        System.out.println("In integration point no.: " + integrationPoint);

        System.out.println("Jacobi matrix");

        for (int i = 0; i < jacobiMatrix.length; i++) {
            for (int j = 0; j < jacobiMatrix[0].length; j++) {
                System.out.print(jacobiMatrix[i][j] + " ");
                if (j == 1) System.out.println();
            }
        }
        System.out.println();
        System.out.println("DET J");
        System.out.println(detJ);
    }

    public void printInverseJacobiMatrix() {

        System.out.println();
        System.out.println("Inverse Jacobi Matrix");

        for (int i = 0; i < inverseJacobiMatrix.length; i++) {
            for (int j = 0; j < inverseJacobiMatrix[0].length; j++) {
                System.out.print(inverseJacobiMatrix[i][j] + " ");
                if (j == 1) System.out.println();
            }
        }
    }

    public void printdNdX() {
        System.out.println("dNdX for " + integrationPoint + " integration point ");
        for (int i = 0; i < dNdX.length; i++)
            System.out.print(dNdX[i]);
    }

    public void printdNdY() {
        System.out.println("dNdY for " + integrationPoint + " integration point ");
        for (int i = 0; i < dNdY.length; i++)
            System.out.print(dNdY[i]);

    }

    public void printdNdXdNdXT() {
        System.out.println();
        System.out.println("dNdX * transposed dNdX for" + integrationPoint + " integration point ");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(dNdXdNdXT[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printdNdYdNdYT() {
        System.out.println();
        System.out.println("dNdY * transposed dNdX for" + integrationPoint + " integration point ");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(dNdYdNdYT[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printKsumDetJ() {
        System.out.println();
        System.out.println("K *( dNdX* transposed dNdX + dNdY * transposed dNdY ) * det J" +
                " for integration point no.:  " + integrationPoint);
        for (int i = 0; i < KsumdetJ[0].length; i++) {
            for (int j = 0; j < KsumdetJ.length; j++) {
                System.out.print(KsumdetJ[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printC() {
        System.out.println();
        System.out.println("C matrix for integration point no.:  " + integrationPoint);
        for (int i = 0; i < C[0].length; i++) {
            for (int j = 0; j < C.length; j++) {
                System.out.print(C[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }


    public double[][] getJacobiMatrix() {
        return jacobiMatrix;
    }
    public double[][] getInverseJacobiMatrix() {
        return inverseJacobiMatrix;
    }
    public double getDetJ() {
        return detJ;
    }
    public double[][] getC() {
        return C;
    }

}
