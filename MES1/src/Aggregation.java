import static java.util.Arrays.sort;

public class Aggregation {

    double[][] globalH;
    double[][] globalH_BC;
    double[][] globalC;
    double[] globalP;
    double[][] A;
    double[] B;
    double[] tempVector;
    double[] tmpForTemp;
    double[] Bpart;

    public Aggregation(globalData global_data, Element[] elements) {
        globalH = new double[(int) (global_data.nL * global_data.nL)][(int) (global_data.nH * global_data.nH)];
        globalH_BC = new double[(int) (global_data.nL * global_data.nL)][(int) (global_data.nH * global_data.nH)];
        globalC = new double[(int) (global_data.nL * global_data.nL)][(int) (global_data.nH * global_data.nH)];
        globalP = new double[(int) (global_data.nL * global_data.nL)];
        tmpForTemp = new double[(int) (global_data.nL * global_data.nL)];

        int k = 0;
        while (k < (((global_data.nL - 1) * (global_data.nH - 1)))) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    globalH[elements[k].nodeID[i] - 1][elements[k].nodeID[j] - 1] += elements[k].H[i][j] + elements[k].H_BC[i][j];
                    globalC[elements[k].nodeID[i] - 1][elements[k].nodeID[j] - 1] += elements[k].bigC[i][j];
                }
                globalP[elements[k].nodeID[i] - 1] += elements[k].P[i];
            }
            k++;
        }
        //printGlobalH();
        //printGlobalP();
        compute(global_data);
    }


    public void compute(globalData global_data){

        A = new double[(int) (global_data.nL*global_data.nL)][(int)(global_data.nH*global_data.nL)];
        B = new double[(int)(global_data.nH*global_data.nL)];
        Bpart = new double[(int)(global_data.nH*global_data.nL)];
        jacobiMethod jacobi;

        tempVector = new double[(int)(global_data.nH*global_data.nL)];
        for (int i = 0; i < tempVector.length; i++) {
            tempVector[i] = global_data.initTemp; }


        System.out.println();
        System.out.println();
        System.out.println("Max and min temperature in further time steps");
        System.out.println("Time[s] \t MinTemp[s] \t MaxTemp[s]");


        for (int i = 0; i < A[0].length; i++) {
            for (int j = 0; j < A.length; j++) {
                A[i][j] = globalH[i][j] + globalC[i][j] / global_data.dT;
            }
        }
       //printA();
        /***************************** MAIN LOOP ******************************************/

        double time = global_data.dT;
        while (time <= global_data.simTime) {

            for (int i = 0; i < Bpart.length; i++) {
                Bpart[i] = 0.0;
            }

            for (int i = 0; i < Bpart.length; i++) {
                for (int j = 0; j < Bpart.length; j++) {
                    Bpart[i] += globalC[i][j] / global_data.dT * tempVector[j];
                }
            }

            for (int i=0; i<B.length; i++){
                B[i] = globalP[i] + Bpart[i];
            }

            jacobi= new jacobiMethod(global_data, A, B);
            for (int j = 0; j < jacobi.getX1().length; j++) {
                tempVector[j] = jacobi.getX1()[j];
                tmpForTemp[j] = jacobi.getX1()[j];
            }

            sort(tmpForTemp);
            System.out.println(time + "\t" + tmpForTemp[0] + "\t" + tmpForTemp[tmpForTemp.length-1]);
            time += global_data.dT;
        }
    }



    public void printGlobalH() {
        System.out.println();
        System.out.println("Global H");
        for (int i = 0; i < globalH[0].length; i++) {
            for (int j = 0; j < globalH.length; j++) {
                {
                    System.out.print(globalH[i][j] + " ");
                    if (j == globalH[0].length - 1) System.out.println();
                }
            }
        }
    }

    public void printGlobalC() {
        System.out.println();
        System.out.println("Global C");
        for (int i = 0; i < globalC[0].length; i++) {
            for (int j = 0; j < globalC.length; j++) {
                {
                    System.out.print(globalC[i][j] + " ");
                    if (j == globalC[0].length - 1) System.out.println();

                }
            }
        }
    }

    public void printGlobalP() {
        System.out.println();
        System.out.println("Global P");
        for (int i = 0; i < globalP.length; i++) {
            System.out.print(globalP[i] + " ");
        }
    }


    public void printA() {
        System.out.println();
        System.out.println("A");
        for (int i = 0; i < A[0].length; i++) {
            for (int j = 0; j < A.length; j++) {
                {
                    System.out.print(A[i][j] + " ");
                    if (j == A[0].length - 1) System.out.println();
                }
            }
        }
    }

    public void printB() {
        System.out.println();
        System.out.println("B");
        for (int i = 0; i < B.length; i++) {
            System.out.print(B[i] + " ");
        }
        System.out.println();
    }

}

