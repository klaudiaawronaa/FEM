import java.util.Scanner;
import java.nio.file.Paths;
import java.nio.file.Path;
import java.io.IOException;

public class Main {

    public static void main(String[] args) {

        /************************ READING FILE *************************/
        String[] splitedArray;
        double[] tmp = new double[13];
        int i = 0;
        try {
            Path path = Paths.get("C:/Users/wrona/OneDrive/Pulpit/data.txt");

            Scanner scanner = new Scanner(path);
            System.out.println("Read text file using Scanner");
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                splitedArray = line.split("=");
                tmp[i] = Double.parseDouble(splitedArray[1]);

                i++;
                System.out.println(line);
            }

            scanner.close();
        } catch (IOException e) {
        }
        globalData dataByUser = new globalData(tmp);

        /******************************** ELEMENTS MATRIX ************************************/
        double iloczyn_double = (dataByUser.nL - 1) * (dataByUser.nH - 1);
        int iloczyn = (int) (iloczyn_double);
        Element[] elements = new Element[iloczyn];
        //Element[] elements = new Element[(dataByUser.nL-1)*(dataByUser.nH-1)];
        for (i = 0; i < ((dataByUser.nL - 1) * (dataByUser.nH - 1)); i++) {
            elements[i] = new Element(dataByUser.nL, dataByUser.nH, dataByUser);

            System.out.println("Nodes - element no. " + i);
            System.out.println(elements[i].nodeID[0]);
            System.out.println(elements[i].nodeID[1]);
            System.out.println(elements[i].nodeID[2]);
            System.out.println(elements[i].nodeID[3]);

            System.out.println("**************************************");
        }

        /**************************** NODES MATRIX ****************************************/

        double x = 0;
        double y = 0;
        int id = 1;

        Node[] nodes = new Node[(int) (dataByUser.nL * dataByUser.nH)];
        for (i = 0; i < dataByUser.nL * dataByUser.nH; i++) {
            nodes[i] = new Node(x, y, id);
            if (y < dataByUser.H) {
                y = y + (dataByUser.H / (dataByUser.nH - 1));
            } else {
                y = 0;
                x = x + (dataByUser.L / (dataByUser.nL - 1));
            }

            id++;
        }

        /************************************* CHECKING NODE 0 ******************************/

        System.out.println("");
        System.out.println("WSPOLRZEDNE NODE 0");
        System.out.println(nodes[0].id);
        System.out.println("X: " + nodes[0].x);
        System.out.println("Y: " + nodes[0].y);


        /********************** GAUSS *****************************/

        System.out.println("");
        System.out.println("Gauss table (points + weights) ");
        for (i = 0; i < dataByUser.gauss.length; i++) {
            for (int j = 0; j < dataByUser.gauss[0].length; j++) {
                System.out.print(dataByUser.gauss[i][j] + " ");
                if (j == dataByUser.gauss[0].length - 1) System.out.println("");
            }
        }

        /************************* COORDINATES TABLE  ***************************/

        System.out.println("");
        System.out.println("Coordinates table: ");

        for (i = 0; i < dataByUser.points.length; i++) {
            for (int j = 0; j < dataByUser.points[0].length; j++) {
                System.out.print(dataByUser.points[i][j] + " ");
                if (j == dataByUser.points[0].length - 1) System.out.println("");
            }

        }

        /************************* Universal element *********************************/
        elementUni elementuni = new elementUni(dataByUser);
        elementuni.printShapeFunction();
        elementuni.printdNdKSI();
        elementuni.printdNdEta();



        /***************************** Jacobian *******************************************/

        Jacobian[] jakobianElementZero= new Jacobian[dataByUser.PointQuantity*dataByUser.PointQuantity];
        for(i=0; i<jakobianElementZero.length; i++) {
            jakobianElementZero[i] = new Jacobian(elements, nodes,
                    elementuni, dataByUser, 0, i);
            jakobianElementZero[i].printdNdXdNdXT();
            jakobianElementZero[i].printKsumDetJ();
            /*jakobianElementZero[i].printMacierzJakobiego();
            jakobianElementZero[i].printdNdY();*/
        }

        for (i =0; i < elements.length; i++ ){
        elements[i].setBigJacobiMatrix(jakobianElementZero);
        elements[0].printBigJacobiMatrix();
        elements[0].printBigReverseJacobiMatrix();
        elements[0].printBigDetJ();
        elements[0].printBigdNdX();
        elements[0].printBigdNdY();
        elements[0].printH();

        jakobianElementZero[0].printC();

        elements[i].setBigC(jakobianElementZero);
        elements[0].printBigC();
        elements[0].printHBC();
        elements[0].printP();}

        Aggregation agr = new Aggregation(dataByUser, elements);
        //agr.printGlobalC();
        //agr.printGlobalH();
        //elements[8].printP();

       //agr.printGlobalP();

        //elements[0].printHBC();

    }

}

