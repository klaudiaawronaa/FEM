public class Node {
    double x;
    double y;
    //double t;
    int id;


    public Node (double x, double y, int id)
    {   this.id=id;
        this.x=x;
        this.y=y;
    }

    public double getX(){
        return x;
    }

    public double getY(){
        return y;
    }

}
