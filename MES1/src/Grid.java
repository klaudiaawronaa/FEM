public class Grid {
    int tmp_nH;
    int tmp_nL;
    int ne = (tmp_nL - 1) * (tmp_nH - 1);
    int nh = tmp_nL * tmp_nH;
    Element[] elementGrid;
    Node[] nodeGrid;

    public Grid(int nL, int nH, Element[] element, Node[] node) {
        tmp_nH = nH;
        tmp_nL = nL;
        elementGrid = element;
        nodeGrid = node;
    }
}
