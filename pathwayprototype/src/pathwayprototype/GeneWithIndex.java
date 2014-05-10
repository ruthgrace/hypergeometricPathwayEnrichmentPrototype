/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pathwayprototype;

/**
 *
 * @author ruthgrace
 */



public class GeneWithIndex {
    private int index;
    private String name;
    public GeneWithIndex(int index, String name) {
        this.index = index;
        this.name = name;
    }
    public double getIndex() {
        return this.index;
    }
    public String getName() {
        return this.name;
    }
}

