/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pathwayprototype;

/**
 *
 * @author ruthgrace
 */

import java.util.HashSet;

public class TestedPathway {
    private double p;
    private String name;
    private HashSet<String> genes;
    public TestedPathway(double p, String name, HashSet<String> genes) {
        this.p = p;
        this.name = name;
        this.genes = new HashSet<String>(genes);
    }
    public double getP() {
        return this.p;
    }
    public String getName() {
        return this.name;
    }
    public HashSet<String> getGenes() {
        return this.genes;
    }
}

