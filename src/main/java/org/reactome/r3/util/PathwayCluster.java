/*
 * Created on Apr 25, 2007
 *
 */
package org.reactome.r3.util;

import java.util.HashSet;
import java.util.Set;

/**
 * This cluster is used to describe a pathway cluster that is grouped from
 * different pathway databases.
 * @author guanming
 *
 */
public class PathwayCluster {
    // index number: actually ids 
    private int index;
    // an double to indicate how good this cluster is
    private double clusterIndex;
    private String label;
    private Set<String> pathways;
    private Set<String> proteins;
    
    /**
     * Defautl constructor.
     *
     */
    public PathwayCluster() {
    }

    public void setClusterIndex(double value) {
        this.clusterIndex = value;
    }
    
    public double getClusterIndex() {
        return this.clusterIndex;
    }
    
    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public String getLabel() {
        return label;
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public Set<String> getPathways() {
        return pathways;
    }

    public void setPathways(Set<String> pathways) {
        this.pathways = pathways;
    }
    
    public void addPathway(String pathway) {
        if (pathways == null)
            pathways = new HashSet<String>();
        pathways.add(pathway);
    }

    public Set<String> getProteins() {
        return proteins;
    }

    public void setProteins(Set<String> proteins) {
        this.proteins = proteins;
    }
    
    public void addProteins(Set<String> set) {
        if (proteins == null)
            proteins = new HashSet<String>();
        proteins.addAll(set);
    }
    
}
