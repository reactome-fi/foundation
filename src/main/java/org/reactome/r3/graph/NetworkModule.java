/*
 * Created on Apr 14, 2010
 *
 */
package org.reactome.r3.graph;

import java.util.Set;

/**
 * A simple data structure to describe a network module generated from network clustering algorithm.
 * @author wgm
 *
 */
public class NetworkModule {
    private int index;
    private Set<String> ids;
    private double modularity;
    private Double weightedModularity;
    private double FDR;
    private Double nomialPvalue;
    
    public NetworkModule() {
    }

    public Double getWeightedModularity() {
        return weightedModularity;
    }

    public void setWeightedModularity(Double weightedModularity) {
        this.weightedModularity = weightedModularity;
    }
    
    public Double getNomialPvalue() {
        return nomialPvalue;
    }

    public void setNomialPvalue(Double nomialPvalue) {
        this.nomialPvalue = nomialPvalue;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public Set<String> getIds() {
        return ids;
    }

    public void setIds(Set<String> ids) {
        this.ids = ids;
    }

    public double getModularity() {
        return modularity;
    }

    public void setModularity(double modularity) {
        this.modularity = modularity;
    }

    public double getFDR() {
        return FDR;
    }

    public void setFDR(double fDR) {
        FDR = fDR;
    }
    
    public int getSize() {
        return ids == null ? 0 : ids.size();
    }
    
}
