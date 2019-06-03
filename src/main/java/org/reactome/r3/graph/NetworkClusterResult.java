/*
 * Created on Jul 23, 2010
 *
 */
package org.reactome.r3.graph;

import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;


@XmlRootElement
public class NetworkClusterResult {
    // class used for network clustering
    private String clsName;
    private Double modularity;
    private List<GeneClusterPair> geneClusterPairs;
    
    public NetworkClusterResult() {
    }

    public String getClsName() {
        return clsName;
    }

    public void setClsName(String clsName) {
        this.clsName = clsName;
    }

    public Double getModularity() {
        return modularity;
    }

    public void setModularity(Double modularity) {
        this.modularity = modularity;
    }

    public List<GeneClusterPair> getGeneClusterPairs() {
        return geneClusterPairs;
    }

    public void setGeneClusterPairs(List<GeneClusterPair> geneClusterPairs) {
        this.geneClusterPairs = geneClusterPairs;
    }
    
    
}
