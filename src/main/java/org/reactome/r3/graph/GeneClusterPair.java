package org.reactome.r3.graph;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * A simple class for XML converting.
 * @author wgm
 */
@XmlRootElement
public class GeneClusterPair {
    private String geneId;
    private Integer cluster;
    
    public GeneClusterPair() {
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public Integer getCluster() {
        return cluster;
    }

    public void setCluster(Integer cluster) {
        this.cluster = cluster;
    }
    
}