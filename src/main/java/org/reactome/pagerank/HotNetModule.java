/*
 * Created on Mar 19, 2013
 *
 */
package org.reactome.pagerank;

import java.util.Set;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * @author gwu
 *
 */
@XmlRootElement
public class HotNetModule {
    private Set<String> genes;
    private Double pvalue;
    private Double fdr;
    
    public HotNetModule() {
    }
    
    public Set<String> getGenes() {
        return genes;
    }
    public void setGenes(Set<String> genes) {
        this.genes = genes;
    }
    public Double getPvalue() {
        return pvalue;
    }
    public void setPvalue(Double pvalue) {
        this.pvalue = pvalue;
    }
    public Double getFdr() {
        return fdr;
    }
    public void setFdr(Double fdr) {
        this.fdr = fdr;
    }
}
