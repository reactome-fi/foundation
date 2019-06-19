/*
 * Created on Apr 23, 2013
 *
 */
package org.reactome.pagerank;

import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * This model class is used to represent the results from one HotNet run.
 * @author gwu
 *
 */
@XmlRootElement
public class HotNetResult {
    private Double delta;
    private Double fdrThreshold;
    private Integer permutation;
    private Boolean useAutoDelta;
    private List<HotNetModule> modules;
    
    public HotNetResult() {
    }

    public Double getDelta() {
        return delta;
    }

    public void setDelta(Double delta) {
        this.delta = delta;
    }

    public List<HotNetModule> getModules() {
        return modules;
    }

    public void setModules(List<HotNetModule> modules) {
        this.modules = modules;
    }

    public Double getFdrThreshold() {
        return fdrThreshold;
    }

    public void setFdrThreshold(Double fdrThreshold) {
        this.fdrThreshold = fdrThreshold;
    }

    public Integer getPermutation() {
        return permutation;
    }

    public void setPermutation(Integer permutation) {
        this.permutation = permutation;
    }

    public Boolean getUseAutoDelta() {
        return useAutoDelta;
    }

    public void setUseAutoDelta(Boolean useAutoDelta) {
        this.useAutoDelta = useAutoDelta;
    }
    
}
