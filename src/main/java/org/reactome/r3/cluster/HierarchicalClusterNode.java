/*
 * Created on Jun 8, 2009
 *
 */
package org.reactome.r3.cluster;

import java.util.Set;

/**
 * This class is refactored from org.reactome.cancer package. Some of member variables
 * are used as public because of many methods in the original package refer them directly.
 * This should be changed soon.
 * @author wgm
 *
 */
public class HierarchicalClusterNode {
    public Set<String> ids;
    public double pathDistance;
    protected HierarchicalClusterNode childNode1;
    protected HierarchicalClusterNode childNode2;
    // For drawing dendrograms
    int x;
    int y;
    
    public HierarchicalClusterNode() {
    }
    
    public Set<String> getIds() {
        return ids;
    }
    
    public HierarchicalClusterNode getChildNode1() {
        return this.childNode1;
    }
    
    public void setChildNode1(HierarchicalClusterNode node) {
        this.childNode1 = node;
    }
    
    public HierarchicalClusterNode getChildNode2() {
        return this.childNode2;
    }
    
    public void setChildNode2(HierarchicalClusterNode node) {
        this.childNode2 = node;
    }
    
}