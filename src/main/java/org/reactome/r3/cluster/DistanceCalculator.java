/*
 * Created on Dec 18, 2009
 *
 */
package org.reactome.r3.cluster;

/**
 * This interface is used to calculate distance between two HierarchicalCluster node.
 * This interface is used in HierarchicalCluster or other similar class.
 * @author wgm
 *
 */
public interface DistanceCalculator {
    
    /**
     * Calculate a distance between two leave nodes represented by ids.
     * @param id1
     * @param id2
     * @return
     */
    public double calculateDistance(String id1, 
                                    String id2);
    
}
