/*
 * Created on Dec 18, 2009
 *
 */
package org.reactome.r3.graph;

import java.util.Map;
import java.util.Set;

import org.reactome.r3.cluster.DistanceCalculator;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to calculate topological overlap as defined in Ravasz et al (2002)
 * Science 297: 1551-1555.
 * @author wgm
 *
 */
public class TopologicalOverlapCalculator implements DistanceCalculator {
    private Map<String, Set<String>> idToPartners;
    
    public TopologicalOverlapCalculator() {
    }
    
    public void setIdToPartners(Map<String, Set<String>> idToPartners) {
        this.idToPartners = idToPartners;
    }

    /**
     * Calculate distance between two ids based on topological overlap.
     */
    public double calculateDistance(String id1, String id2) {
        Set<String> partners1 = idToPartners.get(id1);
        Set<String> partners2 = idToPartners.get(id2);
        int minDegree = Math.min(partners1.size(), partners2.size());
        Set<String> shared = InteractionUtilities.getShared(partners1, partners2);
        int sharedSize = shared.size();
        // Check if id1 and id2 interact. If true, add one to total shared
        if (partners1.contains(id2)) 
            sharedSize ++;
        // We want to calculate distance. For toplogical overlap,
        // the smaller, the better. The max value is 1.0.
        return 1.0 - (double) sharedSize / minDegree;
    }
    
}
