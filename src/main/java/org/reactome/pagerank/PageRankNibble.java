/*
 * Created on Feb 14, 2013
 *
 */
package org.reactome.pagerank;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This is an implementation of PageRank-Nibble in the paper published by 
 * Andersen R et al: http://www.math.ucsd.edu/~fan/wp/localpartition.pdf.
 * Note: this algorithm is just too slow to be useful. Also several numbers
 * used (e.g. 48, 225) are not easily to be understood. 
 * @author gwu
 *
 */
public class PageRankNibble {
    private final boolean debug = true;
    
    public PageRankNibble() {
        
    }
    
    public Set<String> runPageRankNibble(Map<String, Set<String>> nodeToNeighbor,
                                         String startNode,
                                         double phi,
                                         int b) {
        int m = getTotalDegre(nodeToNeighbor);
        if (debug)
            System.out.println("Total edges: " + m);
        int B = (int) Math.ceil(Math.log(m) / Math.log(2)); // change base for log2(m).
        // Make sure passed parameters are correct
        if (phi <= 0.0d || phi > 1.0d)
            throw new IllegalArgumentException("phi should be in (0, 1]");
        if (b < 1 || b > B)
            throw new IllegalArgumentException("b should be in [1, " + B + "]");
        double alpha = phi * phi / (225 * Math.log(100.0d * Math.sqrt(m)));
        if (debug)
            System.out.println("Alpha: " + alpha);
        double epsilon = Math.pow(2.0d, -b) / (48.0d * B);
        if (debug)
            System.out.println("Epsilon: " + epsilon);
        ApproximatePageRank apr = new ApproximatePageRank();
        apr.calculateAPR(nodeToNeighbor,
                         startNode,
                         alpha,
                         epsilon);
        if (debug)
            System.out.println("Finished APR calculation!");
        Map<String, Double> aprScores = apr.getAPR();
        List<String> sweptNodes = sweep(nodeToNeighbor, aprScores);
        // Calculate probability change first since it is the same across all
        List<Integer> volumes = calculateVolumes(nodeToNeighbor, sweptNodes);
        double prob1 = calculateProbability((int)Math.pow(2, b), 
                                            volumes, 
                                            aprScores,
                                            sweptNodes,
                                            nodeToNeighbor);
        double prob2 = calculateProbability((int)Math.pow(2, b - 1), 
                                            volumes, 
                                            aprScores,
                                            sweptNodes,
                                            nodeToNeighbor);
        double diff = prob1 - prob2;
        if (debug)
            System.out.println("Probability change: " + diff);
        if (!(diff > 1.0 / (48.0d * B))) { 
            if (debug)
                System.out.println("Probability Change condition cannot be met!");
            return new HashSet<String>();
        }
                
        double minConductance = Double.MAX_VALUE;
        int minIndex = 0;
        double minVolume = Math.pow(2.0d, b - 1);
        double maxVolume = 2.0d / 3.0d * (2.0d * m);
        for (int i = 0; i < sweptNodes.size(); i++) {
            Set<String> set = new HashSet<String>(sweptNodes.subList(0, i));
            double conductance = calculateConductance(nodeToNeighbor, 
                                                      m, 
                                                      set);
            double volume = calculateVolume(nodeToNeighbor, set);
            // Check three conditions
            if (conductance < phi &&
                volume > minVolume && volume < maxVolume) {
                minConductance = conductance;
                minIndex = i;
                break;
            }
//            if (conductance < minConductance) {
//                minConductance = conductance;
//                minIndex = i;
//            }
        }
        Set<String> rtn = new HashSet<String>(sweptNodes.subList(0, minIndex));
        if (debug) {
            System.out.println("Minimum conductance: " + minConductance);
            System.out.println("Volume: " + calculateVolume(nodeToNeighbor, rtn));
        }
        return rtn;
    }
    
    /**
     * Calculate a spread of a distribution
     * @param k
     * @param volumes
     * @return
     */
    private double calculateProbability(int k,
                                        List<Integer> volumes,
                                        Map<String, Double> nodeToAPR,
                                        List<String> sweptNodes,
                                        Map<String, Set<String>> nodeToNeighbor) {
        int index = Collections.binarySearch(volumes, k);
        double prop = 0.0d;
        if (index >= 0) {
            for (int i = 0; i < index + 1; i++) {
                String node = sweptNodes.get(i);
                prop += nodeToAPR.get(node);
            }
        }
        else { // It is between two numbers
            int insertionIndex = -index - 1; // Volume at insertionIndex should be the first one greater than k
            // Get the first part
            for (int i = 0; i < insertionIndex; i++) {
                String node = sweptNodes.get(i);
                prop += nodeToAPR.get(node);
            }
            // Get the second part
            String currentNode = sweptNodes.get(insertionIndex - 1);
            String nextNode = sweptNodes.get(insertionIndex);
            int degree = nodeToNeighbor.get(currentNode).size();
            double nextProp = nodeToAPR.get(nextNode);
            prop += ((k - volumes.get(insertionIndex - 1)) / degree * nextProp);
        }
        return prop;
    }
    
    private List<Integer> calculateVolumes(Map<String, Set<String>> nodeToNeighbor,
                                           List<String> sweptNodes) {
        List<Integer> volumes = new ArrayList<Integer>();
        for (int i = 0; i < sweptNodes.size(); i++) {
            Set<String> set = new HashSet<String>(sweptNodes.subList(0, i));
            Integer volume = calculateVolume(nodeToNeighbor, set);
            volumes.add(volume);
        }
        return volumes;
    }
    
    private int getTotalDegre(Map<String, Set<String>> nodeToNeighbor) {
        int total = 0;
        for (String node : nodeToNeighbor.keySet()) {
            total += nodeToNeighbor.get(node).size();
        }
        return total / 2;
    }
    
    /**
     * This helper method is used to sweep the APR scores.
     * @param nodeToNeighbor
     * @param nodeToAPR
     * @return
     */
    private List<String> sweep(Map<String, Set<String>> nodeToNeighbor,
                               Map<String, Double> nodeToAPR) {
        if (debug)
            System.out.println("Total nodes: " + nodeToNeighbor.size());
        final Map<String, Double> nodeToRatio = new HashMap<String, Double>();
        for (String node : nodeToAPR.keySet()) {
            Double apr = nodeToAPR.get(node);
            if (apr == 0.0d)  // Nothing there!
                continue;
            double ratio = apr / nodeToNeighbor.get(node).size();
            nodeToRatio.put(node, ratio);
        }
        if (debug)
            System.out.println("Support in APR: " + nodeToRatio.size());
        List<String> sortedNode = new ArrayList<String>(nodeToRatio.keySet());
        Collections.sort(sortedNode, new Comparator<String>() {
            public int compare(String node1, String node2) {
                Double ratio1 = nodeToRatio.get(node1);
                Double ratio2 = nodeToRatio.get(node2);
                return ratio2.compareTo(ratio1);
            }
        });
        return sortedNode;
    }
    
    /**
     * This method is used to calculate conductance for a set, which is defined
     * as edge boundary / min(vol(s), 2m - vol(s)). Here m is total_edges. Edge
     * boundary is the edge number between nodes in the set and nodes not in the set.
     * @param nodeToNeighbor
     * @param set
     * @return
     */
    private double calculateConductance(Map<String, Set<String>> nodeToNeighbor,
                                        int totalEdge,
                                        Set<String> set) {
        // Calculate edge boundary
        int boundary = 0;
        for (String node : set) {
            Set<String> neighbor = nodeToNeighbor.get(node);
            for (String node1 : neighbor) {
                if (!set.contains(node1))
                    boundary ++;
            }
        }
        int volume = calculateVolume(nodeToNeighbor, set);
        return (double) boundary / Math.min(volume, 2.0d * totalEdge - volume);
    }
    
    /**
     * This method is used to calculate volumne for a set, which is defined as
     * sum of node degrees in a set.
     * @param nodeToNeighbor
     * @param set
     * @return
     */
    private int calculateVolume(Map<String, Set<String>> nodeToNeighbor,
                                Set<String> set) {
        int total = 0;
        for (String node : set) {
            total += nodeToNeighbor.get(node).size();
        }
        return total;
    }
    
    @Test
    public void testRun() throws Exception {
        String targetGene = "TP53";
        double phi = 0.25d;
        int b = 9;
        System.out.println("2^b: " + Math.pow(2, b));
        FileUtility fu = new FileUtility();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> nodeToNeighbor = InteractionUtilities.generateProteinToPartners(fis);
        Set<String> neighbor = nodeToNeighbor.get(targetGene);
        System.out.println("Degree for starting node: " + neighbor.size());
        long time1 = System.currentTimeMillis();
        Set<String> cluster = runPageRankNibble(nodeToNeighbor, 
                                                targetGene, 
                                                phi,
                                                b);
        long time2 = System.currentTimeMillis();
        System.out.println("Total running time: " + (time2 - time1));
        System.out.println("Returned cluster: " + cluster.size());
        System.out.println(cluster);
    }
    
}
