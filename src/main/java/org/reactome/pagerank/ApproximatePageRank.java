/*
 * Created on Feb 6, 2013
 *
 */
package org.reactome.pagerank;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;
import org.junit.Test;
import org.reactome.r3.graph.JungGraphUtilities;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

import edu.uci.ics.jung.algorithms.scoring.PageRankWithPriors;
import edu.uci.ics.jung.graph.Graph;

/**
 * This class is used to calculate approximate page rank based on paper by 
 * Andersen R et al: http://www.math.ucsd.edu/~fan/wp/localpartition.pdf.
 * @author gwu
 *
 */
public class ApproximatePageRank {
    private String startNode;
    private double alpha;
    private double epsilon;
    // Hold the actual approximate page ranks
    private Map<String, Double> nodeToRank;
    private Map<String, Double> nodeToResidue;
    
    public ApproximatePageRank() {
    }
    
    /**
     * The actual method that should be called in order to calculate approximate
     * page rank (APR). According to the original paper, the running time should
     * be O(1 / (alapha * epsilon)).
     * @param startNode
     * @param alpha
     * @param epsilon
     */
    public void calculateAPR(Map<String, Set<String>> nodeToNeighbor,
                             String startNode,
                             double alpha,
                             double epsilon) {
        this.alpha = alpha;
        this.epsilon = epsilon;
        this.startNode = startNode;
        initialize();
        List<String> pushQueue = new ArrayList<String>();
        // First node should be the start node since it has the highest residue/degree.
        // Others are zero for this value.
        if (isPushNode(startNode, nodeToNeighbor))
            pushQueue.add(startNode);
        while (pushQueue.size() > 0) {
            String pushNode = pushQueue.get(0);
            push(nodeToNeighbor, pushNode);
            if (!isPushNode(pushNode, nodeToNeighbor)) {
                pushQueue.remove(0); // This node has been handled correctly.
            }
            // Check pushNode's neighbor since their residues have changed.
            for (String node1 : nodeToNeighbor.get(pushNode)) {
                if (isPushNode(node1, nodeToNeighbor) && !pushQueue.contains(node1)) {
                    pushQueue.add(node1);
                }
            }
        }
        double max = Double.MIN_VALUE;
        for (String gene : nodeToNeighbor.keySet()) {
            Set<String> neighbor = nodeToNeighbor.get(gene);
            Double residue = nodeToResidue.get(gene);
            if (residue == null)
                continue;
            double ratio = residue / neighbor.size();
            if (ratio > max)
                max = ratio;
        }
        System.out.println("Maximum residue / degree: " + max);
    }
    
    /**
     * Check if a node has residue/degree >= episilon. If true,
     * it should be a push node.
     * @param node
     * @return
     */
    private boolean isPushNode(String node,
                               Map<String, Set<String>> nodeToNeighbor) {
        double ratio = nodeToResidue.get(node) / nodeToNeighbor.get(node).size();
        return ratio >= epsilon;
    }
    
    private void initialize() {
        if (nodeToRank == null)
            nodeToRank = new HashMap<String, Double>();
        else
            nodeToRank.clear();
        if (nodeToResidue == null)
            nodeToResidue = new HashMap<String, Double>();
        else
            nodeToResidue.clear();
        nodeToResidue.put(startNode, 1.0d);
    }
    
    /**
     * Push method as described in the original paper.
     * @param nodeToNeighbor
     * @param node
     */
    private void push(Map<String, Set<String>> nodeToNeighbor,
                      String node) {
        Double rank = nodeToRank.get(node);
        if (rank == null) {
            rank = 0.0d;
        }
        Double residue = nodeToResidue.get(node);
        if (residue == null)
            residue = 0.0d;
        nodeToRank.put(node, rank + alpha * residue);
        nodeToResidue.put(node, (1.0d - alpha) * residue / 2.0d);
        Set<String> neighbor = nodeToNeighbor.get(node);
        for (String node1 : neighbor) {
            Double residue1 = nodeToResidue.get(node1);
            if (residue1 == null)
                residue1 = 0.0d;
            nodeToResidue.put(node1, residue1 + (1.0d - alpha) * residue / (2 * neighbor.size()));
        }
    }
    
    public Map<String, Double> getAPR() {
        return this.nodeToRank;
    }
    
    public Map<String, Double> getResidueRank() {
        return this.nodeToResidue;
    }

    public String getStartNode() {
        return startNode;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getEpsilon() {
        return epsilon;
    }

    @Test
    public void testRun() throws IOException {
        String targetGene = "TP53";
        double alpha = 0.15d;
        FileUtility fu = new FileUtility();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> nodeToNeighbor = InteractionUtilities.generateProteinToPartners(fis);
        long time1 = System.currentTimeMillis();
        calculateAPR(nodeToNeighbor, 
                     targetGene,
                     alpha,
                     1.1e-6d);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for running APR: " + (time2 - time1));
//        if (true)
//            return;
        List<String> geneList = new ArrayList<String>(nodeToRank.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double rank1 = nodeToRank.get(gene1);
                Double rank2 = nodeToRank.get(gene2);
                return rank2.compareTo(rank1);
            }
        });
        
        // Calculate original page rank using starting from TP53
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        JungGraphUtilities gu = new JungGraphUtilities();
        Graph<String, String> graph = gu.createJungGraph(fis);
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        for (String gene : genes) {
            if (gene.equals(targetGene)) {
                geneToScore.put(gene, 1.0d);
            }
            else
                geneToScore.put(gene, 0.0d);
        }
        Transformer<String, Double> mutSigPriors = TransformerUtils.mapTransformer(geneToScore);
        // Based on the original paper, pagerank used in APR is different from original one in alpha value
        double rAlpha = 2.0d * alpha / (1.0d + alpha);
        System.out.println("Alpha change: " + alpha + " -> " + rAlpha);
        PageRankWithPriors<String, String> priorRanks = new PageRankWithPriors<String, String>(graph, 
                                                                                               mutSigPriors, 
                                                                                               rAlpha);
        long time3 = System.currentTimeMillis();
        priorRanks.evaluate();
        long time4 = System.currentTimeMillis();
        System.out.println("Running original pagerank using JUNG: " + (time4 - time3));
        System.out.println();
        
        System.out.println("Order\tGene\tAPR\tResidue\tSum\tPageRank");
        for (int i = 0; i < 100; i++) {
            String gene = geneList.get(i);
            double sum = nodeToRank.get(gene) + nodeToResidue.get(gene);
            double pageRank = priorRanks.getVertexScore(gene);
            System.out.println((i + 1) + "\t" + gene + "\t" + nodeToRank.get(gene) + "\t" + nodeToResidue.get(gene) + 
                               "\t" + sum + "\t" + pageRank);
        }
    }
}
