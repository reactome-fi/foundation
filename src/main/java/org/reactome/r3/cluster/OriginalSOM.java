/*
 * Created on Jul 17, 2007
 *
 */
package org.reactome.r3.cluster;

import java.util.Iterator;
import java.util.List;

/**
 * This class implements the original incremental SOM algorithm: each input will
 * trigger a learning and adaption process.
 * @author guanming
 *
 */
public class OriginalSOM extends SOM {
    private double learningRate;
    private double learningRateDecay;
    private double updateRadius;
    private double updateRadiusDecay;
    
    public OriginalSOM() {
        init();
        // Initial values
        learningRate = 0.2; // This should be 0.1. Don't change it!
        updateRadius = INIT_WIDTH;
        learningRateDecay = 0.99;
        updateRadiusDecay = 0.999;
    }
    
    /**
     * The original incremental SOM algorithm is implemented in this method. There are
     * two steps in this learning process: find the best match node, and best nodes
     * and its neghbor nodes adaptes after this matching.
     */
    public void learn() {
//        List<String> topics = new ArrayList<String>(pathwayToData.keySet());
//        int size = topics.size();
//        TOTAL_STEP = 100000;
        TOTAL_STEP = 100;
        for (int step = 0; step < TOTAL_STEP; step ++) {
//            int index = (int)(Math.random() * size);
//            String topic = topics.get(index);
//            double[] vector = pathwayToData.get(topic);
//            learn(vector);
            // Print all information for debug
            System.out.println("Step: " + step + 
                               ", width: " + updateRadius + 
                               ", learningRate: " + learningRate);
            // Feeding the pathway data into the reference vectors
            for (Iterator<String> it = pathwayToData.keySet().iterator(); it.hasNext();) {
                String topic = it.next();
                double[] vector = pathwayToData.get(topic);
                learn(vector);
            }
//            learningRate = 0.90 * (1.0 - (double) step / TOTAL_STEP);
//            updateRadius = 1.0 + (INIT_WIDTH - 1.0) * (1.0 - (double) step / TOTAL_STEP);
        }
    }
    
    private void learn(double[] input) {
        ReferenceNode winningNode = referenceNodes.searchMatchNode(input);
        train(winningNode, input);
        // Update learningRate and updateRadius after each learning
        learningRate = 0.1 + ((learningRate - 0.1) * learningRateDecay);
        updateRadius = 1.0 + ((updateRadius - 1.0) * updateRadiusDecay);
    }
    
    private void train(ReferenceNode node, double[] input) {
        int ur = (int) Math.round(updateRadius);
        List<ReferenceNode> updateNodes = referenceNodes.getNodesWithin(node.getLocation(),
                                                                        ur);
        double dist = 0.0;
        for (ReferenceNode tNode : updateNodes) {
            // tNode for temp node
            // compare distance
            dist = node.calculateDistance(tNode);
            if (dist < updateRadius) {
                // Update reference vectors
                double[] refVector = tNode.getReferenceVector();
                for (int i = 0; i < refVector.length; i++) {
                    refVector[i] += learningRate * (input[i] - refVector[i]);
                }
            }
        }
    }
}
