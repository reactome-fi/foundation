/*
 * Created on Apr 19, 2007
 *
 */
package org.reactome.r3.cluster;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Iterator;
import java.util.Map;

import org.junit.Test;

/**
 * This class is used to do work related to SOM initialization, training and labelling.
 * The algorithm implemented here is batch-mode.
 * @author guanming
 */
public class SOM {
    // A boolean flag used to switch is FI or proteins should be used
    private static final boolean USE_FI_DATA = true;
    // Make these number small enough so that pathways can be clustered.
    // If the number is too big, the best clustering is one node for each
    // pathway.
    // Right now: 519 pathways in total (April 20, 2007). Also using even number.
    protected static int X_LENGTH = 24;
    protected static int Y_LENGTH = 24;
    // Total learning step
    protected int TOTAL_STEP = 100;
    // Initial neighbor width
    protected int INIT_WIDTH = X_LENGTH / 3;
    // Using INIT_WIDTH = 0 to use SOM as k-means clustering
    //private static final int INIT_WIDTH = 0;
    // Actually data should be used for training
    //private Map<String, Set<String>> pathwayToData;
    protected Map<String, double[]> pathwayToData;
    // A list of Reference nodes
    protected SOMWeightNodes referenceNodes;
    
    public SOM() {
        init();
    }
    
    protected void init() {
        referenceNodes = new SOMWeightNodes();
        referenceNodes.setUpNodes(X_LENGTH, Y_LENGTH);
        try {
            pathwayToData = loadData();
            referenceNodes.initNodeReferenceVector(pathwayToData);
        }
        catch(IOException e) {
            System.err.println("Cannot load data!");
        }
    }
    
    protected Map<String, double[]> loadData() throws IOException {
        PathwayMatrix matrix = new PathwayMatrix();
        matrix.loadData();
        return matrix.getTopicToVector();
    }
    
    /**
     * A batch map learning is implemented in this method (See section 3.6 in the book
     * Self-Organizing Maps).
     */
    public void learn() {
        int width = 0;
        for (int step = 0; step < TOTAL_STEP; step++) {
            width = getNeighborWidth(step);
            System.out.println("Step: " + step + ", width: " + width);
            // Assign the data points to the reference data
            distributeDataPoints();
            averageReferenceNodes(width);
            //referenceNodes.checkNodes();
        }
    }
    
    private void distributeDataPoints() {
        for (Iterator<String> it = pathwayToData.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            double[] vector = pathwayToData.get(topic);
            // Find the best fit reference node
            ReferenceNode winningNode = referenceNodes.searchMatchNode(vector);
            winningNode.addInputVector(vector);
            winningNode.addLabel(topic);
        }
    }
    
    private void averageReferenceNodes(int width) {
        referenceNodes.recalculateNodesVectors(width);
    }
    
    private int getNeighborWidth(int step) {
        return (int) (INIT_WIDTH * (1 - (double) step / TOTAL_STEP));
    }
    
    /**
     * Label the reference nodes in the SOM with the pathway names.
     */
    public void label() {
        System.out.println("Starting labelling...");
        //referenceNodes.checkNodes();
        for (Iterator<String> it = pathwayToData.keySet().iterator(); it.hasNext();) {
            String pathway = it.next();
            double[] data = pathwayToData.get(pathway);
            ReferenceNode winningNode = referenceNodes.searchMatchNode(data);
            // Add input data as labelling
            winningNode.addInputVector(data);
            double distanceSqr = winningNode.calculateDistanceSqrForLearning(data);
            System.out.println("Distance (sqr) for " + pathway + ": " + distanceSqr);
            winningNode.addLabel(pathway);
        }
    }
    
    /**
     * Output the learned SOM
     * @param output
     */
    public void output(OutputStream output) throws IOException {
        referenceNodes.output(output);
    }
    
    @Test
    public void checkDataPoints() throws Exception {
        SOM som = new SOM();
        Map<String, double[]> topicToVector = som.loadData();
        int lessTen = 0;
        for (Iterator<String> it = topicToVector.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            double[] vector = topicToVector.get(topic);
            int c = 0;
            for (double t : vector) {
                if (t == 0.0d)
                    c ++;
            }
            int proteinNumber = vector.length - c;
            if (proteinNumber < 10)
                lessTen ++;
        }
        System.out.println("Total topics: " + topicToVector.size());
        System.out.println("Less Than 10: " + lessTen);
    }
    
    @Test
    public void testSOM() throws Exception {
        SOM som = new SOM();
        som.learn();
        som.label();
        System.out.println("SOM results:");
        som.output(System.out);
    }
    
    public static void main(String[] args) {
        SOM som = new SOM();
        som.learn();
        som.label();
        try {
            som.output(System.out);
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }
    
}
