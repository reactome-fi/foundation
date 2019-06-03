/*
 * Created on Apr 19, 2007
 *
 */
package org.reactome.r3.cluster;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * This class is used to describe the reference node in a SOM.
 * @author guanming
 *
 */
public class ReferenceNode {
    // Position in SOM
    private Point location;
    // An index to indicate the quality of this ReferenceNode as a cluster
    private double index;
    // A collection data of either FI interactions or Protein ids
    private Set<String> referenceData;
    // value for this ReferenceNode
    private double[] referenceVector;
    // Use input data as reference data
    private List<Set<String>> inputReferenceData;
    // A sublist of the input data
    private List<Set<String>> inputData;
    // A sublist of vectors from input data
    private List<double[]> inputVectors;
    // A list of pathway labels
    private Set<String> labels;
    
    public ReferenceNode() {
    }
    
    public List<Set<String>> getInputReferenceData() {
        return inputReferenceData;
    }
    
    public double[] getReferenceVector() {
        return this.referenceVector;
    }
    
    public void addInputReferenceData(Set<String> input) {
        if (inputReferenceData == null)
            inputReferenceData = new ArrayList<Set<String>>();
        inputReferenceData.add(input);
    }
    
    public void setInputReferenceData(List<Set<String>> data) {
        if (inputReferenceData == null)
            inputReferenceData = new ArrayList<Set<String>>();
        else
            inputReferenceData.clear();
        inputReferenceData.addAll(data);
    }
    
    public void clearInputReferenceData() {
        if (inputReferenceData != null)
            inputReferenceData.clear();
    }
    
    public void addLabel(String label) {
        if (labels == null)
            labels = new HashSet<String>();
        labels.add(label);
    }
    
    public Set<String> getLabels() {
        return labels;
    }
    
    public void setLocation(int x, int y) {
        location = new Point(x, y);
    }
    
    public void setLocation(Point p) {
        this.location = p;
    }
    
    public Point getLocation() {
        return this.location;
    }
    
    public void setReferenceVector(double[] vector) {
        this.referenceVector = vector;
    }
    
    public void setReferenceData(Set<String> data) {
        if (referenceData == null)
            referenceData = new HashSet<String>();
        else
            referenceData.clear();
        if (data != null)
            referenceData.addAll(data);
    }
    
    public Set<String> getReferenceData() {
        return referenceData;
    }
    
    public double calculateDistanceSqrForLearning(double[] data) {
        // reference vector may be null. If it is null,
        // the distance should be 0.0 during learning
        // so that this ReferneceNode can be used.
        if (referenceVector == null) {
            return 0.0d;
        }
        return calculateDistanceSqr(data, referenceVector);
    }
    
    private double calculateDistanceSqr(double[] v1,
                                        double[] v2) {
        double total = 0.0d;
        double diff = 0.0d;
        for (int i = 0; i < v1.length; i++) {
             diff = v1[i] - v2[i];
             total += diff * diff;
        }
        return total;
    }
    
    public double calculateDotProduct(double[] data) {
        return calculateDotProduct(data, referenceVector);
    }
    
    private double calculateDotProduct(double[] v1,
                                       double[] v2) {
        double total = 0.0;
        for (int i = 0; i < v1.length; i++) {
            total += v1[i] * v2[i];
        }
        return total;
    }
    
    /**
     * Similarity between two pathways are calculated as the value calculated by the following
     * formula:
     *       sim = o / min (p1, p2)
     * where o is the number of the overlapping data points (FI or IDs), p1 or p2 are the
     * total data points for pathway p1 or p2.
     * @param data
     * @return
     */
    public double calculateSimilarity(Set<String> data) {
        if (inputReferenceData == null || inputReferenceData.size() == 0)
            return 0.0d;
        int size = inputReferenceData.size();
        double total = 0.0d;
        for (Set<String> tmp : inputReferenceData) {
            total += calculateSimilarity(tmp, data);
        }
        return total / size;
        //if (referenceData == null || referenceData.size() == 0)
        //    return 0.0d;
        //return calculateSimilarity(referenceData, data);
    }
    
    /**
     * Calculate an average similiarty between the specified node and this node
     * based on input node.
     * input data.
     * @param data
     * @return
     */
    public double calculateSimilarity(ReferenceNode node) {
        List<Set<String>> inputData1 = node.getInputData();
        if (inputData == null || inputData.size() == 0 ||
            inputData1 == null || inputData1.size() == 0)
            return 0.0d;
        double total = 0.0d;
        int totalCount = 0;
        for (Set<String> input : inputData) {
            for (Set<String> input1 : inputData1) {
                totalCount ++;
                total += calculateSimilarity(input, input1);
            }
        }
        return total / totalCount;
    }
    
    private double calculateSimilarity(Set<String> data1, Set<String> data2) {
        int overlapping = 0;
        int size1 = data1.size();
        int size2 = data2.size();
        int min = Math.min(size1, size2);
        if (min == size1) {
            for (String tmp : data1) {
                if (data2.contains(tmp))
                    overlapping ++;
            }
        }
        else {
            for (String tmp : data2) {
                if (data1.contains(tmp))
                    overlapping ++;
            }
        }
        return (double) overlapping / min;
    }
    
    public double calculateDistance(ReferenceNode node) {
        double sqr = calculateDistanceSqr(referenceVector, 
                                          node.referenceVector);
        return Math.sqrt(sqr);
    }
    
    private double calculateDistance(double[] v1,
                                     double[] v2) {
        double sqr = calculateDistanceSqr(v1, v2);
        return Math.sqrt(sqr);
    }
    
    public double calculateClusterIndex() {
        if (inputData == null || inputData.size() == 0)
            return -1.0d; // Not a cluster
        if (inputData.size() == 1)
            return 1.0d; // Just one topic
        int c = 0;
        double total = 0.0;
        for (int i = 0; i < inputData.size() - 1; i++) {
            Set<String> data1 = inputData.get(i);
            for (int j = i + 1; j < inputData.size(); j++) {
                Set<String> data2 = inputData.get(j);
                total += calculateSimilarity(data1, data2);
                c ++;
            }
        }
        return total / c;
    }
    
    public double calculateAverageDistance() {
        if (inputVectors == null || inputVectors.size() == 0)
            return -1.0d; // Not a cluster
        if (inputVectors.size() == 1)
            return 0.0d; // Just one topic
        int c = 0;
        double total = 0.0;
        for (int i = 0; i < inputVectors.size() - 1; i++) {
            double[] data1 = inputVectors.get(i);
            for (int j = i + 1; j < inputVectors.size(); j++) {
                double[] data2 = inputVectors.get(j);
                total += calculateDistance(data1, data2);
                c ++;
            }
        }
        return total / c;        
    }
    
    public void resetInputData() {
        if (inputData != null)
            inputData.clear();
    }
    
    public void resetInputVectors() {
        if (inputVectors != null)
            inputVectors.clear();
    }
    
    public void resetLabels() {
        if (labels != null)
            labels.clear();
    }
    
    public void addInputData(Set<String> data) {
        if (inputData == null)
            inputData = new ArrayList<Set<String>>();
        inputData.add(data);
    }
    
    public void addInputVector(double[] vector) {
        if (inputVectors == null)
            inputVectors = new ArrayList<double[]>();
        inputVectors.add(vector);
    }
    
    public List<double[]> getInputVectors() {
        return this.inputVectors;
    }
    
    public List<Set<String>> getInputData() {
        return inputData;
    }
    
    public boolean isEmpty() {
        if (inputData == null || inputData.size() == 0)
            return true;
        return false;
    }
    
    public void setIndex(double index) {
        this.index = index;
    }
    
    public double getIndex() {
        return this.index;
    }
    
}
