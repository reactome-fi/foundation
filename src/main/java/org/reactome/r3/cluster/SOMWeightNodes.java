/*
 * Created on Apr 19, 2007
 *
 */
package org.reactome.r3.cluster;

import java.awt.Point;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

import org.apache.commons.math.random.RandomDataImpl;

/**
 * This is a list of ReferenceNodes so that search-related methods can be grouped
 * in this simple class.
 * @author guanming
 *
 */
public class SOMWeightNodes {
    // Flag to switch checking
    private boolean useDotProduct = false;
    // A list of all ReferenceNodes used in a SOM
    private List<ReferenceNode> nodeList;
    // This map is used to speed up the searching based on position,
    // which is encoded as "x,y"
    private Map<String, ReferenceNode> posToNode;
    
    public SOMWeightNodes() {
        nodeList = new ArrayList<ReferenceNode>();
        posToNode = new HashMap<String, ReferenceNode>();
    }
    
    public void setUpNodes(int xLength, int yLength) {
        for (int i = 0; i < xLength; i++) {
            for (int j = 0; j < yLength; j++) {
                ReferenceNode node = new ReferenceNode();
                node.setLocation(i, j);
                nodeList.add(node);
                posToNode.put(i + "," + j,
                              node);
            }
        }
    }
    
    public List<ReferenceNode> getReferenceNodes() {
        return this.nodeList;
    }
    
    public void addNode(ReferenceNode node) {
        nodeList.add(node);
        Point location = node.getLocation();
        posToNode.put(location.x + "," + location.y,
                      node);
    }
    
    public void initNodeReferenceVector(Map<String, double[]> pathwayToData) {
//        Set<String> keySet = pathwayToData.keySet();
//        RandomData random = new RandomDataImpl(new JDKRandomGenerator());
//        int sampleSize = Math.min(nodeList.size(), keySet.size());
//        Object[] sample = random.nextSample(keySet, sampleSize);
//        for (int i = 0; i < sampleSize; i++) {
//            ReferenceNode node = nodeList.get(i);
//            // Use pathway vectors
//            node.setReferenceVector(pathwayToData.get(sample[i].toString()));
//        }
        // Do another way to use a fixed order 
        List<String> pathways = new ArrayList<String>(pathwayToData.keySet());
        int index = 0;
        for (int i = 0; i < nodeList.size(); i++) {
            if (index == pathways.size())
                index = 0;
            ReferenceNode node = nodeList.get(i);
            String pathway = pathways.get(index);
            node.setReferenceVector(pathwayToData.get(pathway));
            index ++;
        }
        // Use a random data set
        // Need to get the size of an double array
//        int size = 0;
//        for (Iterator<String> it = pathwayToData.keySet().iterator(); it.hasNext();) {
//            String topic = it.next();
//            double[] vector = pathwayToData.get(topic);
//            size = vector.length;
//            break;
//        }
//        for (int i = 0; i < nodeList.size(); i++) {
//            ReferenceNode node = nodeList.get(i);
//            double[] vector = new double[size];
//            for (int j = 0; j < size; j++)
//                vector[j] = Math.random();
//            node.setReferenceVector(vector);
//        }
    }
    
    /**
     * Randomly assign reference data from the assigned map.
     * @param pathwayToData
     */
    public void initNodeReferenceData(Map<String, Set<String>> pathwayToData) {
//        Set<String> keySet = pathwayToData.keySet();
//        RandomData random = new RandomDataImpl(new JDKRandomGenerator());
//        int sampleSize = Math.min(nodeList.size(), keySet.size());
//        Object[] sample = random.nextSample(keySet, sampleSize);
//        for (int i = 0; i < sampleSize; i++) {
//            ReferenceNode node = nodeList.get(i);
//            // Objects in sample are String
//            node.setReferenceData(pathwayToData.get(sample[i].toString()));
//        }
        // Implement a new way to random reference data
        // Find the average size
        Set<String> allDataPoints = new HashSet<String>();
        for (Iterator<String> it = pathwayToData.keySet().iterator(); it.hasNext();) {
            String pathway = it.next();
            Set<String> data = pathwayToData.get(pathway);
            allDataPoints.addAll(data);
        }
        int averageSize = allDataPoints.size() / pathwayToData.size();
        Random seeder = new Random();
        RandomDataImpl sampler = new RandomDataImpl();
        for (ReferenceNode node : nodeList) {
            sampler.reSeedSecure(seeder.nextLong());
            Object[] sample = sampler.nextSample(allDataPoints, averageSize);
            Set<String> set = new HashSet<String>();
            for (Object obj : sample)
                set.add(obj.toString());
            //node.setReferenceData(set);
            node.addInputReferenceData(set);
        }
    }
    
    public ReferenceNode getNodeAt(int x, int y) {
        return posToNode.get(x + "," + y);
    }
    
    public ReferenceNode searchMatchNodeForDist(double[] data) {
        // Try the first node first
        ReferenceNode node = nodeList.get(0);
        double distSqr = node.calculateDistanceSqrForLearning(data);
        ReferenceNode winningNode = node;
        double newDistSqr = 0.0d;
        for (int i = 1; i < nodeList.size(); i++) {
            node = nodeList.get(i);
            newDistSqr = node.calculateDistanceSqrForLearning(data);
            if (newDistSqr < distSqr) {
                winningNode = node;
                distSqr = newDistSqr;
            }
        }
        return winningNode;
    }
    
    public ReferenceNode searchMatchNode(double[] data) {
        if (useDotProduct) // Dot product cannot work. Don't know why.
            return searchMatchNodeForDotProduct(data);
        else
            return searchMatchNodeForDist(data);
    }
    
    public ReferenceNode searchMatchNodeForDotProduct(double[] data) {
        ReferenceNode node = nodeList.get(0);
        double product = node.calculateDotProduct(data);
        ReferenceNode winningNode = node;
        double newValue;
        for (int i = 1; i < nodeList.size(); i++) {
            node = nodeList.get(i);
            newValue = node.calculateDotProduct(data);
            if (newValue > product) {
                product = newValue;
                winningNode = node;
            }
        }
        return winningNode;
    }
    
    public void checkNodes() {
        for (int i = 0; i < nodeList.size(); i++) {
            ReferenceNode node = nodeList.get(i);
            System.out.println(i + ": " + node.getReferenceVector());
        }
    }
    
    /**
     * Search the winning ReferenceNode for the specified data.
     * @param data
     * @return
     */
    public ReferenceNode searchMatchNode(Set<String> data) {
        // Try the first first
        ReferenceNode node = nodeList.get(0);
        double similarity = node.calculateSimilarity(data);
        ReferenceNode winningNode = node;
        double newSim = 0.0;
        for (int i = 1; i < nodeList.size(); i++) {
            node = nodeList.get(i);
            newSim = node.calculateSimilarity(data);
            if (newSim > similarity) {
                winningNode = node;
                similarity = newSim;
            }
            if (newSim == 1.0d)
                break; // The highest
        }
        if (similarity == 0.0d) {
            // Try to find an empty node. Otherwise, pathways will be grouped in the
            // first node
            for (ReferenceNode tmp : nodeList) {
                if (tmp.isEmpty()) {
                    winningNode = tmp;
                    break;
                }
            }
        }
        return winningNode;
    }
    
    /**
     * Calculate reference nodes based on the input data
     *
     */
    public void recalculateNodes(int width) {
        int x, y, minX, minY, maxX, maxY;
        List<Set<String>> inputData = new ArrayList<Set<String>>();
        // use squar of distance
        for (ReferenceNode node : nodeList) {
            Point location = node.getLocation();
            x = location.x;
            y = location.y;
            minX = Math.max(0, x - width);
            minY = Math.max(0, y - width);
            maxX = Math.min(SOM.X_LENGTH - 1, x + width);
            maxY = Math.min(SOM.Y_LENGTH - 1, y + width);
            inputData.clear();
            for (x = minX; x < maxX + 1; x++) {
                for (y = minY; y < maxY + 1; y++) {
                    ReferenceNode neighborNode = getNodeAt(x, y);
                    List<Set<String>> data = neighborNode.getInputData();
                    if (data != null)
                        inputData.addAll(data);
                }
            }
            // Set<String> newData = average(inputData);
            // node.setReferenceData(newData);
            node.setInputReferenceData(inputData);
        }
        // Have to empty input nodes
        for (ReferenceNode node : nodeList)
            node.resetInputData();
    }
    
    /**
     * Calculate reference nodes' vectors based on the input sublist.
     * @param width
     */
    public void recalculateNodesVectors(int width) {
        List<double[]> inputVectors = new ArrayList<double[]>();
        // use squar of distance
        for (ReferenceNode node : nodeList) {
            Point location = node.getLocation();
            List<ReferenceNode> neighborNodes = getNodesWithin(location, width);
            inputVectors.clear();
            for (ReferenceNode neighborNode : neighborNodes) {
                List<double[]> data = neighborNode.getInputVectors();
                if (data != null)
                    inputVectors.addAll(data);
            }
            // If nothing in the inputVectors, keep the original value
            if (inputVectors.size() > 0) {
                double[] newVector = average(inputVectors);
                node.setReferenceVector(newVector);
            }
        }
        // Have to empty input nodes
        for (ReferenceNode node : nodeList) {
            node.resetInputVectors();
            node.resetLabels();
        }
    }
    
    public List<ReferenceNode> getNodesWithin(Point location, int width) {
        List<ReferenceNode> rtn = new ArrayList<ReferenceNode>();
        int x = location.x;
        int y = location.y;
        int minX = Math.max(0, x - width);
        int minY = Math.max(0, y - width);
        int maxX = Math.min(SOM.X_LENGTH - 1, x + width);
        int maxY = Math.min(SOM.Y_LENGTH - 1, y + width);
        for (x = minX; x < maxX + 1; x++) {
            for (y = minY; y < maxY + 1; y++) {
                ReferenceNode neighborNode = getNodeAt(x, y);
                rtn.add(neighborNode);
            }
        }
        return rtn;
    }
    
    public double calculateDistWithNeighbors(ReferenceNode node) {
        Point location = node.getLocation();
        // Need to get all neighbor nodes
        List<ReferenceNode> neighborNodes = getNodesWithin(location, 1);
        neighborNodes.remove(node);
        // calculate average distance
        double total = 0.0d;
        for (ReferenceNode neighbor : neighborNodes) {
            if (neighbor == null)
                System.out.println(node.getLocation());
            total += node.calculateDistance(neighbor);
        }
        return total / neighborNodes.size(); 
    }
    
    public double[] average(List<double[]> vectors) {
        double[] first = vectors.get(0);
        double[] average = new double[first.length];
        for (double[] tmp : vectors) {
            for (int i = 0; i < tmp.length; i++) {
                average[i] += tmp[i];
            }
        }
        int size = vectors.size();
        for (int i = 0; i < average.length; i++) {
            average[i] /= size;
        }
        return average;
    }
    
    /**
     * Calculate referenceData based on input data. The value used for referenceData
     * actually is the median of the size of the input data points.
     *
     */
    public Set<String> median(List<Set<String>> inputData) {
        if (inputData == null || inputData.size() == 0) {
            return null; // Cannot do anything
        }
        Collections.sort(inputData, new Comparator<Set<String>>() {
            public int compare(Set<String> set1, Set<String> set2) {
                int size1 = set1.size();
                int size2 = set2.size();
                return size1 - size2;
            }
        });
        // Take the median
        int size = inputData.size();
        int index = size / 2;
        return inputData.get(index);
    }
    
//    public Set<String> average(List<Set<String>> inputData) {
//        if (inputData == null || inputData.size() == 0) {
//            return null; // Cannot do anything
//        }
//        Set<String> rtn = new HashSet<String>(inputData.get(0));
//        // Get all common data points
//        for (Iterator<String> it = rtn.iterator(); it.hasNext();) {
//            String tmp = it.next();
//            for (int i = 1; i < inputData.size(); i++) {
//                Set<String> set = inputData.get(i);
//                if (!set.contains(tmp)) {
//                    it.remove();
//                    break;
//                }
//            }
//        }
//        // Random pick other remaining data points
//        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
//        Set<String> extra = new HashSet<String>();
//        for (Set<String> set : inputData) {
//            if (set.size() > rtn.size()) {
//                Set<String> copy = new HashSet<String>(set);
//                int diff = copy.size() - rtn.size();
//                int number = diff / inputData.size();
//                if (number == 0)
//                    continue;
//                copy.removeAll(rtn);
//                Object[] sample = randomizer.nextSample(copy, number);
//                for (Object obj : sample)
//                    extra.add(obj.toString());
//            }
//        }
//        rtn.addAll(extra);
//        return rtn;
//    }
    
    public void output(OutputStream os) throws IOException {
        PrintStream ps = new PrintStream(os);
        ps.println("SOM results:");
        for (ReferenceNode node : nodeList) {
            Point p = node.getLocation();
            double index = node.calculateAverageDistance();
            ps.println("location:" + p.x + ", " + p.y);
            ps.println("index:"+  index);
            Set<String> labels = node.getLabels();
            if (labels != null) {
                ps.print("label:");
                for (Iterator<String> it = labels.iterator(); it.hasNext();) {
                    String label = it.next();
                    ps.print(label);
                    if (it.hasNext())
                        ps.print(", ");
                }
                ps.println();
            }
            StringBuilder builder = new StringBuilder();
            for (double d : node.getReferenceVector())
                builder.append(d).append(",");
            ps.println("vector:" + builder.toString());
        }
        ps.close();
        os.close();
    }
}
