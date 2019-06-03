/*
 * Created on Jul 11, 2007
 *
 */
package org.reactome.r3.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to encode a matrix for pathway documents. 
 * @author guanming
 *
 */
public class PathwayMatrix {
    private Map<String, double[]> topicToVector;
    
    public PathwayMatrix() {
    }
    
    private Map<String, Map<String, Integer>> loadTopicToIDNumber() throws IOException {
        Map<String, Map<String, Integer>> topicToIdNumber = new HashMap<String, Map<String, Integer>>();
        FileUtility fu = new FileUtility();
        String fileName = R3Constants.RESULT_DIR + "TopicIDNumber.txt";
        fu.setInput(fileName);
        // Title line
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Map<String, Integer> idToNumber = topicToIdNumber.get(tokens[0]);
            if (idToNumber == null) {
                idToNumber = new HashMap<String, Integer>();
                topicToIdNumber.put(tokens[0],
                                    idToNumber);
            }
            idToNumber.put(tokens[1],
                           Integer.parseInt(tokens[2]));
        }
        fu.close();
        return topicToIdNumber;
    }
    
    public void loadData() throws IOException {
        Map<String, Map<String, Integer>> topicToIdNumber = loadTopicToIDNumber();
        // Used to get a list of all genes used in these topics.
        // These genes are used as a vocabulary in our pathway documents
        Set<String> ids = new HashSet<String>();
        for (Iterator<String> it = topicToIdNumber.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            Map<String, Integer> idToNumber = topicToIdNumber.get(topic);
            ids.addAll(idToNumber.keySet());
        }
        // Sort it to get an ordered view
        List<String> idList = new ArrayList<String>(ids);
        Collections.sort(idList);
        // Want to get IDF (inverse document frequency)
        double[] idfs = new double[idList.size()];
        int df = 0;
        int index = 0;
        for (String id : idList) {
            df = 0;
            for (Iterator<String> it = topicToIdNumber.keySet().iterator(); 
                 it.hasNext();) {
                Map<String, Integer> idToNumber = topicToIdNumber.get(it.next());
                if (idToNumber.containsKey(id))
                    df ++;
            }
            idfs[index] = Math.log10((double)topicToIdNumber.size() / df);
            index ++;
        }
        // Arrange ids in topics as the order of idList
        Map<String, int[]> topicToNumber = new HashMap<String, int[]>();
        for (Iterator<String> it = topicToIdNumber.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            Map<String, Integer> idToNumber = topicToIdNumber.get(topic);
            int[] numberArray = new int[idList.size()];
            topicToNumber.put(topic, numberArray);
            // Fill up the values to numberArray
            for (int i = 0; i < idList.size(); i++) {
                Integer value = idToNumber.get(idList.get(i));
                if (value == null)
                    numberArray[i] = 0;
                else
                    numberArray[i] = value.intValue();
            }
        }
        // Convert topicToNumber as frequencies
        topicToVector = new HashMap<String, double[]>();
        for (Iterator<String> it = topicToNumber.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            int[] numberArray = topicToNumber.get(topic);
            int total = 0;
            for (int number : numberArray)
                total += number;
            // Convert to percentages
            double[] vector = new double[numberArray.length];
            for (int i = 0; i < numberArray.length; i++) {
                // tf * idf (tf is used as frequence: percentage)
                vector[i] = (double) numberArray[i] / total * idfs[i];
            }
            topicToVector.put(topic, vector);
        }
    }
    
    public Map<String, double[]> getTopicToVector() {
        return topicToVector;
    }
    
    public double[] getVector(String topic) {
        if (topicToVector == null)
            return null;
        return topicToVector.get(topic);
    }    
    
    @Test
    public void generateDataFile() throws IOException {
        loadData();
        FileUtility fu = new FileUtility();
        String fileName = R3Constants.RESULT_DIR + "PathwayVectors.txt";
        fu.setOutput(fileName);
        int i = 1;
        StringBuilder builder = new StringBuilder();
        for (Iterator<String> it = topicToVector.keySet().iterator(); it.hasNext();) {
            String topic = it.next();
            double[] vector = topicToVector.get(topic);
            builder.append(i).append("\t");
            for (double d : vector)
                builder.append(d).append("\t");
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkData() throws IOException {
        loadData();
        System.out.println("Total Pathways: " + topicToVector.size());
        for (Iterator<String> it = topicToVector.keySet().iterator(); it.hasNext();)
            System.out.println(it.next());
    }
}
