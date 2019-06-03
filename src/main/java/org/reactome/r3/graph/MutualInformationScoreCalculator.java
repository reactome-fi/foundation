/*
 * Created on Jun 10, 2008
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.reactome.r3.graph.GraphComponent.ScoreCalculator;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to calculate a GraphComponent score based on mutual information.
 * Most of methods in this class are refactored from class BreastCancerArrayDataAnalyzer.
 * @author wgm
 *
 */
public class MutualInformationScoreCalculator implements ScoreCalculator {
    // These cached information is used to calcualte mutual information
    private Map<String, List<Double>> idToValues;
    private Map<String, String> sampleToType;
    private List<String> samples;
    // Used to load file
    private FileUtility fu = new FileUtility();
    
    public MutualInformationScoreCalculator() {
    }
    
    /**
     * To calculate a score based on mutual information. To call this method,
     * this object should be properly initialized by calling two init methods:
     * initIdToValueInfo(String) and initSampelToPhenotypeInfo(String).
     */
    public double calculateScore(GraphComponent comp) {
        List<String> discValues = calculateDiscValues(comp);
        if (discValues == null)
            return Double.NEGATIVE_INFINITY; // No information available
        return calculateScore(discValues);
    }
    
    /**
     * This method is used to calculate the mutable information for a specified gene.
     * @param gene
     * @return
     */
    public double calculateScore(String gene) {
        List<Double> values = idToValues.get(gene);
        if (values == null)
            return Double.NEGATIVE_INFINITY;
        List<String> discValues = discretizeValues(values);
        return calculateScore(discValues);
    }
    
    public double calculateTScore(String gene) throws MathException {
        List<Double> values = idToValues.get(gene);
        if (values == null)
            return Double.NEGATIVE_INFINITY;
        return calculateTScore(values);
    }
    
    public double calculateTScore(Collection<String> genes) throws MathException {
        List<Double> avgValues = averageValues(genes, idToValues);
        if (avgValues == null)
            return 0.0;
        return calculateTScore(avgValues);
    }

    private double calculateTScore(List<Double> values) throws MathException {
        List<Double> sample1 = new ArrayList<Double>();
        List<Double> sample2 = new ArrayList<Double>();
        // Get two samples
        Set<String> types = new HashSet<String>(sampleToType.values());
        List<String> typeList = new ArrayList<String>(types);
        for (int i = 0; i < values.size(); i++) {
            Double value = values.get(i);
            String sample = samples.get(i);
            String type = sampleToType.get(sample);
            int typeIndex = typeList.indexOf(type);
            if (typeIndex == 0)
                sample1.add(value);
            else if (typeIndex == 1)
                sample2.add(value);
        }
        double pvalue = MathUtilities.calculateTTest(sample1, sample2);
        return -Math.log10(pvalue);
    }
    
    /**
     * Calculate score based on a list of discretized values.
     * @param discValues
     * @return
     */
    public double calculateScore(List<String> discValues) {
        Map<String, BinInfo> labelToBinInfo = generateBinInfo(discValues);
        return calculateMutualInformation(labelToBinInfo.values());
    }
    
    /**
     * Generate a list discretized values for a passed GraphComponent.
     * @param comp
     * @return
     */
    public List<String> calculateDiscValues(GraphComponent comp) {
        Set<String> nodes = comp.getAllNodes();
        if (nodes == null || nodes.size() == 0)
            return null; // No any information available
        List<Double> averageValues = averageValues(nodes, idToValues);
        if (averageValues == null)
            return null; // Cannot find score
        List<String> discValues = discretizeValues(averageValues);
        return discValues;
    }
    
    /**
     * Calculate the average expression values for a passed GraphComponent.
     * The returned values are in the order of samples initialized in this calculator.
     * @param comp
     * @return
     */
    public List<Double> averageValues(GraphComponent comp) {
        Set<String> nodes = comp.getAllNodes();
        if (nodes == null || nodes.size() == 0)
            return null; // No any information available
        List<Double> averageValues = averageValues(nodes, idToValues);
        return averageValues;
    }
    
    /**
     * Calculate mututal information based on information in each bins.
     * @param binInfo
     * @param samples
     * @return
     */
    private double calculateMutualInformation(Collection<BinInfo> binInfo) {
        int totalSample = samples.size();
        Map<String, Double> activityToProb = new HashMap<String, Double>();
        Map<String, Double> typeToProb = new HashMap<String, Double>();
        for (BinInfo info : binInfo) {
            Double value = typeToProb.get(info.type);
            if (value == null) 
                typeToProb.put(info.type, (double)info.count);
            else
                typeToProb.put(info.type, info.count + value);
            value = activityToProb.get(info.activity);
            if (value == null)
                activityToProb.put(info.activity, (double)info.count);
            else
                activityToProb.put(info.activity, info.count + value);
        }
        for (Iterator<String> it = activityToProb.keySet().iterator(); it.hasNext();) {
            String activity = it.next();
            Double value = activityToProb.get(activity);
            activityToProb.put(activity,
                               value / totalSample);
        }
        for (Iterator<String> it = typeToProb.keySet().iterator(); it.hasNext();) {
            String type = it.next();
            Double value = typeToProb.get(type);
            typeToProb.put(type,
                           value / totalSample);
        }
        // Calculate mutual information
        double mi = 0.0d;
        for (BinInfo info : binInfo) {
            double pxy = info.count / (double) totalSample;
            double px = activityToProb.get(info.activity);
            double py = typeToProb.get(info.type);
            mi += pxy * Math.log(pxy / (px * py));
        }
        return mi;
    }
    
    /**
     * Generate a map from labels for bins to BinInfo. This is the next to last step
     * to calculate mututal information.
     * @param discValues
     * @return
     */
    private Map<String, BinInfo> generateBinInfo(List<String> discValues) {
        Map<String, BinInfo> labelToBinInfo = new HashMap<String, BinInfo>();
        for (int i = 0; i < discValues.size(); i++) {
            String value = discValues.get(i);
            if (value == null)
                continue;
            String sample = samples.get(i);
            String type = sampleToType.get(sample);
            String label = type + ":" + value;
            BinInfo info = labelToBinInfo.get(label);
            if (info == null) {
                info = new BinInfo();
                info.type = type;
                info.activity = value;
                labelToBinInfo.put(label, info);
            }
            info.count ++;
        }
        return labelToBinInfo;
    }
    
    /**
     * Discretize a list of values in double, encoded in String. 
     * @param values
     * @return
     */
    private List<String> discretizeValues(List<Double> values) {
//        DescriptiveStatistics stat = new DescriptiveStatisticsImpl();
//        for (Double value : values) {
//            stat.addValue(value);
//        }
//        double mean = stat.getMean();
//        double sd = stat.getStandardDeviation();
        double total = 0.0d;
        int count = 0;
        for (Double value : values) {
            if (value == null)
                continue;
            total += value;
            count ++;
        }
        double totalSqr = 0.0d;
        double mean = total / (double) count;
        for (Double value : values) {
            if (value == null)
                continue;
            totalSqr += (value - mean) * (value - mean);
        }       
        double sd = Math.sqrt(totalSqr / (count - 1));
        int binNumber = getBinNumber();
        double[] bins = generateBins(mean, sd, binNumber);
        List<String> discValues = new ArrayList<String>();
        for (Double value : values) {
            if (value == null) {
                discValues.add(null);
                continue;
            }
            Integer bin = null;
            // Find the bin
            for (int j = 0; j < bins.length; j++) {
                if (value < bins[j]) {
                    bin = j;
                    break;
                }
            }
            if (bin == null)
                bin = binNumber - 1; // The largest one
            discValues.add(bin + "");
        }
        return discValues;
    }
    
    /**
     * This method is used to generate bins by using a double array,
     * which is used as threshold values.
     * @param mean
     * @param sd
     * @param binNumber
     * @return
     */
    private double[] generateBins(double mean,
                                  double sd,
                                  int binNumber) {
        double[] bins = new double[binNumber - 1];
        double step = 4.0 * sd / binNumber;
        for (int i = 0; i < binNumber - 1; i++) {
            bins[i] = mean - 2 * sd + step * (i + 1);
        }
        return bins;
    }
    
    /**
     * This method is used to average several row values. The returned values should
     * be ordered based on sample list.
     * @param ids
     * @param idToValues
     * @return
     */
    private List<Double> averageValues(Collection<String> ids,
                                       Map<String, List<Double>> idToValues) {
        List<List<Double>> allData = new ArrayList<List<Double>>();
        // Peek
        List<Double> values = null;
        int lastIndex = 0;
        for (String id : ids) {
            values = idToValues.get(id);
            lastIndex ++;
            if (values != null)
                break;
        }
        if (values == null)
            return null;
        // This means that we reached the end of our ids set to
        // get a non-null value. Don't do average. Just return a copy
        // of values
        if (lastIndex == ids.size()) {
            return new ArrayList<Double>(values);
        }
        // Initialize the above array
        for (int i = 0; i < values.size(); i++)
            allData.add(new ArrayList<Double>());
        // Check the coverage:
        int covered = 0;
        for (String id : ids) {
            values = idToValues.get(id);
            if (values == null)
                continue;
                //throw new IllegalStateException(probeset + " has no data values!");
            for (int i = 0; i < values.size(); i++) {
                List<Double> list = allData.get(i);
                Double value = values.get(i);
                list.add(value);
            }
            covered ++;
        }
        //System.out.println("Covered: " + (double)covered / ids.size());
        List<Double> average = new ArrayList<Double>(allData.size());
        for (List<Double> colValues : allData) {
            double total = 0.0d;
            int c = 0;
            for (Double value : colValues) {
                if (value == null)
                    continue;
                c ++;
                total += value;
            }
            if (c == 0)
                average.add(null);
            else
                average.add(total / c);
        }
        return average;
    }
    
    /**
     * This method is used to load sample to phenotype information from a specified
     * file.
     * @param fileName
     * @throws IOException
     */
    public void initSampleToPhenotypeInfo(String fileName) throws IOException {
        sampleToType = new HashMap<String, String>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#") ||
                line.startsWith("PID")) 
                continue;
            String[] tokens = line.split("\t");
            sampleToType.put(tokens[0],
                             tokens[3]);
        }
        fu.close();
    }
    
    /**
     * Set the sample to phenotype information.
     * @param sampleToType
     */
    public void setSampleToPhenotypeInfo(Map<String, String> sampleToType) {
        this.sampleToType = sampleToType;
    }
    
    /**
     * This method is used to load an id to value (expression values, should be discretized)
     * information from a specified file.
     * @param fileName
     * @throws IOException
     */
    public void initIdToValueInfo(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        extractSample(line);
        // Need to load gene expression values
        idToValues = new HashMap<String, List<Double>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            List<Double> values = new ArrayList<Double>();
            for (int i = 1; i < tokens.length; i++) { 
                if (tokens[i].length() == 0)
                    values.add(null);
                else
                    values.add(new Double(tokens[i]));
            }
            idToValues.put(tokens[0],
                           values);
        }
        fu.close();
    }
    
    public Map<String, List<Double>> getIdToValues() {
        return this.idToValues;
    }
    
    public void setIdToValues(Map<String, List<Double>> idToValues) {
        this.idToValues = idToValues;
    }
    
    public List<String> getSamples() {
        return this.samples;
    }
    
    public void setSamples(List<String> samples) {
        this.samples = samples;
    }
    
    public Map<String, String> getSampleToType() {
        return this.sampleToType;
    }
    
    /**
     * A helper method to extract sample information from a title line.
     */
    private void extractSample(String titleLine) {
        String[] tokens = titleLine.split("\t");
        samples = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++)
            samples.add(tokens[i]);
    }
    
    /**
     * A helper method to calculate the bin number for discretizing.
     * @return
     */
    private int getBinNumber() {
        double log = Math.log(samples.size()) / Math.log(2.0d);
        return (int)log + 1;
    }
    
    /**
     * A very simple data strucutre to catch bin information.
     * @author wgm
     *
     */
    private class BinInfo {
       
        private String type;
        private String activity;
        private int count;
        
        public String label() {
            return type + ":" + activity;
        }
        
    }
    
}
