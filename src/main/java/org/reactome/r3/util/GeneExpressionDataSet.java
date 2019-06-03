/*
 * Created on Jan 21, 2011
 *
 */
package org.reactome.r3.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * This helper class is used to model a gene expression data set
 * @author wgm
 *
 */
public class GeneExpressionDataSet {
    // Gene names or probeset names
    private List<String> featureList;
    // List of sample names
    private List<String> sampleList;
    // Actual values: List<Double> should be values for samples in each features.
    private List<List<Double>> values;
    // Platform used to generate this GeneExpressionDataSet
    private String platform;
    
    public GeneExpressionDataSet() {
    }
    
    public void setPlatform(String platform) {
        this.platform = platform;
    }
    
    public String getPlatform() {
        return this.platform;
    }

    public List<String> getFeatureList() {
        return featureList;
    }

    public void setFeatureList(List<String> featureList) {
        this.featureList = featureList;
    }

    public List<String> getSampleList() {
        return sampleList;
    }
    
    public void addSample(String sample) {
        if (sampleList == null)
            sampleList = new ArrayList<String>();
        sampleList.add(sample);
    }
    
    public void addFeature(String feature) {
        if (featureList == null)
            featureList = new ArrayList<String>();
        featureList.add(feature);
    }

    public void setSampleList(List<String> sampleList) {
        this.sampleList = sampleList;
    }

    public List<List<Double>> getValues() {
        return values;
    }

    public void setValues(List<List<Double>> values) {
        this.values = values;
    }
    
    /**
     * Drop a feature from this GeneExpressionDataSet object.
     * @param feature
     */
    public void dropFeature(String feature) {
        int featureIndex = featureList.indexOf(feature);
        dropFeature(featureIndex);
    }
    
    public void dropFeature(int featureIndex) {
        values.remove(featureIndex);
        featureList.remove(featureIndex);
    }
    
    public void replaceFeatureName(String oldName, String newName) {
        int featureIndex = featureList.indexOf(oldName);
        featureList.set(featureIndex, newName);
    }
    
    public Double getValue(String feature, String sample) {
        int featureIndex = featureList.indexOf(feature);
        int sampleIndex = sampleList.indexOf(sample);
        List<Double> list = values.get(featureIndex);
        return list.get(sampleIndex);
    }
    
    public void addValuesForFeature(List<Double> featureValues) {
        if (values == null)
            values = new ArrayList<List<Double>>();
        values.add(featureValues);
    }
    
    public List<Double> getValuesForFeature(int featureIndex) {
        return values.get(featureIndex);
    }
    
    /**
     * Export this object into a local file.
     * @param fileName
     * @throws IOException
     */
    public void export(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        // First header line
        StringBuilder builder = new StringBuilder();
        builder.append("Feature");
        for (String sample : sampleList)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (int i = 0; i < featureList.size(); i++) {
            String feature = featureList.get(i);
            List<Double> featureValues = values.get(i);
            builder.append(feature);
            for (Double value : featureValues) {
                builder.append("\t");
                if (value == null)
                    builder.append("");
                else
                    builder.append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    public void logTransformation() {
        for (List<Double> list : values) {
            for (int i = 0; i < list.size(); i++) {
                Double value = list.get(i);
                if (value != null) {
                    value = Math.log(value);
                    list.set(i, value);
                }
            }
        }
    }
    
    public void sampleWiseZscoreTransformation() {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (int i = 0; i < sampleList.size(); i++) {
            stat.clear();
            for (List<Double> list : values) {
                Double value = list.get(i);
                if (value != null)
                    stat.addValue(value);
            }
            double mean = stat.getMean();
            double sd = stat.getStandardDeviation();
            for (List<Double> list : values) {
                Double value = list.get(i);
                if (value != null) {
                    value = (value - mean) / sd;
                    list.set(i, value);
                }
            }
        }
    }
    
    public void sampleWiseRankTransformation() {
    	List<Double> listOfSample = new ArrayList<Double>();
    	for (int i = 0; i < sampleList.size(); i++) {
    		for (List<Double> list : values) {
    			Double value = list.get(i);
    			if (value != null)
    				listOfSample.add(value);
    		}
    		Collections.sort(listOfSample, new Comparator<Double>() {
    			public int compare(Double value1, Double value2) {
    				return value2.compareTo(value1); // A a reverse order
    			}
    		});
    		for (List<Double> list : values) {
    			Double value = list.get(i);
    			if (value != null) {
    				int index = listOfSample.indexOf(value);
    				double rank = (double) (index + 1) / listOfSample.size();
    				list.set(i, rank);
    			}
    		}
    		listOfSample.clear();
    	}
    }
    
    public void geneWiseMedianTransformation() {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (int i = 0; i < values.size(); i++) {
            List<Double> geneValues = values.get(i);
            // Try to get median
            stat.clear();
            for (Double d : geneValues) {
                if (d != null) {
                    stat.addValue(d);
                }
            }
            double median = stat.getPercentile(0.50d);
            for (int j = 0; j < geneValues.size(); j++) {
                Double d = geneValues.get(j);
                if (d == null)
                    continue;
                geneValues.set(j, 
                               d - median);
            }
        }
    }
    
    public void zscoreTansformation() {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (List<Double> list : values) {
            for (int i = 0; i < list.size(); i++) {
                Double value = list.get(i);
                if (value != null) {
                    stat.addValue(value);
                }
            }
        }
        double mean = stat.getMean();
        double sd = stat.getStandardDeviation();
        System.out.println("Mean: " + mean);
        System.out.println("SD: " + sd);
        for (List<Double> list : values) {
            for (int i = 0; i < list.size(); i++) {
                Double value = list.get(i);
                if (value != null) {
                    value = (value - mean) / sd;
                    list.set(i, value);
                }
            }
        }
    }
    
    /**
     * This method is used to merge two gene expression data set together. The merge data set 
     * should have the same sample list.
     * @param another
     */
    public void merge(GeneExpressionDataSet another) {
        if (!sampleList.equals(another.getSampleList()))
            throw new IllegalArgumentException("Two merging GeneExpressionDataSet have different sample lists!");
        featureList.addAll(another.getFeatureList());
        values.addAll(another.getValues());
    }
    
    /**
     * Return a subset of values from this GeneExpressionDataSet.
     * @param features
     * @return
     */
    public GeneExpressionDataSet selectFeatures(List<String> features) {
        GeneExpressionDataSet rtn = new GeneExpressionDataSet();
        rtn.platform = this.platform;
        rtn.featureList = new ArrayList<String>(features);
        rtn.sampleList = new ArrayList<String>(sampleList);
        List<List<Double>> rtnValues = new ArrayList<List<Double>>();
        // Need to copy values
        for (String feature : features) {
            int index = this.featureList.indexOf(feature);
            if (index == -1)
                throw new IllegalArgumentException(feature + " cannot be found!");
            List<Double> values = this.values.get(index);
            rtnValues.add(new ArrayList<Double>(values));
        }
        rtn.values = rtnValues;
        return rtn;
    }
    
    /**
     * Use this method to conver this GeneExpressionDataSet into a map of sampleToFeatureToValue.
     * @return
     */
    public Map<String, Map<String, Double>> convertToSampleToFeatureToValueMap() {
        Map<String, Map<String, Double>> sampleToFeatureToValue = new HashMap<String, Map<String,Double>>();
        for (int i = 0; i < sampleList.size(); i++) {
            String sample = sampleList.get(i);
            Map<String, Double> featureToValue = new HashMap<String, Double>();
            sampleToFeatureToValue.put(sample, featureToValue);
            for (int j = 0; j < featureList.size(); j++) {
                String feature = featureList.get(j);
                List<Double> sampleValues = values.get(j);
                Double value = sampleValues.get(i);
                featureToValue.put(feature, value);
            }
        }
        return sampleToFeatureToValue;
    }
    
    /**
     * This class method is used to convert a sampleToFeatureToValue map into a GeneExpressionDataSet object.
     */
    public static GeneExpressionDataSet convertSampleToFeatureToValueMap(Map<String, Map<String, Double>> sampleToFeatureToValue) {
        GeneExpressionDataSet dataset = new GeneExpressionDataSet();
        // Get sample list
        List<String> sampleList = new ArrayList<String>(sampleToFeatureToValue.keySet());
        // Get feature list: To make sure all features should be listed, we need to check all maps.
        Set<String> featureSet = new HashSet<String>();
        for (Map<String, Double> featureToValue : sampleToFeatureToValue.values()) {
            for (String feature : featureToValue.keySet())
                featureSet.add(feature);
        }
        List<String> featureList = new ArrayList<String>(featureSet);
        dataset.setSampleList(sampleList);
        dataset.setFeatureList(featureList);
        // Generate values
        List<List<Double>> allValues = new ArrayList<List<Double>>();
        for (String feature : featureList) {
            List<Double> values = new ArrayList<Double>();
            for (String sample : sampleList) {
                Map<String, Double> featureToValue = sampleToFeatureToValue.get(sample);
                Double value = featureToValue.get(feature);
                values.add(value);
            }
            allValues.add(values);
        }
        dataset.setValues(allValues);
        return dataset;
    }
    
}
