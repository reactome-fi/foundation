/*
 * Created on Apr 3, 2012
 *
 */
package org.reactome.r3.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

// Need to perform more tests for using this. It seems that this one cannot produce the same
// results from R.
//import javastat.survival.regression.CoxRegression;

/**
 * This class is used CoxRegression class provided by http://www2.thu.edu.tw/~wenwei/ for doing
 * CoxPH analysis to avoid slow R based analysis.
 * @author gwu
 *
 */
public class JStatCoxRegressionWrapper {
    
    public JStatCoxRegressionWrapper() {
    }
    
    /**
     * Run CoxPH analysis for a list of covariates in a univariate way.
     * @param sampleToValues
     * @param sampleToTime
     * @param sampleToEvent
     * @return
     */
    public List<Double> doSurvivalAnalyses(List<Map<String, Double>> sampleToValues,
                                           Map<String, Double> sampleToTime,
                                           Map<String, Double> sampleToEvent) {
        List<Double> rtn = new ArrayList<Double>();
        for (Map<String, Double> sampleToValue : sampleToValues) {
            Double pvalue = doSurvivalAnalysis(sampleToValue, sampleToTime, sampleToEvent);
            rtn.add(pvalue);
        }
        return rtn;
    }
    
    /**
     * Do a survival analysis and return p-value from survival analysis.
     * @param sampleToValue
     * @param sampleToTime
     * @param sampleToEvent
     * @return
     */
    public double doSurvivalAnalysis(Map<String, Double> sampleToValue,
                                     Map<String, Double> sampleToTime,
                                     Map<String, Double> sampleToEvent) {
        // Get shared samples
        Set<String> samples = new HashSet<String>(sampleToValue.keySet());
        samples.retainAll(sampleToTime.keySet());
        samples.retainAll(sampleToEvent.keySet());
        // Create three arrays required by CoxRegression
        double[] times = new double[samples.size()];
        double[] censors = new double[samples.size()];
        double[] values = new double[samples.size()];
        int index = 0;
        for (String sample : samples) {
            times[index] = sampleToTime.get(sample);
            censors[index] = sampleToEvent.get(sample);
            values[index] = sampleToValue.get(sample);
            index ++;
        }
        return -1.0d; 
//        CoxRegression regression = new CoxRegression(times, censors, values);
//        return regression.pValue[0];
    }
    
    /**
     * This method is used to load clinical information
     * @param sampleToEvent
     * @param sampleToTime
     * @param clinFileName
     * @throws IOException
     */
    public Map<String, Double> loadSampleToEvent(String clinFileName,
                                                  String colName) throws IOException {
        Map<String, Double> sampleToEvent = new HashMap<String, Double>();
        Map<String, String> sampleToValue = loadSampleToValue(clinFileName, colName);
        for (String sample : sampleToValue.keySet()) {
            sampleToEvent.put(sample, new Double(sampleToValue.get(sample)));
        }
        return sampleToEvent;
    }
    
    private Map<String, String> loadSampleToValue(String clinFileName,
                                                  String colName) throws IOException {
        Map<String, String> sampleToValue = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        fu.setInput(clinFileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        // Get the index
        int index = 0;
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals(colName)) {
                index = i;
                break;
            }
        }
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (tokens[index].length() > 0 && !tokens[index].equals("NA"))
                sampleToValue.put(tokens[0], tokens[index]);
        }
        fu.close();
        return sampleToValue;
    }
    
    /**
     * Load sample to survival time.
     * @param clinFileName
     * @param colName
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadSampleToSurvival(String clinFileName,
                                                    String colName) throws IOException {
        Map<String, Double> sampleToTime = new HashMap<String, Double>();
        Map<String, String> sampleToValue = loadSampleToValue(clinFileName, colName);
        for (String sample : sampleToValue.keySet()) {
            sampleToTime.put(sample, new Double(sampleToValue.get(sample)));
        }
        return sampleToTime;
    }
    
}
