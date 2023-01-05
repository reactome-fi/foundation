/*
 * Created on Jan 11, 2011
 *
 */
package org.reactome.cancer.base;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

/**
 * This class group a set of methods that are related to gene expression data analysis.
 * @author wgm
 *
 */
public class CancerGeneExpressionCommon {
	private FileUtility fu = new FileUtility();

	public CancerGeneExpressionCommon() {
	}

	public Map<String, Map<String, Double>> loadGeneExp(String fileName) throws IOException {
		return loadGeneExp(fileName, true);
	}

	/**
	 * Load gene expression data per gene. The value is stored as gene to sample to value.
	 * @return
	 * @throws IOException
	 */
	private Map<String, Map<String, Double>> loadGeneExp(String fileName,
	                                                     boolean escapeRowContainNA) throws IOException {
		Map<String, Map<String, Double>> geneToData = new HashMap<String, Map<String,Double>>();
		fu.setInput(fileName);
		// Sample list
		String line = fu.readLine();
		List<String> sampleList = new ArrayList<String>();
		String[] tokens = line.split("\t");
		for (String token : tokens) {
			String sample = token.replaceAll("\"", "");
			//if (sample.length() > 12)
			//    sample = sample.substring(0, 12);
			sampleList.add(sample);
		}
		// Starting parsing
		while ((line = fu.readLine()) != null) {
			//            if (line.contains("NA"))
			//                continue; // Don't want any genes containing "NA".
			//
			int index = line.indexOf("\t");
			if (line.substring(index + 1).contains("NA") && escapeRowContainNA)
				continue; // Don't want any genes with values containing "NA".
			//            System.out.println(line);
			tokens = line.split("\t");
			// The first one is gene name
			String gene = tokens[0].replace("\"", "");
			Map<String, Double> sampleToValue = new HashMap<String, Double>();
			for (int i = 1; i < tokens.length; i++) {
				if (tokens[i].equals("NA") || tokens[i].length() == 0)
					continue; // Just escape these values
				String sample = sampleList.get(i); // The first sample has been checked.
				Double value = new Double(tokens[i]);
				sampleToValue.put(sample, value);
			}
			if (geneToData.containsKey(gene)) {
				System.out.println("Duplicated gene: " + gene);
			}
			geneToData.put(gene, sampleToValue);
		}
		fu.close();
		// These two genes having NA values in exp
		//geneToData.remove("C1ORF129");
		//geneToData.remove("LRRC50");
		return geneToData;
	}

	/**
	 * This is the actual method to calculate Pearson correlation for FIs based on provided
	 * gene gene data set in a Map.
	 * @param geneToSampleToValue
	 * @param fis
	 * @param useAbsoluteValue
	 * @param sampleEscapePattern
	 * @return
	 */
	public Set<String> calculateGeneExpCorrForFIs(Map<String, Map<String, Double>> geneToSampleToValue,
	                                              Set<String> fis,
	                                              Boolean useAbsoluteValue,
	                                              String sampleEscapePattern) {
		Set<String> fisWithCorrs = new HashSet<String>();
		int index = 0;
		String gene1, gene2;
		List<Double> values1 = new ArrayList<Double>();
		List<Double> values2 = new ArrayList<Double>();
		//        long time1 = System.currentTimeMillis();
		for (String fi : fis) {
			//            fisWithCorrs.add(fi + "\t1.0");
			//            if (true)
			//                continue;
			index = fi.indexOf("\t");
			gene1 = fi.substring(0, index);
			Map<String, Double> sampleToValue1 = geneToSampleToValue.get(gene1);
			if (sampleToValue1 == null)
				continue;
			gene2 = fi.substring(index + 1);
			Map<String, Double> sampleToValue2 = geneToSampleToValue.get(gene2);
			if (sampleToValue2 == null)
				continue;
			for (String sample : sampleToValue1.keySet()) {
				if (sampleEscapePattern != null && sample.contains(sampleEscapePattern))
					continue; // Escape normal samples
				Double value1 = sampleToValue1.get(sample);
				Double value2 = sampleToValue2.get(sample);
				// Always escape if there is any null value
				if (value1 == null || value2 == null)
					continue;
				values1.add(value1);
				values2.add(value2);
			}
			double corr = MathUtilities.calculatePearsonCorrelation(values1, values2);
			// The following check is very slow. Need to depends a pre-processing to 
			// remove any no correct values in the array data set.
			//            // Just in case
			//            if ((corr + "").equals(nan)) // This is much faster than Double.isNaN(double)
			//                continue;
			if (useAbsoluteValue)
				fisWithCorrs.add(fi + "\t" + Math.abs(corr));
			else
				fisWithCorrs.add(fi + "\t" + corr);
			values1.clear();
			values2.clear();
		}
		//        long time2 = System.currentTimeMillis();
		//        System.out.println("Time for correaltion caluation: " + (time2 - time1));
		return fisWithCorrs;
	}
}
