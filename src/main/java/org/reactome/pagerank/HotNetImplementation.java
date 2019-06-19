/*
 * Created on Feb 25, 2013
 *
 */
package org.reactome.pagerank;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.util.InteractionUtilities;

import cern.colt.matrix.DoubleMatrix2D;

/**
 * This class implements the HotNet method published in this paper:
 * DISCOVERY OF MUTATED SUBNETWORKS ASSOCIATED WITH CLINICAL DATA IN CANCER.
 * @author gwu
 *
 */
public class HotNetImplementation {
    private final boolean DEBUG = false;
    
    private EdgeScoreMethod edgeScoreMethod = EdgeScoreMethod.MAX;
    private int permutation = 1000;
    private Double fdrCutoff = 0.25d;
    
    public HotNetImplementation() {
    }
    
    public Double getFdrCutoff() {
        return fdrCutoff;
    }

    public void setFdrCutoff(Double fdrCutoff) {
        this.fdrCutoff = fdrCutoff;
    }
    
    public int getPermutation() {
        return permutation;
    }

    public void setPermutation(int permutation) {
        this.permutation = permutation;
    }

    public EdgeScoreMethod getEdgeScoreMethod() {
        return edgeScoreMethod;
    }

    public void setEdgeScoreMethod(EdgeScoreMethod edgeScoreMethod) {
        this.edgeScoreMethod = edgeScoreMethod;
    }

    /**
     * Run the HotNet method. Make sure the genes sorted in the sortedGenes should be the
     * same as listed in the pre-generated heat kernel file.
     */
    public List<HotNetModule> searchForModules(String heatKernelFile,
                                              Set<String> fis,
                                              Map<String, Double> geneToScore,
                                              Double delta) throws IOException {
        HeatKernelCalculator helper = new HeatKernelCalculator();
        DoubleMatrix2D heatKernel = helper.loadHeatKernelMatrix(heatKernelFile);
        if (DEBUG)
            System.out.println("Pre-generated heat kernel is loaded!");
        
        return searchForModules(heatKernel, 
                                fis, 
                                geneToScore, 
                                delta);
    }



    public List<HotNetModule> searchForModules(DoubleMatrix2D heatKernel,
                                               Set<String> fis,
                                               Map<String, Double> geneToScore,
                                               Double delta) {
        // Make sure delta is not null
        if (delta == null)
            throw new IllegalArgumentException("delta should not be null!");
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> sortedGenes = new ArrayList<String>(genes);
        Collections.sort(sortedGenes);
        if (sortedGenes.size() != heatKernel.rows())
            throw new IllegalArgumentException("The number of genes in the fis is not the same as in the pre-generated HeatKernel!");
        // Use a map for quick fetching of gene's index in the list. This can reduce time for seraching
        // gene's index in the list
        Map<String, Integer> geneToIndex = new HashMap<String, Integer>();
        for (int i = 0; i < sortedGenes.size(); i++) {
            geneToIndex.put(sortedGenes.get(i), i);
        }
        // This should be a sorted list ordered based on size.
        List<Set<String>> modules = searchModules(heatKernel,
                                                  geneToIndex,
                                                  fis, 
                                                  geneToScore,
                                                  delta);
        if (DEBUG)
            System.out.println("Found modules: " + modules.size());
        // Permutation tests for calculating P-values and FDRs
        List<List<Integer>> nullDist = permutationTest(heatKernel,
                                                       geneToIndex,
                                                       fis,
                                                       geneToScore,
                                                       delta);
        List<HotNetModule> rtn = new ArrayList<HotNetModule>();
        for (int i = 0; i < modules.size(); i++) {
            Set<String> module = modules.get(i);
            HotNetModule hotNetModule = new HotNetModule();
            hotNetModule.setGenes(module);
            double pvalue = calculatePValue(module, modules, nullDist);
            hotNetModule.setPvalue(pvalue);
            double fdr = calculateFDR(module, modules, nullDist);
            hotNetModule.setFdr(fdr);
            rtn.add(hotNetModule);
        }
        return rtn;
    }
    
    private double calculatePValue(Set<String> module,
                                   List<Set<String>> modules,
                                   List<List<Integer>> nullDist) {
        // Real number with size > module.size();
        int realNumber = 0;
        for (Set<String> module1 : modules) {
            if (module1.size() >= module.size())
                realNumber ++;
        }
        List<Integer> filteredSizes = new ArrayList<Integer>();
        for (List<Integer> list : nullDist) {
            int count = 0;
            for (Integer size : list) {
                if (size >= module.size())
                    count ++;
            }
            filteredSizes.add(count);
        }
        Collections.sort(filteredSizes, new Comparator<Integer>() {
            public int compare(Integer size1, Integer size2) {
                return size2 - size1;
            }
        });
        // Pick the first number
        Integer size = filteredSizes.get(0);
        if (size < realNumber)
            return 0.0d;
        // Get p-value
        for (int i = 0; i < filteredSizes.size(); i++) {
            size = filteredSizes.get(i);
            if (size < realNumber)
                return (double) (i + 1) / filteredSizes.size();
        }
        return 1.0d; // All random list has more modules
    }
    
    /**
     * A simple way to calculate FDR by the following formula based on module sizes:
     * (nullNumber / total_null_Number)/(realNumber / total_real_module)
     * @param module
     * @param modules
     * @param nullDist
     * @return
     */
    private double calculateFDR(Set<String> module,
                                List<Set<String>> modules,
                                List<List<Integer>> nullDist) {
        // Real number with size > module.size();
        int realNumber = 0;
        for (Set<String> module1 : modules) {
            if (module1.size() >= module.size())
                realNumber ++;
        }
        int totalNullNumber = 0;
        int nullNumber = 0;
        for (List<Integer> list : nullDist) {
            for (Integer count : list) {
                totalNullNumber ++;
                if (count >= module.size())
                    nullNumber ++;
            }
        }
        double fdr = ((double) nullNumber / totalNullNumber) / ((double)realNumber / modules.size());
        if (fdr > 1.0d)
            fdr = 1.0d; // Should be the largest
        return fdr;
    }
    
//    private double selectDelta(DoubleMatrix2D matrix,
//                               Map<String, Integer> geneToIndex,
//                               Set<String> fis,
//                               Map<String, Double> geneToScore) {
//        SummaryStatistics stat = new SummaryStatistics();
//        int size = 5;
//        double chosen = 0.0d;
//        double max = Double.MIN_NORMAL;
//        // check 10 delta values
//        for (double delta = 0.01d; delta <= 0.1d; delta += 0.01) {
//            // Use 100 permutations
//            for (int i = 0; i < 100; i++) {
//                // Since we need to use binary search for performance reason,
//                // we shuffle gene scores here
//                Map<String, Double> randomGeneToScore = generateRandomGeneToScore(geneToScore);
//                List<Set<String>> modules = searchModules(matrix, 
//                                                          geneToIndex,
//                                                          fis, 
//                                                          randomGeneToScore,
//                                                          delta);
//                int total = 0;
//                for (Set<String> module : modules) {
//                    if (module.size() >= size)
//                        total ++;
//                }
//                stat.addValue(total);
//            }
//            if (stat.getMean() > max) {
//                chosen = delta;
//                max = stat.getMean();
//            }
//            System.out.println(delta + ": " + stat.getMean());
//            stat.clear();
//        }
//        System.out.println("Selected delta: " + chosen);
//        return chosen;
//    }
    
    /**
     * This implementation is based on the original HotNet Python code. But using
     * a fixed size cutoff, which is 3, only.
     * @param heatKernel
     * @param fis
     * @param geneToScore
     * @return
     */
    public Double selectDeltaViaPreList(DoubleMatrix2D heatKernel,
                                        Set<String> fis,
                                        Map<String, Double> geneToScore) {
        int oldPermutation = this.permutation;
        // Use 100 permutaiton for delta chosen regardless of what the user's choice
        this.permutation = 100;
        int numberOfDelta = 20;
        int maxCount = 0;
        double chosenDelta = 0.0d;
        for (int i = 0; i < numberOfDelta; i++) {
            // For the FI network, it seems that we may need a much smaller delta
            // So add * 0.01
            double checkingDelta = (i + 1) / (double) numberOfDelta * 0.01d;
            if (DEBUG)
                System.out.println("\tchecking delta: " + checkingDelta);
            List<HotNetModule> modules = searchForModules(heatKernel, 
                                                          fis, 
                                                          geneToScore,
                                                          checkingDelta);
            int count = 0;
            for (HotNetModule module : modules) {
                if (module.getFdr() <= fdrCutoff)
                    count ++;
            }
            if (DEBUG)
                System.out.println("\tSignificant modules: " + count);
            if (count > maxCount) {
                maxCount = count;
                chosenDelta = checkingDelta;
            }
        }
        if (DEBUG) {
            System.out.println("Chosen delta: " + chosenDelta);
            System.out.println("Maximum signficant modules (FDR cutoff = 0.25): " + maxCount);
        }
        // Reset back to the old permutation
        this.permutation = oldPermutation; 
        return chosenDelta;
    }
                                        
    
    public Double selectDeltaViaRealData(DoubleMatrix2D heatKernel,
                                         Set<String> fis,
                                         Map<String, Double> geneToScore) {
        // We want to use 100 permutations only for quick performance
        int oldPermutation = this.permutation;
        this.permutation = 100;
        // Use this loop to get an optimum delta for getting the largest number of network modules
        Double delta = 0.0d;
        // Want to run five times in order to get an ten thousandth values
        double currentStep = 1.0d;
        for (int i = 0; i < 5; i++) {
            currentStep /= 10;
            Double currentDelta = selectDeltaViaRealData(heatKernel,
                                                         fis,
                                                         geneToScore,
                                                         delta,
                                                         currentStep);
            if (currentDelta == null || currentDelta == delta)
                break;
            delta = currentDelta;
        }
        this.permutation = oldPermutation;
        return delta;
    }
    
    private Double selectDeltaViaRealData(DoubleMatrix2D heatKernel,
                                          Set<String> fis,
                                          Map<String, Double> geneToScore,
                                          double delta,
                                          double step) {
        if (DEBUG) {
            System.out.println("delta: " + delta);
        }
        Double chosenDelta = 0.0d;
        int maxCount = 0;
        for (int i = -10; i <= 10; i++) {
            double checkingDelta = delta + i * step;
            if (checkingDelta < 0.0d)
                continue;
            if (DEBUG)
                System.out.println("\tchecking delta: " + checkingDelta);
            List<HotNetModule> modules = searchForModules(heatKernel, 
                                                          fis, 
                                                          geneToScore,
                                                          checkingDelta);
            int count = 0;
            for (HotNetModule module : modules) {
                if (module.getFdr() <= fdrCutoff)
                    count ++;
            }
            if (DEBUG)
                System.out.println("\tSignificant modules: " + count);
            if (count > maxCount) {
                maxCount = count;
                chosenDelta = checkingDelta;
            }
        }
        if (DEBUG) {
            System.out.println("Chosen delta: " + chosenDelta);
            System.out.println("Maximum signficant modules: " + maxCount);
        }
        return chosenDelta;
    }
    
    private Map<String, Double> generateRandomGeneToScore(Map<String, Double> geneToScore) {
        List<String> geneList = new ArrayList<String>(geneToScore.keySet());
        Collections.shuffle(geneList);
        List<Double> values = new ArrayList<Double>(geneToScore.values());
        Collections.shuffle(values);
        Map<String, Double> randomGeneToScore = new HashMap<String, Double>();
        for (int j = 0; j < geneList.size(); j++) {
            randomGeneToScore.put(geneList.get(j), values.get(j));
        }
        return randomGeneToScore;
    }

    private List<Set<String>> searchModules(DoubleMatrix2D matrix,
                                            Map<String, Integer> geneToIndex,
                                            Set<String> fis,
                                            Map<String, Double> geneToScore,
                                            double delta) {
        Set<String> keptFIs = new HashSet<String>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            Double score1 = geneToScore.get(gene1);
            if (score1 == null)
                continue;
            String gene2 = fi.substring(index + 1);
            Double score2 = geneToScore.get(gene2);
            if (score2 == null)
                continue;
            Double edgeScore = null;
            if (edgeScoreMethod == EdgeScoreMethod.AVERAGE)
                edgeScore = (score1 + score2) / 2.0d;
            else
                edgeScore = Math.max(score1, score2); // Default method
            int index1 = geneToIndex.get(gene1);
            int index2 = geneToIndex.get(gene2);
            // In our case, these two values are the same
            double weight = Math.min(matrix.get(index1, index1),
                                     matrix.get(index2, index1));
            weight *= edgeScore;
            if (weight >= delta)
                keptFIs.add(fi);
        }
//        System.out.println("Kept FIs: " + keptFIs.size());
        if (keptFIs.size() == 0)
            return new ArrayList<Set<String>>();
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        List<Set<String>> modules = graphAnalyzer.calculateGraphComponents(keptFIs);
        return modules;
    }
    
    private List<List<Integer>> permutationTest(DoubleMatrix2D matrix,
                                                Map<String, Integer> geneToIndex,
                                                Set<String> fis,
                                                Map<String, Double> geneToScore,
                                                double delta) {
//        int permutation = 100;
        List<List<Integer>> rtn = new ArrayList<List<Integer>>();
        long time1 = System.currentTimeMillis();
        for (int i = 0; i < permutation; i++) {
//            if (DEBUG)
//                System.out.println("Permutation " + i + "...");
            Map<String, Double> randomGeneToScore = generateRandomGeneToScore(geneToScore);
            List<Set<String>> randomModules = searchModules(matrix,
                                                            geneToIndex,
                                                            fis, 
                                                            randomGeneToScore,
                                                            delta);
            List<Integer> randomSizes = new ArrayList<Integer>();
            for (Set<String> module : randomModules)
                randomSizes.add(module.size());
            rtn.add(randomSizes);
        }
        long time2 = System.currentTimeMillis();
        if (DEBUG)
            System.out.println("Time for permutation tests: " + (time2 - time1));
        return rtn;
    }
    
    /**
     * Method to get edge score based on two gene scores.
     * @author gwu
     *
     */
    static enum EdgeScoreMethod {
        MAX, // Use the maximum score of two genes in an edge. The original HotNet method
        AVERAGE // Use the average number of two genes in an edge.
    }
}
