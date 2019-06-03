/*
 * Created on Apr 14, 2010
 *
 */
package org.reactome.r3.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to calculate network modularity.
 * @author wgm
 *
 */
public class NetworkModularityCalculator {
    private final boolean debug = false;
    
    public NetworkModularityCalculator() {
    }
    
    public List<Set<String>> convertModuleObjectToSet(List<NetworkModule> modules) {
        List<Set<String>> list = new ArrayList<Set<String>>();
        for (NetworkModule module : modules) {
            list.add(new HashSet<String>(module.getIds()));
        }
        // Sorting from bigger to smaller
        Collections.sort(list, new Comparator<Set<String>>() {
            public int compare(Set<String> set1, Set<String> set2) {
                return set2.size() - set1.size();
            }
        });
        return list;
    }
    
    /**
     * This method is used to do network clustering and calculate modularity and FDR for each individual cluster.
     * SpectralpartitionNetworkCluster class is used for fast performance. The returned FDR values in NetworkModule
     * objects are based on two criteria: module size and modularity. 
     * A random network module should be both have bigger size and higher modularity to be regarded as "better" than
     * an original module.
     * @param fis
     * @param totalGenes gene pool that is used for permutation test.
     * @param permutationNumber
     * @throws Exception
     */
    public List<NetworkModule> calculateClusterModularityWithFDR(Set<String> genes,
                                                                 boolean needRemoveZNF,
                                                                 int permutationNumber) throws Exception {
        // Load background information
        FileUtility fu = new FileUtility();
        Set<String> allFIs = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        if (debug)
            System.out.println("Total FIs: " + allFIs.size());
        if (needRemoveZNF) {
            // Used to remove ZNF cliques
            Collection<String> znfClique = new GraphAnalyzer().searchZNFClique(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
            Set<String> znfFIs = InteractionUtilities.getFIs(znfClique, allFIs);
            allFIs.removeAll(znfFIs);
            if (debug) {
                System.out.println("ZNF clique size: " + znfClique.size());
                System.out.println("FIs involved with ZNF clique: " + znfFIs.size());
                System.out.println("After removing FIs with ZNF clique: " + allFIs.size());
            }
        }
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(allFIs);
        // Get the total number of genes
        genes.retainAll(allGenes); // Only want to consider genes in the FI network
        Set<String> fis = InteractionUtilities.getFIs(genes, allFIs);
        int size = genes.size();
        // Set up required objects
        SpectralPartitionNetworkCluster clustering = new SpectralPartitionNetworkCluster();
        // Handle the original clusters
        List<Set<String>> originalClusters = clustering.cluster(fis);
        List<NetworkModule> originalModules = new ArrayList<NetworkModule>();
        for (int i = 0; i < originalClusters.size(); i++) {
           NetworkModule module = new NetworkModule();
           module.setIndex(i);
           module.setIds(originalClusters.get(i));
           double modularity = calculateModularity(module.getIds(), fis);
           module.setModularity(modularity);
           //module.setWeightedModularity(calculateWeightedModularity(module.getIds(), 
           //                                                         fis,
           //                                                         genes));
           originalModules.add(module);
        }
        if (permutationNumber == 0)
            return originalModules;
        // Start permutation test
        RandomData randomizer = new RandomDataImpl();
        List<NetworkModule> permModules = new ArrayList<NetworkModule>();
        for (int j = 0; j < permutationNumber; j++) {
            if (debug)
                System.out.println("Permutation " + j + "...");
            Set<String> permGenes = MathUtilities.randomSampling(allGenes, 
                                                                 size,
                                                                 randomizer);
            Set<String> fisInSampleGenes = InteractionUtilities.getFIs(permGenes, 
                                                                       allFIs);
            List<Set<String>> permuClusters = clustering.cluster(fisInSampleGenes);
            for (int i = 0; i < permuClusters.size(); i++) {
                Set<String> permCluster = permuClusters.get(i);
                NetworkModule module = new NetworkModule();
                module.setIndex(i);
                module.setIds(permCluster);
                // Don't want to calculate modularity for performance reason!
                //if (debug)
                module.setModularity(calculateModularity(permCluster, fisInSampleGenes));
                //module.setWeightedModularity(calculateWeightedModularity(permCluster, 
                //                                                         fisInSampleGenes, 
                //                                                         permGenes));
                permModules.add(module);
            }
        }
        // calculate FDRs for original modules.
        // Note: originalModule has been sorted based on size
        // Sort original clusters based on weighted modularities
//        Comparator<NetworkModule> comparator = new Comparator<NetworkModule>() {
//            public int compare(NetworkModule module1, NetworkModule module2) { 
//                return module2.getWeightedModularity().compareTo(module1.getWeightedModularity());
//            }
//         };
//        Collections.sort(originalModules, comparator);
//        Collections.sort(permModules, comparator);
        if (debug) {
            String fileName = null;
            if (needRemoveZNF)
                fileName = "tmp/IndividualModuleSizeAndModualarity_" + permutationNumber + "Perm_with_ZNF_removed.txt";
            else
                fileName = "tmp/IndividualModuleSizeAndModualarity_" + permutationNumber + "Perm.txt";
            fu.setOutput(fileName);
            fu.printLine("Index\tSize\tModularity");
            for (int i = 0; i < permModules.size(); i++) {
                NetworkModule module = permModules.get(i);
                fu.printLine(i + "\t" + module.getSize() + "\t" + module.getModularity());
                //if (module.getModularity() >= 0.20)
                //    fu.printLine(module.getIds() + "");
            }
            fu.close();
        }
        int count = 0;
        for (int i = 0; i < originalModules.size(); i++) {
            NetworkModule orig = originalModules.get(i);
            count = 0;
            for (int j = 0; j < permModules.size(); j++) {
                NetworkModule perm = permModules.get(j);
                if (perm.getSize() >= orig.getSize() && perm.getModularity() >= orig.getModularity())
                    count ++;
            }
            double pvalue = (double) count / permModules.size();
            orig.setNomialPvalue(pvalue);
        }
        // Convert pvalue to FDRs
        Collections.sort(originalModules, new Comparator<NetworkModule>() {
            public int compare(NetworkModule module1, NetworkModule module2) {
                return module1.getNomialPvalue().compareTo(module2.getNomialPvalue());
            }
        });
        for (int i = 0; i < originalModules.size(); i++) {
            NetworkModule module = originalModules.get(i);
            double fdr = module.getNomialPvalue() * originalModules.size() / (i + 1);
            if (fdr > 1.0)
                fdr = 1.0;
            module.setFDR(fdr);
        }
        return originalModules;
    }
    
    /**
     * This method is used to calculate modularity based on paper: Newman, MEG & Girvan, M.
     * Physical Review E 69, 026113 (2004). Finding and evaluating community structure in networks.
     * In this paper, the modularity is defined as: Q = Sigma (e - a*a), which can be converted as
     * Q = Sigma (Eii / Et - (Ei / Et) * (Ei / Et)). For a much cleared definition, see Newman, MEJ 
     * PNAS 103: 8577-8582 (2006).
     * @param clusters
     * @param fis
     * @return
     */
    public double calculateModularity(Collection<Set<String>> clusters,
                                      Set<String> fis) {
        return calculateModularity(clusters, fis, false);
    }
    
    /**
     * Calculate the modularity for a specified cluster. This cluster should be 
     * generated from the passed FIs.
     * @param cluster
     * @param fis
     * @return
     */
    public double calculateModularity(Set<String> cluster,
                                      Set<String> fis) {
        return calculateModularity(cluster, fis, null);
    }
    
    private double calculateModularity(Set<String> cluster,
                                       Set<String> fis,
                                       Set<String> genes) {
        int index;
        int insideEdge = 0;
        int allEdge = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            if (cluster.contains(gene1) && cluster.contains(gene2)) {
                insideEdge ++;
                allEdge ++; // Inside edge should count twice
                allEdge ++;
            }
            else if (cluster.contains(gene1) || cluster.contains(gene2)) {
                allEdge ++;
            }
        }
        double tmp = (double) allEdge / (2.0 * fis.size());
        double score = (double) insideEdge / fis.size() - tmp * tmp;
        if (genes != null) {
            // A flag to need to calculate weighted modularity
            score *= (double) cluster.size() / genes.size();
        }
        return score;
    }
    
    public double calculateWeightedModularity(Set<String> cluster,
                                              Set<String> fis,
                                              Set<String> genes) {
        return calculateModularity(cluster, fis, genes);
    }
    
    private double calculateModularity(Collection<Set<String>> clusters,
                                       Set<String> fis,
                                       boolean useWeight) {
        int total = 0;
        if (useWeight) {
            Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
            total = allGenes.size();
        }
        double score = 0.0;
        int insideEdge = 0;
        int allEdge = 0;
        int index = 0;
        for (Set<String> cluster : clusters) {
            // Get edges in this cluster
            insideEdge = 0;
            allEdge = 0;
            for (String fi : fis) {
                index = fi.indexOf("\t");
                String gene1 = fi.substring(0, index);
                String gene2 = fi.substring(index + 1);
                if (cluster.contains(gene1) && cluster.contains(gene2)) {
                    insideEdge ++;
                    allEdge ++; // Inside edge should count twice
                    allEdge ++;
                }
                else if (cluster.contains(gene1) || cluster.contains(gene2)) {
                    allEdge ++;
                }
            }
            double tmp = (double) allEdge / (2.0 * fis.size());
            if (useWeight) {
                double weight = (double) cluster.size() / total;
                score += ((double) insideEdge / fis.size() - tmp * tmp) * weight;
            }
            else
                score += (double) insideEdge / fis.size() - tmp * tmp;
        }
        return score;
    }
    
    /**
     * Based on method calculateModularity() but considering the sized of modules.
     * @param clusters
     * @param fis
     * @return
     */
    public double calculateWeightedModularity(Collection<Set<String>> clusters,
                                              Set<String> fis) {
        return calculateModularity(clusters, fis, true);
    }    
}
