/*
 * Created on Jan 31, 2013
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

import edu.uci.ics.jung.algorithms.cluster.EdgeBetweennessClusterer;
import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;
import edu.uci.ics.jung.algorithms.importance.Ranking;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
//import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;

/**
 * This utility class group all generic graph related methods powered by the Jung
 * graph library.
 * @author gwu
 *
 */
public class JungGraphUtilities {
    
    /**
     * This method is used to label network clusters based on the highest centrality in 
     * each cluster.
     * @param clusters
     * @return
     */
    public List<String> labelNetworkClusters(List<Set<String>> clusters) throws IOException {
        FileUtility fu = new FileUtility();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        List<String> labels = new ArrayList<String>();
        for (Set<String> cluster : clusters) {
            Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
            Map<String, Double> geneToCluster = calculateCentralities(cluster, fisInCluster);
            // Find the highest gene
            String label = null;
            Double highest = Double.NEGATIVE_INFINITY;
            for (String gene : geneToCluster.keySet()) {
                Double centrality = geneToCluster.get(gene);
                if (centrality > highest) {
                    label = gene;
                    highest = centrality;
                }
                else if (centrality == highest) {
                    // Want to keep all genes with the same centrality
                    label += "," + gene;
                }
            }
            labels.add(label);
        }
        return labels;
    }
    
    /**
     * This method is used to calculate centrality for a set of genes based on the provide
     * FIs among these genes. The provided FIs are used to construct the network. EdgeBetweenness
     * algorithm is used. The centralities have been normalized.
     * @param genes
     * @param fisInGenes
     * @return
     */
    public Map<String, Double> calculateCentralities(Set<String> genes,
                                                     Set<String> fisInGenes) {
        Map<String, Double> geneToCentality = new HashMap<String, Double>();
        Graph<String, String> graph = createJungGraph(genes, fisInGenes);
        BetweennessCentrality<String, String> centrality = new BetweennessCentrality<String, String>(graph, true, true);
        centrality.evaluate();
        List<Ranking<?>> rankings = centrality.getRankings();
        double normalizer = (genes.size() - 1) * (genes.size() - 2) / 2.0;
        for (Ranking<?> ranking : rankings) {
            String object = (String) ranking.getRanked();
            if (genes.contains(object)) {
                geneToCentality.put(object, ranking.rankScore / normalizer);
            }
        }
        return geneToCentality;
    }
    
    public List<Set<String>> cluster(Collection<String> genes,
                                     double removeEdgeRatio, 
                                     Set<String> fis) {
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // Want to use genes in our network only
        Set<String> geneCopy = new HashSet<String>(genes);
        geneCopy.retainAll(fiGenes);
        Set<String> fisInGenes = InteractionUtilities.getFIs(geneCopy, 
                                                             fis);
        EdgeBetweennessClusterer<String, String> clusterer = new EdgeBetweennessClusterer<String, String>((int)(fisInGenes.size() * removeEdgeRatio));
        Set<String> connectedGenes =  InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        Graph<String, String> graph = createJungGraph(connectedGenes, fisInGenes);
        Collection<Set<String>> clusters = clusterer.transform(graph);
        List<Set<String>> clusterList = new ArrayList<Set<String>>(clusters);
        Collections.sort(clusterList, new Comparator<Set<String>>() {
            public int compare(Set<String> cluster1, Set<String> cluster2) {
                return cluster2.size() - cluster1.size();
            }
        });
        return clusterList;
    }
    
    public List<Set<String>> cluster(Collection<String> genes, Set<String> fis) {
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // Want to use genes in our network only
        Set<String> geneCopy = new HashSet<String>(genes);
        geneCopy.retainAll(fiGenes);
        Set<String> fisInGenes = InteractionUtilities.getFIs(geneCopy, 
                                                             fis);
        Collection<Set<String>> optimalClusters = null;
        // Don't want to use genes that are not connected to another other genes
        Set<String> connectedGenes = InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        System.out.println("Connected genes: " + connectedGenes.size());
        Graph<String, String> graph = createJungGraph(connectedGenes, fisInGenes);
        double maxModularity = Double.MIN_VALUE;
        Set<Integer> checkedNumbers = new HashSet<Integer>();
        NetworkModularityCalculator modularityCalculator = new NetworkModularityCalculator();
        for (double ratio = 0.0; ratio <= 0.3; ratio += 0.005) {
            int number = (int) (fisInGenes.size() * ratio);
            if (checkedNumbers.contains(number))
                continue;
            checkedNumbers.add(number);
            EdgeBetweennessClusterer<String, String> clusterer = new EdgeBetweennessClusterer<String, String>(number);
            Collection<Set<String>> clusters = clusterer.transform(graph);
            double modularity = modularityCalculator.calculateModularity(clusters,
                                                                         fisInGenes);
            if (modularity > maxModularity) {
                maxModularity = modularity;
                optimalClusters = clusters;
            }
            System.out.println(ratio + ": " + modularity);
            //break;
        }
        List<Set<String>> clusterList = new ArrayList<Set<String>>(optimalClusters);
        Collections.sort(clusterList, new Comparator<Set<String>>() {
            public int compare(Set<String> cluster1, Set<String> cluster2) {
                return cluster2.size() - cluster1.size();
            }
        });
        return clusterList;
    }
    
    /**
     * Create a set of genes and their interactions into a Jung Graph.
     * @param genes
     * @param fisInGenes
     * @return
     */
    public Graph<String, String> createJungGraph(Set<String> genes,
                                                 Set<String> fisInGenes) {
        Graph<String, String> graph = new UndirectedSparseGraph<String, String>();
        for (String gene : genes)
            graph.addVertex(gene);
        int index = 0;
        for (String fi : fisInGenes) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            graph.addEdge(fi, gene1, gene2, EdgeType.UNDIRECTED);
        }
        return graph;
    }
    
    /**
     * Create a Jung graph from a set of FIs.
     * @param fis
     * @return
     */
    public Graph<String, String> createJungGraph(Set<String> fis) {
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        return createJungGraph(genes, fis);
    }
    
}
