/*
 * Created on Jun 28, 2010
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.*;

import org.junit.Test;
import org.reactome.r3.cluster.DistanceCalculator;
import org.reactome.r3.cluster.HierarchicalCluster;
import org.reactome.r3.cluster.HierarchicalClusterNode;
import org.reactome.r3.graph.BreadthFirstSearch.Edge;
import org.reactome.r3.graph.BreadthFirstSearch.TreeNode;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to build an FI network among a gene set. If there are some genes that cannot be linked,
 * linker genes will be used from the whole FI network.
 * @author wgm
 *
 */
public class NetworkBuilderForGeneSet {
    private String fiFileName;
    private String pathwayFiGeneFileName;
    private FileUtility fu;
    // Cache the whole FI set to avoid reloading
    private Set<String> allFIs;
    // Cache
    private BreadthFirstSearch bfs;
    
    public NetworkBuilderForGeneSet() {
        fu = new FileUtility();
        bfs = new BreadthFirstSearch();
    }
    
    public String getPathwayFiGeneFileName() {
        return pathwayFiGeneFileName;
    }

    public void setPathwayFiGeneFileName(String pathwayFiGeneFileName) {
        this.pathwayFiGeneFileName = pathwayFiGeneFileName;
    }

    public String getFiFileName() {
        return fiFileName;
    }
    
    public void setFiFileName(String fiFileName) {
        this.fiFileName = fiFileName;
    }
    
    private Set<String> getAllFIs() throws IOException {
        if (allFIs == null)
            allFIs = fu.loadInteractions(fiFileName);
        return allFIs;
    }
    
    public void setAllFIs(Set<String> allFIs) {
        this.allFIs = allFIs;
    }
    
    public Set<String> constructFINetworkForGeneSet(Collection<String> geneSet) throws IOException {
        return constructFINetworkForGeneSet(geneSet, new HashMap<>());
    }
    
    public Set<String> constructFINetworkForGeneSet(Collection<String> geneSet,
                                                    Map<String, Double> geneToScore) throws IOException {
        List<String> geneList = new ArrayList<String>(geneSet);
        // Filter genes that are not in the FI network
        Set<String> fis = getAllFIs();
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        geneList.retainAll(genesInFIs);
        Map<String, Integer> pairToDistance = calculateDistances(geneList);
        List<HierarchicalClusterNode> newClusters = hierarchicalCluster(geneList, 
                                                                        pairToDistance, 
                                                                        true);
        // Calculate total cluster weight
        HierarchicalClusterNode newFirstNode = newClusters.get(0);
        Set<HierarchicalClusterNode> allNodes = new HashSet<HierarchicalClusterNode>();
        HierarchicalCluster hclust = new HierarchicalCluster();
        hclust.grepAllClusters(newFirstNode, allNodes);
        // These are used to calculate the actual path
        Set<String> spanningFIs = generateSpanFromClusters(geneList,
                                                           allNodes,
                                                           geneToScore);
        // Want to generate all fis
        return generateFIsForGeneSet(fis, 
                                     spanningFIs, 
                                     geneList);
    }
    
    /**
     * This is the main method to construct a FI network for a gene set. Extra linker genes will be 
     * added by using this class.
     * @param geneSet
     * @param geneExpFileName
     * @return
     * @throws IOException
     */
    public Set<String> constructFINetworkForGeneSet(Collection<String> geneSet,
                                                    String geneExpFileName) throws IOException {
        if (geneExpFileName == null)
            return constructFINetworkForGeneSet(geneSet);
        Map<String, Double> geneToTScore = loadGeneExpTValue(geneExpFileName);
        return constructFINetworkForGeneSet(geneSet, geneToTScore);
    }
    
    @Test
    public void testConstructFINetworkForGeneSet() throws IOException {
        setFiFileName(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        String[] genes = {"ARID2", "GPC3"};
        Set<String> fis = constructFINetworkForGeneSet(Arrays.asList(genes));
        System.out.println(fis);
    }
    
//    @Test
//    public void testConstruction() throws Exception {
//        String dir = "/Users/wgm/Documents/Irina/";
//        String fileName = dir + "genesInFI.txt";
//        Set<String> genes = fu.loadInteractions(fileName);
//        String text = InteractionUtilities.joinStringElements(",", genes);
//        System.out.println(text);
//        if (true)
//            return;
//        setFiFileName(R3Constants.GENE_FI_FILE_NAME);
//        constructFINetworkForGeneSet(genes, null);
//    }
    
    private Set<String> generateSpanFromClusters(List<String> ids,
                                                 Set<HierarchicalClusterNode> allClusters,
                                                 Map<String, Double> geneToScore) throws IOException {
        // These are used to calculate shortest path
        Set<String> fis = getAllFIs();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        // Want to sort edges for BFS tree
        //FIFileAnalyzer fiAnalyzer = new FIFileAnalyzer();
        //Set<String> pathwayFIs = fiAnalyzer.loadPathwayAndTFTargetFIs();
        //Set<String> pathwayFIs = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> pathwayFIs = null;
        if (pathwayFiGeneFileName != null)
            pathwayFIs = fu.loadInteractions(pathwayFiGeneFileName);
        else
            pathwayFIs = new HashSet<String>();
        if (geneToScore == null) // Use an empty map for easy API calling
            geneToScore = new HashMap<String, Double>();
        EdgeSorter edgeSorter = new EdgeSorter();
        edgeSorter.geneToPartners = geneToPartners;
        edgeSorter.pathwayFIs = pathwayFIs;
        edgeSorter.geneToTValue = geneToScore; 
        edgeSorter.mutatedGenes = ids;
        //edgeSorter.createLocalParntersMap();
        bfs.setEdgeSorter(edgeSorter);
        Map<String, List<String>> pairToPath = bfs.generateShortestPath(ids, 
                                                                        nodeToEdges);
        return generateSpanFromClusters(allClusters, 
                                        pairToPath,
                                        geneToPartners,
                                        pathwayFIs,
                                        geneToScore,
                                        ids);
    }
    
    private Set<String> generateSpanFromClusters(Set<HierarchicalClusterNode> clusters,
                                                 Map<String, List<String>> pairToPath,
                                                 Map<String, Set<String>> geneToPartners,
                                                 Set<String> pathwayFIs,
                                                 Map<String, Double> geneToTValue,
                                                 List<String> ids) {
        Set<String> spanningFIs = new HashSet<String>();
        PathSorter pathSorter = new PathSorter();
        pathSorter.geneToPartners = geneToPartners;
        pathSorter.pathwayFIs = pathwayFIs;
        pathSorter.mutatedGenes = new HashSet<String>(ids);
        Set<String> linkerGenes = new HashSet<String>();
        pathSorter.linkerGenes = linkerGenes;
        pathSorter.createLocalParntersMap();
        pathSorter.geneToTValue = geneToTValue;
        List<HierarchicalClusterNode> clusterList = new ArrayList<HierarchicalClusterNode>(clusters);
        Collections.sort(clusterList, new Comparator<HierarchicalClusterNode>() {
            public int compare(HierarchicalClusterNode node1, HierarchicalClusterNode node2) {
                return node2.ids.size() - node1.ids.size();
            }
        });
        
        for (HierarchicalClusterNode node : clusterList) {
            if (node.getChildNode1() == null || node.getChildNode2() == null)
                continue; // The bottom node
            List<String> shortestPath = calculateShortestPath(node.getChildNode1(),
                                                              node.getChildNode2(), 
                                                              pairToPath,
                                                              pathSorter);
            List<String> fis1 = convertPathToFIs(shortestPath);
            spanningFIs.addAll(fis1);
        }
        return spanningFIs;
    }
    
    private List<String> calculateShortestPath(HierarchicalClusterNode node1,
                                               HierarchicalClusterNode node2,
                                               Map<String, List<String>> pairToPath,
                                               PathSorter pathSorter) {
        List<String> miniPath = null;
        for (String id1 : node1.ids) {
            for (String id2 : node2.ids) {
                String key = generateSortedKey(id1, id2);
                List<String> shortestPath = pairToPath.get(key);
                if (miniPath == null)
                    miniPath = shortestPath;
                else {
                    int compare = pathSorter.compare(miniPath, shortestPath);
                    if (compare > 0)
                        miniPath = shortestPath;
                }
            }
        }
        return miniPath;
    }
    
    private String generateSortedKey(String id1, String id2) {
        int compare = id1.compareTo(id2);
        if (compare < 0)
            return id1 + "\t" + id2;
        return id2 + "\t" + id1;
    }
    
    private List<String> convertPathToFIs(List<String> path) {
        List<String> fis = new ArrayList<String>();
        for (int i = 0; i < path.size() - 1; i ++) {
            String id1 = path.get(i);
            String id2 = path.get(i + 1);
            // Want to make sure all FIs in in prefixed order
            int compare = id1.compareTo(id2);
            if (compare < 0)
                fis.add(id1 + "\t" + id2);
            else
                fis.add(id2 + "\t" + id1);
        }
        return fis;
    }
    
    public Map<String, Double> loadGeneExpTValue(String fileName) throws IOException {
        fu.setInput(fileName);
        Map<String, Double> rtn = new HashMap<String, Double>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            rtn.put(tokens[0],
                    new Double(tokens[1])); 
        }
        return rtn;
    }
    
    private List<HierarchicalClusterNode> hierarchicalCluster(Collection<String> ids,
                                                              Map<String, Integer> pairToPath, // Use String - Integer or String - double
                                                              boolean useShortestPath) {
        // Make sure pairToPath should be String - Double
        final Map<String, Double> pairToPath1 = new HashMap<String, Double>();
        for (Object key : pairToPath.keySet()) {
            Object value = pairToPath.get(key);
            pairToPath1.put(key.toString(),
                            new Double(value.toString()));
        }
        DistanceCalculator distanceCalculator = new DistanceCalculator() {
            
            public double calculateDistance(String id1, String id2) {
                return distance(id1, id2, pairToPath1);
            }
        };
        
        HierarchicalCluster cluster = new HierarchicalCluster();
        cluster.setDistanceCalculator(distanceCalculator);
        if (useShortestPath)
            cluster.setMethod(HierarchicalCluster.ClusterDistanceMethod.SINGLE);
        else
            cluster.setMethod(HierarchicalCluster.ClusterDistanceMethod.AVERAGE);
        HierarchicalClusterNode top = cluster.cluster(ids);
        List<HierarchicalClusterNode> rtn = new ArrayList<HierarchicalClusterNode>();
        rtn.add(top);
        return rtn;
    }
    
    private double distance(String id1,
                            String id2,
                            Map<String, Double> pairToPath) {
        String key = id1 + "\t" + id2;
        Double path = pairToPath.get(key);
        if (path == null) {
            key = id2 + "\t" + id1;
            path = pairToPath.get(key);
        }
        return path;
    }
    
    private Map<String, Integer> calculateDistances(List<String> genes) throws IOException {
        Set<String> fis = getAllFIs();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        return calculateDistances(genes, bfs, geneToPartners);
    }
    
    /**
     * Calculate average distance for a list of genes.
     * @param genes
     * @param bfs
     * @param geneToPartners
     * @return
     */
    private Map<String, Integer> calculateDistances(List<String> genes,
                                                    BreadthFirstSearch bfs,
                                                    Map<String, Set<String>> geneToPartners) {
        Map<String, Integer> pairToDistance = new HashMap<String, Integer>();
        List<String> targets = new ArrayList<String>();
        for (int i = 0; i < genes.size() - 1; i++) {
            String id1 = genes.get(i);
            if (!(geneToPartners.containsKey(id1)))
                continue;
            for (int j = i + 1; j < genes.size(); j++) {
                String id2 = genes.get(j);
                if (!(geneToPartners.containsKey(id2)))
                    continue;
                targets.add(id2);
                //int length = bfs.getDistance(id1, id2, geneToPartners);
                //pairToDistance.put(id1 + "\t" + id2, length);
            }
            Map<String, Integer> targetToDist = bfs.getDistances(id1,
                                                                 targets,
                                                                 geneToPartners);
            for (String target : targetToDist.keySet()) {
                Integer dist = targetToDist.get(target);
                pairToDistance.put(id1 + "\t" + target, dist);
            }
            targets.clear();
        }
        return pairToDistance;
    }
    
    private Set<String> generateFIsForGeneSet(Set<String> fis,
                                              Set<String> spanFIs,
                                              List<String> cancerGenes) throws IOException {
        // Get all genes
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(spanFIs);
        List<String> geneList = new ArrayList<String>(genes);
        Set<String> rtn = new HashSet<String>();
        for (int i = 0; i < geneList.size() - 1; i++) {
            String gene1 = geneList.get(i);
            for (int j = i + 1; j < geneList.size(); j++) {
                String gene2 = geneList.get(j);
                // Check if this pair should have FIs
                if (isInteracting(gene1, gene2, fis))
                    rtn.add(gene1 + "\t" + gene2);
            }
        }
        return rtn;
    }
    
    private boolean isInteracting(String id1,
                                  String id2,
                                  Set<String> fis) {
        String key = id1 + "\t" + id2;
        if (fis.contains(key))
            return true;
        // Try reverse way
        key = id2 + "\t" + id1;
        return fis.contains(key);
    }
    
    private class PathSorter extends GeneSorter implements Comparator<List<String>> {
        private Set<String> linkerGenes;
        
        public PathSorter() {
        }
        
        private int getNumberOfFIInPath(List<String> path,
                                        Set<String> pathwayFIs) {
            int total = 0;
            for (int i = 0; i < path.size() - 1; i++) {
                String gene1 = path.get(i);
                String gene2 = path.get(i + 1);
                int compare = gene1.compareTo(gene2);
                String fi = null;
                if (compare < 0)
                    fi = gene1 + "\t" + gene2;
                else
                    fi = gene2 + "\t" + gene1;
                if (pathwayFIs.contains(fi))
                    total ++;
            }
            return total;
        }
        
        public int compare(List<String> path1, List<String> path2) {
            if (path1.size() < path2.size())
                return -1;
            if (path1.size() > path2.size())
                return 1;
            // Need to find a more meaningful shortest path
            // Check the total degrees
            List<String> copy1 = new ArrayList<String>(path1);
            copy1.removeAll(mutatedGenes); // Want to check linkers
            List<String> copy2 = new ArrayList<String>(path2);
            copy2.removeAll(mutatedGenes);
            Collections.sort(copy1);
            Collections.sort(copy2);
            if (copy1.equals(copy2))
                return 0; // Really don't care: same linker genes
            // Less is better
            int rtn = copy1.size() - copy2.size();
            if (rtn != 0)
                return rtn;
            // If both are zero. Fine
            if (copy1.size() == 0)
                return 0;
            int total1 = 0;
            int total2 = 0;
            // Use genes used already
            for (String gene : copy1) {
                if (linkerGenes.contains(gene))
                    total1++ ;
            }
            for (String gene : copy2) {
                if (linkerGenes.contains(gene))
                    total2 ++;
            }
            rtn = total2 - total1;
            if (rtn != 0)
                return rtn;
            total1 = total2 = 0;
            // Check the links with mutated genes
            for (String gene : copy1)
                total1 += geneToMutatedPartners.get(gene).size();
            for (String gene : copy2)
                total2 += geneToMutatedPartners.get(gene).size();
            rtn = total2 - total1;
            if (rtn != 0)
                return rtn;
            total1 = total2 = 0;
            for (String gene : copy1)
                total1 += geneToPartners.get(gene).size();
            for (String gene : copy2)
                total2 += geneToPartners.get(gene).size();
            rtn = total2 - total1;
            if (rtn != 0)
                return rtn;
            total1 = getNumberOfFIInPath(path1, pathwayFIs);
            total2 = getNumberOfFIInPath(path2, pathwayFIs);
            if (total1 > total2)
                return -1;
            if (total1 < total2)
                return 1;
            // Based on gene expression
            Double geneExp1 = 0.0d;
            Double geneExp2 = 0.0d;
            for (String gene : copy1) {
                Double tmp = geneToTValue.get(gene);
                if (tmp != null)
                    geneExp1 += tmp;
            }
            for (String gene : copy2) {
                Double tmp = geneToTValue.get(gene);
                if (tmp != null)
                    geneExp2 += tmp;
            }
            rtn = geneExp2.compareTo(geneExp1);
            if (rtn != 0)
                return rtn;
            System.out.println("Cannot see diff between two paths: " + 
                               path1 + ", " + 
                               path2);
            return 0;
        }
        
    }
    
    private class GeneSorter {
        protected Map<String, Set<String>> geneToPartners;
        protected Map<String, Set<String>> geneToMutatedPartners;
        protected Collection<String> mutatedGenes;
        protected Set<String> pathwayFIs;
        protected Map<String, Double> geneToTValue;
        
        public GeneSorter() {
        }
        
        public void createLocalParntersMap() {
            if (geneToPartners == null || mutatedGenes == null)
                throw new IllegalStateException("geneToParnters or mutatedGenes is null!");
            geneToMutatedPartners = new HashMap<String, Set<String>>();
            for (String gene : geneToPartners.keySet()) {
                Set<String> set = geneToPartners.get(gene);
                Set<String> setCopy = new HashSet<String>(set);
                setCopy.retainAll(mutatedGenes);
                geneToMutatedPartners.put(gene, setCopy);
            }
        }
    }
    
    private class EdgeSorter extends GeneSorter implements Comparator<Edge> {
        
        public EdgeSorter() {
        }
        
        public int compare(Edge edge1, Edge edge2) {
            // Find which node should be used
            // Should use the nodes that have been labeled as anchors
            TreeNode node1, node2;
            if (edge1.getNode1().getLabel() == null)
                node1 = edge1.getNode2();
            else
                node1 = edge1.getNode1();
            if (edge2.getNode1().getLabel() == null)
                node2 = edge2.getNode2();
            else
                node2 = edge2.getNode1();
            if (node1 == node2) {// They are the same
                return 0;
            }
            if (geneToMutatedPartners != null) {
                Set<String> mutatedParnters1 = geneToMutatedPartners.get(node1.getId());
                Set<String> mutatedParnters2 = geneToMutatedPartners.get(node2.getId());
                int rtn = mutatedParnters2.size() - mutatedParnters1.size();
                if (rtn != 0)
                    return rtn;
            }
            if (geneToPartners != null) {
                Set<String> partners1 = geneToPartners.get(node1.getId());
                Set<String> partners2 = geneToPartners.get(node2.getId());
                // If want to use less degree genes
                int rtn = partners1.size() - partners2.size();
                // If want to use more degree genes
                //int rtn = parnters2.size() - partners1.size();
                if (rtn != 0)
                    return rtn; // The second should be the first
            }
            // If one if in pathwayFIs and another not, use the first one
            if (pathwayFIs != null) {
                String fi1 = getFIFromEdge(edge1);
                String fi2 = getFIFromEdge(edge2);
                if (pathwayFIs.contains(fi1) && !pathwayFIs.contains(fi2))
                    return -1;
                if (!pathwayFIs.contains(fi1) && pathwayFIs.contains(fi2))
                    return 1; // The second should be the first
                if (mutatedGenes.contains(node1.getId()) && !mutatedGenes.contains(node2.getId()))
                    return -1;
                if (mutatedGenes.contains(node2.getId()) && !mutatedGenes.contains(node1.getId()))
                    return 1;
            }
            if (geneToTValue != null) {
                // Check t-value
                Double tvalue1 = geneToTValue.get(node1.getId());
                Double tvalue2 = geneToTValue.get(node2.getId());
                if (tvalue1 != null && tvalue2 != null) {
                    tvalue1 = Math.abs(tvalue1);
                    tvalue2 = Math.abs(tvalue2);
                    int rtn = tvalue2.compareTo(tvalue1);
                    if (rtn != 0)
                        return rtn;
                }
                else if (tvalue1 != null) // Pick up genes having expression values
                    return -1;
                else if (tvalue2 != null)
                    return 1;
            }
            //System.out.println(node1.getId() + ", " + node2.getId() + " cannot be ranked!");
            // This is completely random based on names
            return node1.getId().compareTo(node2.getId());
        }
        
        private String getFIFromEdge(Edge edge) {
            String gene1 = edge.getNode1().getId();
            String gene2 = edge.getNode2().getId();
            int compare = gene1.compareTo(gene2);
            if (compare < 0)
                return gene1 + "\t" + gene2;
            return gene2 + "\t" + gene1;
        }
    }
}
