/*
 * Created on Apr 27, 2007
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomDataImpl;
import org.jgrapht.Graph;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.r3.graph.BreadthFirstSearch.Edge;
import org.reactome.r3.graph.BreadthFirstSearch.TreeNode;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

public class ShortestPathAnalyzer {
    //public static final String BIG_COMP_INT_FILE = "FI73_042108_BigComp.txt";
    public static final String BIG_COMP_INT_FILE = R3Constants.INTERACTION_BIG_COMP_FILE_NAME;
    private GraphAnalyzer graphAnalyzer;
    
    public ShortestPathAnalyzer() {
        graphAnalyzer = new GraphAnalyzer();
    }
    
    /**
     * Calculate shortest pathway for a collection of genes.
     * @param genes
     * @return
     */
    public double calculateShortestPath(Collection<String> genes) throws IOException {
        // Check average shortest path
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        FileUtility fu = new FileUtility();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<TreeNode, List<Edge>> geneToPartners = bfs.initGraph(fis);
        List<String> geneList = new ArrayList<String>(genes);
        geneList.retainAll(InteractionUtilities.grepIDsFromInteractions(fis));
        double averagePath = calculateShortestPath(geneList, 
                                                   bfs, 
                                                   geneToPartners);
        return averagePath;
    }
    
    public double calculateShortestPath(Collection<String> ids,
                                         BreadthFirstSearch bfs,
                                         Map<TreeNode, List<Edge>> nodeToEdges) {
        Map<String, List<String>> pairToPath = bfs.generateShortestPath(new ArrayList<String>(ids),
                                                                        nodeToEdges);
        int total = 0;
        for (String pair : pairToPath.keySet()) {
            total += pairToPath.get(pair).size() - 1;
        }
        double value = (double) total / pairToPath.size();
        return value;
    }
    
//    /**
//     * This method is used to check average shortest path for genes in chromosome
//     * to see if there are any correlation.
//     * @throws IOException
//     */
//    @Test
//    public void checkAverageShortestPathInChromosome() throws IOException {
//        UCSCDataAnalyzer ucscAnalyzer = new UCSCDataAnalyzer();
//        List<RefSeqInfo> refSeqInfos = ucscAnalyzer.loadHumanRefSeqInfos();
//        Map<String, List<RefSeqInfo>> chrToRefSeqs = ucscAnalyzer.convertRefSeqListToMap(refSeqInfos);
//        Map<String, Integer> chrToLength = ucscAnalyzer.loadChromosomeSizes();
//        ucscAnalyzer.sortRefSeqsInChr(chrToRefSeqs);
//        Map<String, List<String>> chrToGeneNames = ucscAnalyzer.convertRefSeqToGeneInOrder(chrToRefSeqs);
//        Map<String, Object[]> genePosMap = generateGenePosMap(chrToGeneNames);
//        
//        List<String> chrList = new ArrayList<String>(chrToGeneNames.keySet());
//        Collections.sort(chrList);
//        int totalGenes = 0;
//        for (String chr : chrList) {
//            if (chr.endsWith("_random") || chr.contains("hap"))
//                continue;
//            List<String> genes = chrToGeneNames.get(chr);
//            System.out.println(chr + "\t" + genes.size());
//            totalGenes += genes.size();
//        }
//        System.out.println("Total genes: " + totalGenes);
//        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
//        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        //BreadthFirstSearch bfs = new BreadthFirstSearch();
//        //Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(fis);
//        // Check correlation
//        // Check if these genes are adjacent in the chromosome
//        int count = 0;
//        for (String fi : fis) {
//            int index = fi.indexOf("\t");
//            String gene1 = fi.substring(0, index);
//            String gene2 = fi.substring(index + 1);
//            if(isChromosomeAdjacent(genePosMap, 
//                                    gene1,
//                                    gene2))
//                count ++;
//        }
//        System.out.println();
//        System.out.println("FIs: " + count + " " + (double)count / fis.size());
//        List<String> geneList = new ArrayList<String>(fiGenes);
//        Random random = new Random();
//        List<Double> randomCountList = new ArrayList<Double>();
//        int permutation = 1000;
//        for (int j = 0; j < permutation; j++) {
//            int randomCount = 0;
//            for (int i = 0; i < fis.size(); i++) {
//                int index1 = random.nextInt(geneList.size());
//                String gene1 = geneList.get(index1);
//                int index2 = random.nextInt(geneList.size());
//                String gene2 = geneList.get(index2);
//                if(isChromosomeAdjacent(genePosMap, 
//                                        gene1,
//                                        gene2))
//                    randomCount ++;
//            }
//            System.out.println("Random count: " + randomCount + " " + (double)randomCount / fis.size());
//            randomCountList.add((double)randomCount/fis.size());
//        }
//        Collections.sort(randomCountList);
//        for (Double v : randomCountList)
//            System.out.println(v);
////        System.out.println("\nCorrelation based on ranks:");
////        for (String chr : chrList) {
////            if (chr.endsWith("_random") || chr.contains("hap"))
////                continue;
////            List<String> genes = chrToGeneNames.get(chr);
////            genes.retainAll(fiGenes);
////            List<Double> chrDistList = new ArrayList<Double>();
////            List<Double> fiDistList = new ArrayList<Double>();
////            for (int i = 0; i < permutation; i++) {
////                int index1 = random.nextInt(genes.size());
////                int index2 = random.nextInt(genes.size());
////                String gene1 = genes.get(index1);
////                String gene2 = genes.get(index2);
////                int chrDist = Math.abs(index1 - index2);
////                chrDistList.add(new Double(chrDist));
////                int fiDist = bfs.getDistance(gene1, gene2, idToPartners);
////                fiDistList.add(new Double(fiDist));
////            }
////            double cor = MathUtilities.calculatePearsonCorrelation(chrDistList, fiDistList);
////            System.out.println(chr + "\t" + cor);
////        }
////        System.out.println("\nCorrelation based on positions:");
////        for (String chr : chrList) {
////            if (chr.endsWith("_random") || chr.contains("hap"))
////                continue;
////            List<RefSeqInfo> refSeqList = chrToRefSeqs.get(chr);
////            // Clean up a little bit
////            Set<String> refGenes = new HashSet<String>();
////            List<RefSeqInfo> filtered = new ArrayList<RefSeqInfo>();
////            for (RefSeqInfo seq : refSeqList) {
////                if (refGenes.contains(seq.getGeneName()))
////                    continue;
////                if (!fiGenes.contains(seq.getGeneName()))
////                    continue;
////                filtered.add(seq);
////            }
////            List<Double> chrDistList = new ArrayList<Double>();
////            List<Double> fiDistList = new ArrayList<Double>();
////            int length = chrToLength.get(chr);
////            for (int i = 0; i < permutation; i++) {
////                int index1 = random.nextInt(filtered.size());
////                int index2 = random.nextInt(filtered.size());
////                RefSeqInfo refseq1 = filtered.get(index1);
////                RefSeqInfo refseq2 = filtered.get(index2);
////                int chrDist = Math.abs(refseq1.getStart() - refseq2.getStart());
////                chrDistList.add((double) chrDist / length);
////                int fiDist = bfs.getDistance(refseq1.getGeneName(),
////                                             refseq2.getGeneName(),
////                                             idToPartners);
////                fiDistList.add(new Double(fiDist));
////            }
////            double cor = MathUtilities.calculatePearsonCorrelation(chrDistList, fiDistList);
////            System.out.println(chr + "\t" + cor);
////        }
//    }

    private boolean isChromosomeAdjacent(Map<String, Object[]> geneToPosInfo,
                                         String gene1,
                                         String gene2) {
        Object[] info1 = geneToPosInfo.get(gene1);
        Object[] info2 = geneToPosInfo.get(gene2);
        if (info1 == null || info2 == null)
            return false;
        if (info1[0].equals(info2[0])) { // Make sure they are in the same chromosome
            Integer index1 = (Integer) info1[1];
            Integer index2 = (Integer) info2[1];
            if (Math.abs(index1 - index2) < 2)
                return true;
        }
        return false;
    }
    
    private Map<String, Object[]> generateGenePosMap(Map<String, List<String>> chrToGeneNames) {
        Map<String, Object[]> rtn = new HashMap<String, Object[]>();
        for (String chr : chrToGeneNames.keySet()) {
            List<String> geneNames = chrToGeneNames.get(chr);
            for (int i = 0; i < geneNames.size(); i++) {
                String gene = geneNames.get(i);
                Object[] tmp = new Object[2];
                tmp[0] = chr;
                tmp[1] = i;
                rtn.put(gene, tmp);
            }
        }
        return rtn;
    }
    
//    /**
//     * This method is used to split the random generated shortest path files
//     * into two parts: proteins pairs are in the same pathways, and proteins
//     * pairs are not in the same pathways.
//     */
//    @Test
//    public void splitShortestPathFile() throws IOException {
//        String fileName = R3Constants.RESULT_DIR + "ShortestPathFromRandomPairs.txt";
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        List<String> inPathwayLines = new ArrayList<String>();
//        List<String> nonInPathwayLines = new ArrayList<String>();
//        Map<String, Set<String>> topicToIds = new TopicAnalyzer().getTopicToIdMap();
//        String line = null;
//        boolean inPathway;
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("[ \t]");
//            // Check if tokens[0] and tokens[1] are in the same pathways
//            inPathway = false;
//            for (Iterator<String> it = topicToIds.keySet().iterator(); it.hasNext();) {
//                String topic = it.next();
//                Set<String> idSet = topicToIds.get(topic);
//                if (idSet.contains(tokens[0]) && idSet.contains(tokens[1])) {
//                    inPathway = true;
//                    break;
//                }
//            }
//            if (inPathway)
//                inPathwayLines.add(line);
//            else
//                nonInPathwayLines.add(line);
//        }
//        System.out.println("In pathways: " + inPathwayLines.size());
//        System.out.println("Not in Pathways: " + nonInPathwayLines.size());
//        fileName = R3Constants.RESULT_DIR + "ShortestPathFromRandomPairsInPathways.txt";
//        fu.setOutput(fileName);
//        for (String line1 : inPathwayLines)
//            fu.printLine(line1);
//        fu.close();
//        fileName = R3Constants.RESULT_DIR + "ShortestPathFromRandomPairsNonInPathways.txt";
//        fu.setOutput(fileName);
//        for (String line1 : nonInPathwayLines)
//            fu.printLine(line1);
//        fu.close();
//    }
//    
    /**
     * Remove lines that have path length is -1.
     * @throws IOException
     */
    @Test
    public void cleanUpShortestPath() throws IOException {
        String fileName = R3Constants.RESULT_DIR + "ShortestPathFromRandomPairs.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // Hold all data
        List<String> lines = new ArrayList<String>();
        String line = null;
        int index = 0;
        int path = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            path = Integer.parseInt(line.substring(index + 1));
            if (path < 1) // self or none
                continue;
            lines.add(line);
        }
        fu.close();
        fu.setOutput(fileName);
        for (String tmp : lines) {
            fu.printLine(tmp);
        }
        fu.close();
    }
    
//    public void analyzeShortestPath() throws Exception {
//        // Data from Reactome
//        String reactomeFileName = R3Constants.RESULT_DIR + "ReactomeInteractions.txt";
//        FileUtility fu = new FileUtility();
//        Set<String> reactomeInteractions = fu.loadInteractions(reactomeFileName);
//        Set<String> reactomeIds = InteractionUtilities.grepIDsFromInteractions(reactomeInteractions);
//        Graph<String, DefaultEdge> reactomeGraph = graphAnalyzer.createGraph(reactomeIds, reactomeInteractions);
//        // Data from PPI
//        HumanPsiMiInteractionAnalyzer interactionAnalyzer = new HumanPsiMiInteractionAnalyzer();
//        //Set<String> ppiInteractions = interactionAnalyzer.loadHumanInteractions();
//        Set<String> ppiInteractions = interactionAnalyzer.loadNonY2HInteraction();
//        Set<String> ppiIds = InteractionUtilities.grepIDsFromInteractions(ppiInteractions);
//        Graph<String, DefaultEdge> ppiGrah = graphAnalyzer.createGraph(ppiIds, ppiInteractions);
//        // Find overlapping ids
//        Set<String> overlappingIds = new HashSet<String>(reactomeIds);
//        overlappingIds.retainAll(ppiIds);
//        System.out.println("Total overlapping: " + overlappingIds.size());
//        // Check the first 300
//        int c = 0;
//        List<String> idList = new ArrayList<String>();
//        for (String id : overlappingIds) {
//            idList.add(id);
//            c ++;
//            //if (c > 300)
//              //  break;
//        }
//        fu.setOutput(R3Constants.RESULT_DIR + "PathComp.txt");
//        fu.printLine("Reactome PPI");
//        c = 0;
//        for (int i = 0; i < idList.size() - 1; i++) {
//            String id1 = idList.get(i);
//            for (int j = i + 1; j < idList.size(); j++) {
//                String id2 = idList.get(j);
//                List<DefaultEdge> reactomePath = DijkstraShortestPath.findPathBetween(reactomeGraph, id1, id2);
//                if (reactomePath == null)
//                    continue;
//                List<DefaultEdge> ppiPath = DijkstraShortestPath.findPathBetween(ppiGrah, id1, id2);
//                if (ppiPath == null)
//                    continue;
//                fu.printLine(reactomePath.size() + " " + ppiPath.size());
//                c ++;
//                //if (c % 10 == 0)
//                  //  System.out.println(c);
//            }
//        }
//        System.out.println("Total Paths: " + c);
//        fu.close();
//    }    
    
    /**
     * Jung is much slower than JGraphT. Probably it should not be used. 
     * @throws IOException
     */
//    @Test
//    public void generateShortestPathsMatrixByJung() throws IOException {
//        FileUtility fu = new FileUtility();
//        String fileName = GRAPH_DIR  + "PathOfTenOfRandomAndTenOfBCR(N)FromJung.txt";
//        fu.setOutput(fileName);
//        Map<String, Vertex> idToVertex = new HashMap<String, Vertex>();
//        edu.uci.ics.jung.graph.Graph graph = createJungGraph(R3Constants.INTERACTION_FILE_NAME, idToVertex);
//        Set<String> pathwaySet = getProteinIds();
//        //Set<String> randomSet = getRandomProteinIds();
//        Set<String> randomSet = new HashSet<String>();
//        Set<String> proteinIdSet = new HashSet<String>();
//        proteinIdSet.addAll(pathwaySet);
//        proteinIdSet.addAll(randomSet);
//        fu.printLine("Test Protein Ids: " + proteinIdSet);
//        // Do a pairwise shortest path calculation
//        List<String> proteinIds = new ArrayList<String>(proteinIdSet);
//        Collections.sort(proteinIds);
//        int size = proteinIds.size();
//        // Check time
//        long time1 = System.currentTimeMillis();
//        UnweightedShortestPath pathEngine = new UnweightedShortestPath(graph);
//        List<Integer> paths = new ArrayList<Integer>();
//        for (int i = 0; i < size - 1; i++) {
//            String id1 = proteinIds.get(i);
//            Vertex v1 = idToVertex.get(id1);
//            for (int j = i + 1; j < size; j++) {
//                String id2 = proteinIds.get(j);
//                Vertex v2 = idToVertex.get(id2);
//                Number dist = pathEngine.getDistance(v1, v2);
//                if (dist != null)
//                    paths.add(dist.intValue());
//            }
//        }
//        System.out.println("Values in paths: " + paths.size());
//        System.out.println(paths);
//        fu.close();
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for ten x ten shortest path: " + (time2 - time1));
//    }
    
//  private edu.uci.ics.jung.graph.Graph createJungGraph(String intFileName,
//  Map<String, Vertex> idToVertex) throws IOException {
//FileUtility fu = new FileUtility();
//Set<String> interactions = fu.loadInteractions(intFileName);
//Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
//edu.uci.ics.jung.graph.Graph graph = createJungGraph(ids, interactions, idToVertex);
//return graph;
//}
//
///**
//* @param ids
//* @param interactions
//* @param idToVertex
//* @return
//*/
//private edu.uci.ics.jung.graph.Graph createJungGraph(Set<String> ids,
//  Set<String> interactions,
//  Map<String, Vertex> idToVertex) {
//int index = 0;
//// Just want to use String as edge class
//UndirectedSparseGraph graph = new UndirectedSparseGraph();
//for (String id : ids) {
//Vertex v = new SimpleUndirectedSparseVertex();
//graph.addVertex(v);
//idToVertex.put(id, v);
//}
//for (String pair : interactions) {
//index = pair.indexOf(" ");
//Vertex v1 = idToVertex.get(pair.substring(0, index));
//Vertex v2 = idToVertex.get(pair.substring(index + 1));
//UndirectedSparseEdge edge = new UndirectedSparseEdge(v1, v2);
//graph.addEdge(edge);
//}
//System.out.printf("Graph: vertics %d edges %d%n", 
//graph.numVertices(),
//graph.numEdges());
//return graph;
//}
//    
//    @Test
//    public void generateShortestPathMatrix() throws IOException {
//        FileUtility fu = new FileUtility();
//        String fileName = GraphAnalyzer.GRAPH_DIR  + "FiftyOfRandom.txt";
//        fu.setOutput(fileName);
//        Graph<String, DefaultEdge> graph = graphAnalyzer.createGraph(R3Constants.INTERACTION_FILE_NAME);
//        Set<String> pathwaySet = getProteinIds();
//        Set<String> randomSet = getRandomProteinIds();
//        Set<String> proteinIdSet = new HashSet<String>();
//        proteinIdSet.addAll(randomSet);
//        //proteinIdSet.addAll(pathwaySet);
//        fu.printLine("Test Protein Ids: " + proteinIdSet);
//        // Do a pairwise shortest path calculation
//        List<String> proteinIds = new ArrayList<String>(proteinIdSet);
//        Collections.sort(proteinIds);
//        int size = proteinIds.size();
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        Map<String, Set<String>> fiToTopics = topicAnalyzer.loadFIToTopics();
//        // Check time
//        long time1 = System.currentTimeMillis();
//        for (int i = 0; i < size - 1; i++) {
//            String id1 = proteinIds.get(i);
//            for (int j = i + 1; j < size; j++) {
//                String id2 = proteinIds.get(j);
//                List<DefaultEdge> path = DijkstraShortestPath.findPathBetween(graph, 
//                                                                               id1, 
//                                                                               id2);
//                fu.printLine(id1 + "->" + id2 + ":");
//                fu.printLine("\t" + path);
//                fu.printLine("\t" + annotatePath(path, fiToTopics, fu));
////                System.out.println(id1 + "->" + id2 + ":");
////                System.out.println("\t" + path);
////                System.out.println("\t" + annotatePath(path, fiToTopics, fu));
//            }
//        }
//        long time2 = System.currentTimeMillis();
//        //System.out.println("Time for ten x ten shortest path: " + (time2 - time1));
//        fu.printLine("Time (s): " + (time2 - time1) / 1000);
//        fu.close();
//    }
//    
    @Test
    public void calculateAverageShortestPathForRandomGraph() throws IOException {
        FileUtility fu = new FileUtility();
        //String intFileName = R3Constants.INTERACTION_FILE_NAME;
        //String intFileName = R3Constants.RESULT_DIR + "FIInteractions60_121707_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + "FIInteractions73_021308_BigComp.txt";
        String intFileName = R3Constants.RESULT_DIR + BIG_COMP_INT_FILE;
        Set<String> interactions = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        //String outputFileName = "AveragePathForRandomGraph_1000Sampl.txt";
        //String outputFileName = "AveragePathForRandomGraph73_021308_1000Sampl.txt";
        String outputFileName = "AveragePathForRandomGraphBasedOn" + BIG_COMP_INT_FILE;
        fu.setOutput(R3Constants.RESULT_DIR + outputFileName);
        for (int i = 0; i < 6; i++) {
            fu.printLine("Randome graph sampling " + i + "...");
            Graph<String, DefaultEdge> graph = graphAnalyzer.createRandomGraph(ids, interactions.size());
            fu.printLine("Graph: " + graph.vertexSet().size() + " vertices, " + 
                               graph.edgeSet().size() + " edges");
            // Need to get the interctions
            Set<String> newInteractions = new HashSet<String>();
            int index;
            for (DefaultEdge edge : graph.edgeSet()) {
                //edge is output as: (Q6DT37 : Q1WIR0)
                String text = edge.toString();
                index = text.indexOf(":");
                String id1 = text.substring(1, index - 1);
                String id2 = text.substring(index + 2, text.length() - 1);
                int compare = id1.compareTo(id2);
                if (compare < 0)
                    newInteractions.add(id1 + " " + id2);
                else
                    newInteractions.add(id2 + " " + id1);
            }
            calculateAverageShortestPathFromRandomPairs(fu, 
                                                        ids,
                                                        newInteractions,
                                                        1000,
                                                        1);
        }
        fu.close();
    }
    
    /**
     * This method is used to generate a shortest path values as a distributions.
     * @throws Exception
     */
    @Test
    public void generateShortesPathDistribution() throws Exception {
        FileUtility fu = new FileUtility();
        String intFileName = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
        Set<String> interactions = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        //bfs.setInteractions(interactions);
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(interactions);
        int totalSample = 10000;
        Set<String> checkedPairs = new HashSet<String>();
        List<String> idList = new ArrayList<String>(ids);
        int totalSize = idList.size();
        int index1, index2;
        String key = null;
        String id1, id2;
        List<Integer> pathList = new ArrayList<Integer>();
        int path = 0;
        for (int i = 0; i < totalSample; i++) {
            index1 = (int) (Math.random() * totalSize);
            index2 = (int) (Math.random() * totalSize);
            if (index1 == index2)
                continue;
            if (index1 < index2)
                key = index1 + " " + index2;
            else
                key = index2 + " " + index1;
            if (checkedPairs.contains(key))
                continue;
            checkedPairs.add(key);
            id1 = idList.get(index1);
            id2 = idList.get(index2);
            path = bfs.getDistance(id1, id2, idToPartners);
            pathList.add(path);
        }
        System.out.println("Path: " + pathList.size());
        Collections.sort(pathList);
        String outputFileName = R3Constants.RESULT_DIR + "ShortestPathFromFIsInGene_041709_BigComp.txt";
        fu.setOutput(outputFileName);
        for (Integer p : pathList)
            fu.printLine(p + "");
        fu.close();
    }
    
    @Test
    public void calculateAverageShortestPath() throws IOException {
        FileUtility fu = new FileUtility();
        //String intFileName = R3Constants.INTERACTION_FILE_NAME;
        //String intFileName = R3Constants.RESULT_DIR + "FIInteractions60_121707_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + "FIInteractions73_021308_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + BIG_COMP_INT_FILE ;
//        String intFileName = R3Constants.RESULT_DIR + "FI73_042108_BigComp.txt";
//        String intFileName = R3Constants.RESULT_DIR + "FIs_042109_BigComp.txt";
        String intFileName = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
        Set<String> interactions = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
//        Graph<String, DefaultEdge> graph = graphAnalyzer.createGraph(intFileName);
//        System.out.println("Graph: " + graph.vertexSet().size() + " vertices, " + 
//                           graph.edgeSet().size() + " edges");
        //String outputFileName = "AveragePathForFIInteractions60_121707_BigComp_1000Sampl.txt";
        //String outputFileName = "AveragePathForFIInteractions73_021308_BigComp_1000Sampl.txt";
        //String outputFileName = "AveragePathFrom1000SamplFor" + BIG_COMP_INT_FILE;
        //String outputFileName = "AveragePathFrom1000SamplForFI73InGene_061008_BigComp.txt";
        //String outputFileName = "AveragePathFrom1000SamplForFI73_042108_BigComp.txt";
        //String outputFileName = "AveragePathFrom1000SamplForFIs_042109_BigComp.txt";
        String outputFileName = "AveragePathFrom1000SamplForFIs_In_GENE_BigComp.txt";
        fu.setOutput(R3Constants.RESULT_DIR + outputFileName);
        calculateAverageShortestPathFromRandomPairs(fu, 
                                                    ids,
                                                    interactions,
                                                    1000,
                                                    6);
        fu.close();
//        calculateAverageShortestPath(fu,
//                                     interactions, 
//                                     "AveragePathForFIInteractions60_121707_BigComp.txt");
    }
    
    private void calculateAverageShortestPathFromRandomPairs(FileUtility fu,
                                                             Set<String> ids,
                                                             Set<String> interactions,
                                                             int pairSampleNumber,
                                                             int totalSample) throws IOException {
        double totalAverage = 0;
        // Use list for sampling
        List<String> idList = new ArrayList<String>(ids);
        Collections.sort(idList);
        Set<String> sampledPairs = new HashSet<String>();
        int index1, index2;
        String id1, id2;
        String key;
        int totalSize = idList.size();
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        //bfs.setInteractions(interactions);
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(interactions);
        for (int i = 0; i < totalSample; i++) {
            System.out.println("Sampling " + i + "...");
            long time1 = System.currentTimeMillis();
            int totalPath = 0;
            // Get the total pair first
            sampledPairs.clear();
            while (sampledPairs.size() < pairSampleNumber) {
                index1 = (int) (Math.random() * totalSize);
                index2 = (int) (Math.random() * totalSize);
                if (index1 == index2)
                    continue;
                if (index1 < index2)
                    key = index1 + " " + index2;
                else
                    key = index2 + " " + index1;
                if (sampledPairs.contains(key))
                    continue;
                sampledPairs.add(key);
                id1 = idList.get(index1);
                id2 = idList.get(index2);
                totalPath += bfs.getDistance(id1, id2, idToPartners);
            }
            double average = (double) totalPath / sampledPairs.size();
            long time2 = System.currentTimeMillis();
            System.out.println("Time: " + (time2 - time1));
            System.out.println("Sample " + i + ": " + average);
            fu.printLine("Sample " + i + ": " + average);
            totalAverage += average;
        }
        fu.printLine("Average: " + totalAverage / totalSample);
    }
    
    private void calculateAverageShortestPath(FileUtility fu,
                                               Set<String> interactions,
                                               String outputFileName) throws IOException {
        FloydWarshall fw = new FloydWarshall();
        fw.calculateAllPairs(interactions);
        int[][] dists = fw.getAllPairDistances();
        int size = dists.length; 
        int total = 0;
        int count = 0;
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                if (dists[i][j] < Integer.MAX_VALUE) {
                    total = dists[i][j];
                    count ++;
                }
            }
        }
        fu.setOutput(R3Constants.RESULT_DIR + outputFileName);
        fu.printLine("Total count: " + count);
        fu.printLine("Total path length: " + total);
        fu.printLine("Average path: " + (double)total/count);
        fu.close();
    }

    private void calculateAverageShortestPath(FileUtility fu,
                                              Set<String> ids, 
                                              Graph<String, DefaultEdge> graph,
                                              String outputFileName) throws IOException {
        List<String> idList = new ArrayList<String>(ids);
        int total = 0;
        int count = 0;
        for (int i = 0; i < idList.size() - 1; i++) {
            String id1 = idList.get(i);
            for (int j = i + 1; j < idList.size(); j++) {
                String id2 = idList.get(j);
                List<DefaultEdge> path = DijkstraShortestPath.findPathBetween(graph, id1, id2);
                if (path == null)
                    continue;
                count ++;
                total += path.size();
            }
        }
        fu.setOutput(R3Constants.RESULT_DIR + outputFileName);
        fu.printLine("Total count: " + count);
        fu.printLine("Total path length: " + total);
        fu.printLine("Average path: " + (double)total/count);
        fu.close();
    }
    
    private Set<String> annotatePath(List<DefaultEdge> path,
                                     Map<String, Set<String>> fiToTopics,
                                     FileUtility fu) throws IOException {
        Set<String> topics = new HashSet<String>();
        if (path == null)
            return topics;
        for (DefaultEdge edge : path) {
            String edgeName = edge.toString();
            String fi = graphAnalyzer.convertEdgeNameToFI(edgeName);
            Set<String> fiTopics = fiToTopics.get(fi);
            if (fiTopics == null) {
                fu.printLine("Cannot find annotations: " + fi);
                continue;
            }
            topics.addAll(fiTopics);
        }
        return topics;
    }
//    
//    private Set<String> getProteinIds() throws IOException {
//        TopicAnalyzer analyzer = new TopicAnalyzer();
//        Map<String, Set<String>> topicToIds = analyzer.getTopicToIdMap();
//        // Try Reactome Apoptosis 
//        //String topic = "Apoptosis(R)";
//        //String topic = "Integrin signalling pathway(P)";
//        //String topic = "Wnt signaling pathway(P)";
//        //String topic = "HIV Infection(R)";
//        //String topic = "EGFR1(C)";
//        //String topic = "TGFBR(C)";
//        String topic = "BCR Signaling Pathway(N)";
//        Set<String> ids = topicToIds.get(topic);
//        RandomData randomData = new RandomDataImpl();
//        Object[] tenIds = randomData.nextSample(ids, 10);
//        Set<String> rtn = new HashSet<String>();
//        for (Object obj : tenIds)
//            rtn.add(obj.toString());
//        return rtn;
//    }
    
    private Set<String> getRandomProteinIds() throws IOException {
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(R3Constants.INTERACTION_FILE_NAME);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        Object[] tenIds = new RandomDataImpl().nextSample(ids, 50);
        Set<String> rtn = new HashSet<String>();
        for (Object obj : tenIds)
            rtn.add(obj.toString());
        return rtn;
    }
//    
//    @Test
//    public void generateShortestPathFromRandomPairs() throws IOException {
//        System.out.println("Starting running...");
//        long time1 = System.currentTimeMillis();
//        int total = 20000;
//        // For comparision
//        FIFileAnalyzer fiAnalyzer = new FIFileAnalyzer();
//        Map<String, Set<String>> idToInteractions = fiAnalyzer.loadIdToPartners();
//        List<String> ids = new ArrayList<String>(idToInteractions.keySet());
//        FileUtility fu = new FileUtility();
//        String fileName = R3Constants.RESULT_DIR + "ShortestPathFromRandomPairs.txt";
//        fu.setOutput(fileName);
//        // calculate in total
//        int c = 0;
//        int index = 0;
//        int listSize = ids.size();
//        String id1 = null;
//        String id2 = null;
//        int pathLength = 0;
//        int comp = 0;
//        while (c < total) {
//            index = (int) (Math.random() * listSize);
//            id1 = ids.get(index);
//            index = (int) (Math.random() * listSize);
//            id2 = ids.get(index);
//            pathLength = calculateShortestPath(id1, id2, idToInteractions);
//            if (pathLength < 1)
//                continue; // Escape it
//            comp = id1.compareTo(id2);
//            if (comp < 0)
//                //System.out.println(id1 + " " + id2 + "\t" + pathLength);
//                fu.printLine(id1 + " " + id2 + "\t" + pathLength);
//            else
//                //System.out.println(id2 + " " + id1 + "\t" + pathLength);
//                fu.printLine(id2 + " " + id1 + "\t" + pathLength);
//            c ++;
//        }
//        fu.close();
//        long time2 = System.currentTimeMillis();
//        System.out.println("Done: " + (time2 - time1));
//    }
//    
    private int calculateShortestPath(String source, 
                                      String target, 
                                      Map<String, Set<String>> idToPartners) {
        int length = 0;
        if (source.equals(target))
            return length; // Just in case
        Set<String> checked = new HashSet<String>();
        Set<String> next = new HashSet<String>();
        Set<String> current = new HashSet<String>();
        current.add(source);
        boolean isFound = false;
        while (current.size() > 0) {
            length ++;
            for (String checking : current) {
                Set<String> partners = idToPartners.get(checking);
                if (partners.contains(target)) {
                    isFound = true;
                    break;
                }
                next.addAll(partners);
                checked.add(checking);
            }
            if (isFound)
                break;
            current.clear();
            current.addAll(next);
            // Remove ids that have been checked already. Otherwise,
            // an infinity will occur.
            current.removeAll(checked);
            next.clear();
        }
        if (!isFound)
            length = -1;
        return length;
    }

}
