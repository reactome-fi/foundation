/*
 * Created on Jun 28, 2006
 *
 */
package org.reactome.r3.graph;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.alg.BronKerboschCliqueFinder;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.junit.Test;
import org.reactome.r3.graph.ga.Chromosome;
import org.reactome.r3.graph.ga.FastTSPSolverViaGA;
import org.reactome.r3.graph.ga.Gene;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to do graph related analysis. It is based on JGraphT library.
 * @author guanming
 *
 */
public class GraphAnalyzer {
    //private final String RESULT_DIR = "results/interaction/";
    //private final String RESULT_DIR = "results/v2/mcl/";
    private final String RESULT_DIR = R3Constants.RESULT_DIR;
    public static final String GRAPH_DIR = R3Constants.RESULT_DIR + "graph/";
    
    public GraphAnalyzer() {
    }
    
    public void getHubNonHubProteins(List<String> hubProteins, 
                                     List<String> nonHubProteins,
                                     String intFileName) throws IOException {
        //String intFileName = "results/v2/FI73_042108.txt";
        // Get hub and non-hub proteins based on the connection degree
        final Map<String, Integer> proteinToDegree = generateProteinToDegree(intFileName);
        getHubNonHubProteins(hubProteins, 
                             nonHubProteins, 
                             proteinToDegree,
                             0.05d);
    }

    public void getHubNonHubProteins(List<String> hubProteins,
                                     List<String> nonHubProteins,
                                     final Map<String, Integer> proteinToDegree,
                                     double percentile) {
        hubProteins.clear();
        nonHubProteins.clear();
        List<String> proteinList = new ArrayList<String>(proteinToDegree.keySet());
        // Sort proteins based on connection degrees
        Collections.sort(proteinList, new Comparator<String>() {
           public int compare(String protein1, String protein2) { 
               Integer degree1 = proteinToDegree.get(protein1);
               Integer degree2 = proteinToDegree.get(protein2);
               return degree2 - degree1;
           }
        });
        // Pick the first 5% as the hub proteins
        int hubSize = (int) (proteinList.size() * percentile);
        int nonHubSize = proteinList.size() - hubSize;
        for (int i = 0; i < proteinList.size(); i ++) {
            String protein = proteinList.get(i);
//            if (protein.matches("^ZNF(\\d)+")) {
//                System.out.println("Zinc finger proteins: " + protein);
//                nonHubProteins.add(protein);
//                continue;
//            }
            if (i < hubSize)
                hubProteins.add(protein);
            else
                nonHubProteins.add(protein);
            //System.out.println(protein + ": " + proteinToDegree.get(protein));
        }
    }
    
    @Test
    public void removeTopOnePercentHubs() throws IOException {
        String fileName = RESULT_DIR + "FI73InGene_102908.txt";
        final Map<String, Integer> proteinToDegree = generateProteinToDegree(fileName);
        // Sort the file
        List<String> proteins = new ArrayList<String>(proteinToDegree.keySet());
        Collections.sort(proteins, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                int degree1 = proteinToDegree.get(gene1);
                int degree2 = proteinToDegree.get(gene2);
                return degree2 - degree1;
            }
        });
        int cutoff = (int)(proteins.size() * 0.01);
        System.out.println("Remove proteins: " + cutoff);
        for (int i = 0; i < cutoff; i++)
            proteins.remove(0);
        FileUtility fu = new FileUtility();
        Set<String> fis = fu.loadInteractions(fileName);
        String output = RESULT_DIR + "FI73InGene_Remove_Hubs_110308.txt";
        fu.setOutput(output);
        for (String fi : fis) {
            int index = fi.indexOf("\t");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            if (proteins.contains(id1) && proteins.contains(id2)) {
                fu.printLine(fi);
            }
        }
        fu.close();
    }
    
    @Test
    public void calculateAverageDegree() throws IOException {
//        String fileName = R3Constants.INTERACTION_FILE_NAME;
        //String fileName = R3Constants.RESULT_DIR + "FI73InGene_102908_BigComp.txt";
        //String fileName = "results/v2/FI73InGene_102908_BigComp.txt";
        //String fileName = R3Constants.RESULT_DIR + "FIs_042109_BigComp.txt";
//        String fileName = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
//        String fileName = R3Constants.RESULT_DIR + "FIsInGene_071012_BigComp.txt";
        String fileName = "results/v3/FIsInGene_041709_BigComp.txt";
        Map<String, Integer> proteinToDegree = generateProteinToDegree(fileName);
        int total = 0;
        int highest = 0;
        String highestProtein = null;
        for (Iterator<String> it = proteinToDegree.keySet().iterator(); it.hasNext();) {
            String protein = it.next();
            Integer degree = proteinToDegree.get(protein);
            total += degree;
            if (degree > highest) {
                highest = degree;
                highestProtein = protein;
            }
        }
        System.out.println("Average degree: " + total + "/" + proteinToDegree.size() + " = " + ((double)total / proteinToDegree.size()));
        System.out.println("Highest degree: " + highest + " for protein " + highestProtein);
    }
    
    public Map<String, Integer> generateProteinToDegree(String fileName) throws IOException {
        //String fileName = "/Users/wgm/Documents/xin/ReactomeInteractions_GENEID.txt";
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(fileName);
        return InteractionUtilities.generateProteinToDegree(interactions);
    }
    
    /**
     * This method is used to create a map from protein to connection degrees, and generate
     * interaction degree frequency file to calcualte scale-free property for the FI network.
     * @throws IOException
     */
    @Test
    public void countInteractionNumbersForProteins() throws IOException {
        String fileName = R3Constants.INTERACTION_FILE_NAME;
        //String fileName = R3Constants.RESULT_DIR + "FI73InGene_102908_BigComp.txt";
        fileName = R3Constants.RESULT_DIR + "FIsInGene_071012.txt";
        final Map<String, Integer> proteinToInNumber = generateProteinToDegree(fileName);
        FileUtility fu = new FileUtility();
        //fu.setOutput("results/v2/ProtToIntNum121707.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "ProtToIntNum042108_Predicated.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "GeneToIntNum122107.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "GeneToIntNum111308.txt");
//        fu.setOutput(R3Constants.RESULT_DIR + "ProteinToInNum072809.txt");
        fu.setOutput(R3Constants.RESULT_DIR + "GeneToIntNumber091912.txt");
        // Create a header
        fu.printLine("Protein\tInteractions");
        List<String> ids = new ArrayList<String>(proteinToInNumber.keySet());
        Collections.sort(ids, new Comparator<String>() {
            public int compare(String id1, String id2) {
                Integer c1 = proteinToInNumber.get(id1);
                Integer c2 = proteinToInNumber.get(id2);
                return -c1.compareTo(c2);
            }
        });
        for (String id : ids) {
            fu.printLine(id + "\t" + proteinToInNumber.get(id));
        }
        fu.close();
        // The following method is used to calculate interaction number
        // distribution
        Map<Integer, Integer> intNumberToCount = new HashMap<Integer, Integer>();
        for (Iterator<String> it = proteinToInNumber.keySet().iterator(); it.hasNext();) {
            String protein = it.next();
            Integer intNumber = proteinToInNumber.get(protein);
            Integer count = intNumberToCount.get(intNumber);
            if (count == null)
                intNumberToCount.put(intNumber, 1);
            else
                intNumberToCount.put(intNumber, ++count);
        }
        //fu.setOutput(R3Constants.RESULT_DIR + "IntNumberToCount121707.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "IntNumberToCount122107.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "IntNumberToCount042108_Predicated.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "IntNumberToCount111308.txt");
        //fu.setOutput(R3Constants.RESULT_DIR + "IntNumberToCount072809.txt");
        fu.setOutput(R3Constants.RESULT_DIR + "IntNumberToCount091912.txt");
        fu.printLine("Interaction\tFrequence");
        // Try to sorting a little bit
        int index = ids.size() - 1;
        int min = proteinToInNumber.get(ids.get(index));
        int max = proteinToInNumber.get(ids.get(0));
        for (int i = min; i < max + 1; i++) {
            Integer count = intNumberToCount.get(i);
            if (count == null)
                continue;
            fu.printLine(i + "\t" + count);
        }
        fu.close();
    }    
    
//    @Test
//    public void annotateViaIDNumbers() throws IOException {
//        // These files should be annotated
//        String[] fileNames = new String[] {
//                "TenOfRandomAndTenOfBCR(N).txt",
//                "TwentyOfRandom.txt",
//                "FortyOfRandomAndTenOfBCR.txt",
//                "FiftyOfRandom.txt",
//                "NinetyOfRandomAndTenOfBCR(N).txt",
//                "OneHundredOfRandom.txt"
//        };
//        FileUtility fu = new FileUtility();
//        PathwayAnnotator analyzer = new PathwayAnnotator();
//        String line = null;
//        int index1, index2;
//        for (String fileName : fileNames) {
//            fu.setInput(GRAPH_DIR + fileName);
//            line = fu.readLine();
//            // Test Protein Ids: [P37840, O95558, ..., Q92945, Q96JU7, Q13136]
//            index1 = line.indexOf("[");
//            index2 = line.lastIndexOf("]");
//            line = line.substring(index1 + 1, index2);
//            // Have to add "+". Otherwise, an extract space between "," and " "
//            // will be inserted
//            String[] ids = line.split("(,| )+");
//            Set<String> idSet = new HashSet<String>();
//            for (String id : ids)
//                idSet.add(id);
//            System.out.println("File Name: " + fileName);
//            analyzer.topicAnnotate(idSet, System.out);
//            System.out.println();
//            fu.close();
//        }
//    }
    
    private String generateKey(String id1, String id2) {
        int compare = id1.compareTo(id2);
        if (compare < 0)
            return id1 + "->" + id2;
        else
            return id2 + "->" + id1;
    }
    
//    @Test
//    public void calculateAndAnnotateTSPForSingleFile() throws Exception {
//        TopicAnalyzer analyzer = new TopicAnalyzer();
//        Map<String, Set<String>> fiToTopics = analyzer.loadFIToTopics();
//        String[] fileNames = new String[] {
//                "OneHundredOfRandom.txt", "AnnotateTSPOfOneHundredOfRandom.txt",
//                "TwentyOfRandom.txt", "AnnotateTSPOfTwentyOfRandom.txt",
//                "FiftyOfRandom.txt", "AnnotateTSPOfFiftyOfRandom.txt",
//                "NinetyOfRandomAndTenOfBCR(N).txt", "AnnotateNinetyOfRandomAndTenOfBCR(N).txt",
//                "TenOfRandomAndTenOfBCR(N).txt", "AnnotateTenOfRandomAndTenOfBCR(N).txt",
//                "FortyOfRandomAndTenOfBCR.txt", "AnnotateFortyOfRandomAndTenOfBCR.txt"
//        };
//        for (int k = 0; k < fileNames.length; k++) {
//            //String fileName = GRAPH_DIR + "OneHundredOfRandom.txt";
//            //String fileName = GRAPH_DIR + "TwentyOfRandom.txt";
//            String fileName = GRAPH_DIR + fileNames[k];
//            // Annotate this tsp
//            Map<String, Set<String>> fiToPathMap = loadFIToPath(fileName);
//            //String outFileName = GRAPH_DIR + "AnnotateTSPOfOneHundredOfRandom.txt";
//            //String outFileName = GRAPH_DIR + "AnnotateTSPOfTwentyOfRandom.txt";
//            String outFileName = GRAPH_DIR + fileNames[++k];
//            FileUtility fu = new FileUtility();
//            fu.setOutput(outFileName);
//            // Do five times
//            for (int j = 0; j < 5; j++) {
//                System.out.println("Sample : " + j);
//                List<Gene> tsp = calculateTSPForSingleFile(fileName);
//                List<String> paths = new ArrayList<String>();
//                for (int i = 0; i < tsp.size() - 1; i++) {
//                    String id1 = tsp.get(i).getId();
//                    String id2 = tsp.get(i + 1).getId();
//                    String key = generateKey(id1, id2);
//                    Set<String> step = fiToPathMap.get(key);
//                    if (step == null)
//                        continue;
//                    paths.addAll(step);
//                }
//                // Need to calculate the tail to the head
//                String id1 = tsp.get(tsp.size() - 1).getId();
//                String id2 = tsp.get(0).getId();
//                String key = generateKey(id1, id2);
//                Set<String> step = fiToPathMap.get(key);
//                if (step != null)
//                    paths.addAll(step);
//                Map<String, Integer> topicToNumber = new HashMap<String, Integer>();
//                for (String path : paths) {
//                    String fi = convertEdgeNameToFI(path);
//                    Set<String> topics = fiToTopics.get(fi);
//                    if (topics == null)
//                        continue;
//                    for (String topic : topics) {
//                        Integer c = topicToNumber.get(topic);
//                        if (c == null)
//                            topicToNumber.put(topic, 1);
//                        else
//                            topicToNumber.put(topic, ++c);
//                    }
//                }
//                int pathLength = paths.size();
//                PrintWriter writer = fu.getPrintWriter();
//                writer.print("Path: ");
//                for (Gene g : tsp)
//                    writer.print(g.getId() + ", ");
//                writer.println();
//                outputPathAnnotations(topicToNumber, 
//                                      pathLength, 
//                                      fu.getPrintWriter());
//                writer.println();
//                writer.println();
//            }
//            fu.close();
//        }
//    }
    
    private List<Gene> calculateTSPForSingleFile(String fileName) throws Exception {
        FileUtility fu =  new FileUtility();
        fu.setInput(fileName);
        // Need to load fi -> path 
        Map<String, String[]> fiPairToPath = new HashMap<String, String[]>();
        Map<String, Integer> fiPairToLength = new HashMap<String, Integer>();
        String line = null;
        String pair = null;
        while ((line = fu.readLine()) != null) {
            // Check the title line
            if (line.contains("->") && line.endsWith(":")) {
                pair = line.substring(0, line.length() - 1);
                // Next line is about path
                line = fu.readLine();
                line = line.trim();
                if (line.equals("null")) {
                    fiPairToPath.put(pair, null);
                    // Just a big number
                    fiPairToLength.put(pair, 100000);
                }
                else {
                    // Remove []
                    line = line.substring(1, line.length() - 1);
                    String[] tokens = line.split(",");
                    fiPairToPath.put(pair, tokens);
                    fiPairToLength.put(pair, tokens.length);
                }
            }
        }
        fu.close();
        List<String> ids = extractIds(fiPairToLength);
        FastTSPSolverViaGA tspSolver = new FastTSPSolverViaGA(ids,
                                                              fiPairToLength);
        Chromosome optimal = tspSolver.calculateOptimal();
        System.out.println("Solution: " + optimal.getCost());
        List<Gene> genes = optimal.getGenes();
        return genes;
    }
    
    private List<String> extractIds(Map<String, Integer> distMap) {
        Set<String> set = new HashSet<String>();
        for (Iterator<String> it = distMap.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            int index = key.indexOf("->");
            set.add(key.substring(0, index));
            set.add(key.substring(index + 2));
        }
        return new ArrayList<String>(set);
    }
    
    private Map<String, Set<String>> loadFIToPath(String fileName) throws IOException {
        FileUtility fu =  new FileUtility();
        fu.setInput(fileName);
        Map<String, Set<String>> pathToTopics = new HashMap<String, Set<String>>();
        // Load annoations
        String line = null;
        while ((line = fu.readLine()) != null) {
            // Check the title line
            if (line.contains("->") && line.endsWith(":")) {
                String path = line.substring(0, line.length() - 1);
                // Next line is about path
                line = fu.readLine();
                if (line.contains("null"))
                    continue;
                line = line.trim();
                // Remove []
                line = line.substring(1, line.length() - 1);
                if (line.length() == 0)
                    continue;
                String[] tokens = line.split(",");
                Set<String> topics = new HashSet<String>();
                for (String t : tokens)
                    topics.add(t.trim());
                pathToTopics.put(path, topics);
            }
        }
        fu.close();
        return pathToTopics;
    }
    
//    @Test
//    public void calculateAnnotationsForSingleFile() throws IOException {
//        String fileName = GRAPH_DIR + "NinetyOfRandomAndTenOfBCR(N).txt";
//        //String fileName = GRAPH_DIR + "OneHundredOfRandom.txt";
//        FileUtility fu =  new FileUtility();
//        fu.setInput(fileName);
//        Map<String, Integer> topicToCount = new HashMap<String, Integer>();
//        // Load annoations
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        Map<String, Set<String>> fiToTopics = topicAnalyzer.loadFIToTopics();
//        String line = null;
//        int totalPath = 0;
//        int totalLength = 0;
//        while ((line = fu.readLine()) != null) {
//            // Check the title line
//            if (line.contains("->") && line.endsWith(":")) {
//                totalPath ++;
//                // Next line is about path
//                line = fu.readLine();
//                if (line.contains("null"))
//                    continue;
//                line = line.trim();
//                // Remove []
//                line = line.substring(1, line.length() - 1);
//                String[] tokens = line.split("[,]");
//                totalLength += tokens.length;
//                for (String edge : tokens) {
//                    String fi = convertEdgeNameToFI(edge.trim());
//                    Set<String> topics = fiToTopics.get(fi);
//                    if (topics == null)
//                        continue;
//                    for (String topic : topics) {
//                        Integer c = topicToCount.get(topic);
//                        if (c == null)
//                            topicToCount.put(topic, 1);
//                        else
//                            topicToCount.put(topic, ++c);
//                    }
//                }
//            }
//        }
//        fu.close();
//        PrintWriter writer = new PrintWriter(new OutputStreamWriter(System.out));
//        outputPathAnnotations(topicToCount, totalLength, writer);
//    }
//
//    private void outputPathAnnotations(Map<String, Integer> topicToCount, 
//                                       int totalLength,
//                                       PrintWriter out) throws IOException {
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        Map<String, Double> topicToFiRatio = topicAnalyzer.loadTopicToFIRatio();
//        //DistributionFactory factory = DistributionFactory.newInstance();
//        for (Iterator<String> it = topicToCount.keySet().iterator(); it.hasNext();) {
//            String topic = it.next();
//            double o = topicToFiRatio.get(topic);
//            if (o < 0.001d)
//                continue;
//            Integer c = topicToCount.get(topic);
//            double r = (double) c / totalLength;
//            double enrich = r / o;
//            double pValue = calculatePValue(o,
//                                            totalLength,
//                                            c); //, 
//                                            //factory);
//            out.printf("%s\t%d\t%f6\t%f6\t%f6\t%12.10g%n",
//                              topic,
//                              c,
//                              r,
//                              enrich,
//                              o,
//                              pValue); // Check enrichment
//        }
//        out.println("Total Length: " + totalLength);
//    }
    
    private double calculatePValue(double ratio,
                                   int pathLength,
                                   int both) {
        //return MathUtilities.calOneTailedBinomPValue(ratio, pathLength, both);                           
        return MathUtilities.calculateBinomialPValue(ratio, pathLength, both);
    }
                                   
    
    @Test
    public void extractPathLengthsForSingleFile() throws IOException {
        String fileName = GRAPH_DIR + "FiftyOfRandom.txt";
        String outName = GRAPH_DIR + "FiftyOfRandomPath.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        List<Integer> lengths = new ArrayList<Integer>();
        int totalPath = 0;
        while ((line = fu.readLine()) != null) {
            // Check the title line
            if (line.contains("->") && line.endsWith(":")) {
                totalPath ++;
                // Next line is about path
                line = fu.readLine();
                if (line.contains("null"))
                    continue;
                lengths.add(line.split(",").length);
            }
        }
        fu.close();
        fu.setOutput(outName);
        for (Integer i : lengths)
            fu.printLine(i.toString());
        fu.close();
        System.out.println("Total Path: " + totalPath);
    }
    
    @Test
    public void extractPathAnnotationsFromSingleFile() throws IOException {
        String fileName = GRAPH_DIR + "NinetyOfRandomAndTenOfBCR(N).txt";
        String outName = GRAPH_DIR + "NinetyOfRandomAndTenOfBCR(N)Annot..txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        final Map<String, Integer> topicToNumber = new HashMap<String, Integer>();
        int totalPath = 0;
        while ((line = fu.readLine()) != null) {
            // Check the title line
            if (line.contains("->") && line.endsWith(":")) {
                totalPath ++;
                // Next line is about path
                line = fu.readLine();
                if (line.contains("null"))
                    continue;
                do {
                    line = fu.readLine();
                    if (line.startsWith("Cannot find"))
                        continue;
                    line = line.trim();
                    if (line.startsWith("[")) {
                        // remove []
                        line = line.substring(1, line.length() - 1);
                        if (line.length() == 0)
                            break;
                        String[] tokens = line.split(",");
                        for (String topic : tokens) {
                            topic = topic.trim();
                            Integer c = topicToNumber.get(topic);
                            if (c == null)
                                topicToNumber.put(topic, 1);
                            else
                                topicToNumber.put(topic, ++c);
                        }
                        break;
                    }
                }
                while (true);
            }
        }
        fu.close();
        List<String> topics = new ArrayList<String>(topicToNumber.keySet());
        Collections.sort(topics, new Comparator<String>() {
            public int compare(String t1, String t2) {
                int c1 = topicToNumber.get(t1);
                int c2 = topicToNumber.get(t2);
                return c2 - c1;
            }
        });
        System.out.println("Starting out...");
        for (String topic : topics) {
            Integer c = topicToNumber.get(topic);
            System.out.println(topic + "\t" + c);
        }
        //fu.setOutput(outName);
        //fu.close();
        System.out.println("Total Path: " + totalPath);
    }
    
    @Test
    public void extractPathLengths() throws IOException {
        List<Integer> randomList = new ArrayList<Integer>();
        List<Integer> pathwayList = new ArrayList<Integer>();
        FileUtility fu = new FileUtility();
        File dir = new File(GRAPH_DIR);
        File[] files = dir.listFiles();
        String line = null;
        boolean isRandom = false;
        for (File file : files) {
            fu.setInput(file.getAbsolutePath());
            if (file.getName().startsWith("TenOfRandom"))
                isRandom = true;
            else
                isRandom = false;
            while ((line = fu.readLine()) != null) {
                // Check the title line
                if (line.contains("->") && line.endsWith(":")) {
                    // Next line is about path
                    line = fu.readLine();
                    if (line.contains("null"))
                        continue;
                    if (isRandom)
                        randomList.add(line.split(",").length);
                    else
                        pathwayList.add(line.split(",").length);
                }
            }
            fu.close();
        }
        System.out.println("Total Random: " + randomList.size());
        String outFileName = GRAPH_DIR + "Random.txt";
        fu.setOutput(outFileName);
        for (Integer i : randomList)
            fu.printLine(i.toString());
        fu.close();
        System.out.println("Total Pathway: " + pathwayList.size());
        outFileName = GRAPH_DIR + "Pathway.txt";
        fu.setOutput(outFileName);
        for (Integer i : pathwayList)
            fu.printLine(i.toString());
        fu.close();
    }
    
    protected String convertEdgeNameToFI(String edgeName) {
        // (P01375 : P55957)
        int index = edgeName.indexOf(":");
        // Have to remove the space after ":"
        return edgeName.substring(1, index) + edgeName.substring(index + 2, edgeName.length() - 1);
    }
    
    @Test
    public void generateFilesForCW() throws IOException {
        FileUtility fu = new FileUtility();
        String interactionFile = RESULT_DIR + "FIIntInBigComp.txt";
        // Protein index files
        Set<String> interactions = fu.loadInteractions(interactionFile);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        List<String> idList = new ArrayList<String>(ids);
        // Generate node list
        String nodeFile = RESULT_DIR + "FINodes.txt";
        fu.setOutput(nodeFile);
        int index = 1;
        Map<String, Integer> idToIndex = new HashMap<String, Integer>();
        for (String id : idList) {
            fu.printLine(index + "\t" + id);
            idToIndex.put(id, index);
            index ++;
        }
        fu.close();
        // interaction index files
        String edgeFile = RESULT_DIR  + "FIEdges.txt";
        fu.setOutput(edgeFile);
        int ind1, ind2;
        for (String in : interactions) {
            index = in.indexOf(" ");
            ind1 = idToIndex.get(in.substring(0, index));
            ind2 = idToIndex.get(in.substring(index + 1));
            fu.printLine(ind1 + "\t" + ind2 + "\t" + "1"); // Always 1
        }
        fu.close();
    }
    
    /**
     * This method is used to find a clique around ZNF genes.
     * @param fiFileName
     * @return
     * @throws IOException
     */
    public Collection<String> searchZNFClique(String fiFileName) throws IOException {
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(fiFileName);
        Set<String> totalFIs = new HashSet<String>(interactions);
        // Want to do a filtering based on ZNF genes
        for (Iterator<String> it = interactions.iterator(); it.hasNext();) {
            String fi = it.next();
            int index = fi.indexOf("\t");
            String name1 = fi.substring(0, index);
            String name2 = fi.substring(index + 1);
            if (!name1.startsWith("ZNF") || !name2.startsWith("ZNF")) {
                it.remove();
            }
        }
        SimpleCliqueSearcher cliqueSearcher = new SimpleCliqueSearcher();
        Collection<String> foundClique = cliqueSearcher.searchClique(interactions);
        cliqueSearcher.validateClique(foundClique, totalFIs);
        return foundClique;
    }
    
    /**
     * This method is used to do cliques searching.
     * @throws Exception
     */
    @Test
    public void analyzeCliques() throws Exception {
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalFIs = new HashSet<String>(interactions);
        System.out.println("Total interactions: " + interactions.size());
        // Want to do a filtering based on ZNF genes
        for (Iterator<String> it = interactions.iterator(); it.hasNext();) {
            String fi = it.next();
            int index = fi.indexOf("\t");
            String name1 = fi.substring(0, index);
            String name2 = fi.substring(index + 1);
            if (!name1.startsWith("ZNF") || !name2.startsWith("ZNF")) {
                it.remove();
            }
        }
        System.out.println("Total interactions among ZNF genes: " + interactions.size());
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        System.out.println("Total genes in checking FIs: " + ids.size());
        SimpleGraph<String, DefaultEdge> graph = (SimpleGraph<String, DefaultEdge>) JGraphTUtilities.createGraph(ids, interactions);
        
        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<String, DefaultEdge>(graph);
        List<Set<String>> components = inspector.connectedSets();
        System.out.println("Components: " + components.size());
        for (int i = 0; i < components.size(); i++) {
            Set<String> component = components.get(i);
            System.out.println(i + ": " + component.size() + " " + component);
        }
        
        // Check edges
        Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(interactions);
        for (String id : idToPartners.keySet()) {
            Set<String> partners = idToPartners.get(id);
            System.out.println(id + "\t" + partners.size());
            Integer count = counts.get(partners.size());
            if (count == null)
                counts.put(partners.size(), 1);
            else
                counts.put(partners.size(), ++count);
        }
        for (Integer c : counts.keySet())
            System.out.println(c + "\t" + counts.get(c));
        // The following code cannot work!
//        // This is a very simple way to find a clique for ZNF: based on the ouptut from the above,
//        // 359 genes should be linked around ZNF317. Check all genes linked to ZNF317 to see if they are mutally linked
//        Set<String> clique = idToPartners.get("ZNF317");
//        clique.add("ZNF317");
//        for (String id1 : clique) {
//            for (String id2 : clique) {
//                int compare = id1.compareTo(id2);
//                String fi = null;
//                if (compare > 0)
//                    fi = id2 + "\t" + id1;
//                else if (compare < 0)
//                    fi = id1 + "\t" + id2;
//                if (fi == null)
//                    continue;
//                if (!interactions.contains(fi)) {
//                    System.out.println("No interactions: " + fi);
//                }
//            }
//        }
        SimpleCliqueSearcher cliqueSearcher = new SimpleCliqueSearcher();
        Collection<String> foundClique = cliqueSearcher.searchClique(interactions);
        System.out.println("Size of clique: " + foundClique.size());
        cliqueSearcher.validateClique(foundClique, totalFIs);
        System.out.println("Size after validating: " + foundClique.size());
        System.out.println(foundClique);
        if (true)
            return;
        // Re-generate graph from the biggest component
        Set<String> comp = components.get(0);
        Set<String> fisInComp = InteractionUtilities.getFIs(comp, interactions);
        graph = (SimpleGraph<String, DefaultEdge>) JGraphTUtilities.createGraph(comp, fisInComp);
        long time1 = System.currentTimeMillis();
        BronKerboschCliqueFinder<String, DefaultEdge> cliqueFinder = new BronKerboschCliqueFinder<String, DefaultEdge>(graph);
        Collection<Set<String>> cliques = cliqueFinder.getAllMaximalCliques();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for finding: " + (time2 - time1));
        int index = 1;
        for (Set<String> clique1 : cliques) {
            System.out.println("Clique " + index + ": " + clique1.size());
            index ++;
        }
    }
    
    /**
     * Analyze graph components and output the biggest graph component.
     * @throws Exception
     */
    @Test
    public void analyzeComponents() throws Exception {
        //String interactionFile = RESULT_DIR + "FIInteractions67.txt";
        //String interactionFile = RESULT_DIR + "CHAGO-K-1_Int_Or.txt";
        //String interactionFile = R3Constants.INTERACTION_FILE_NAME;
        //String interactionFile = R3Constants.RESULT_DIR + "FI73InGene_102908.txt";
        //String interactionFile = R3Constants.RESULT_DIR + "FI73InGene_Remove_Hubs_110308.txt";
        //String interactionFile = R3Constants.RESULT_DIR + "FI73InGene_062008_NO_P_I.txt";
        //String interactionFile = R3Constants.RESULT_DIR + "FI73InGene_061008.txt";
        // All names are in upper case
        //String interactionFile = R3Constants.RESULT_DIR + "FI73InGeneUpperCase_111208.txt";
//        String interactionFile = R3Constants.RESULT_DIR + "FIsInGene_041709.txt";
//        String interactionFile = R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt";
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
//        FIFileAnalyzer fiAnalyzer = new FIFileAnalyzer();
//        Set<String> interactions = fiAnalyzer.loadFIs();
        System.out.println("Total interactions: " + interactions.size());
        List<Set<String>> componentList = calculateGraphComponents(interactions);
        System.out.println("Total components: " + componentList.size());
        int index = 0;
        for (Set<String> comp : componentList) {
            System.out.println(index + ": " + comp.size());
            index ++;
        }
        Set<String> biggestComp = componentList.get(0);
        // Want to print out the biggest components
        Set<String> interactionsInBiggest = new HashSet<String>();
        index = 0;
        for (String in : interactions) {
            index = in.indexOf("\t");
            if (index < 0)
                index = in.indexOf(" ");
            String id1 = in.substring(0, index);
            String id2 = in.substring(index + 1);
            if (biggestComp.contains(id1) &&
                biggestComp.contains(id2))
                interactionsInBiggest.add(in);
        }
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "CHAGO-K-1_Int_Or_Big_Comp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FIInteractions60_121707_BigComp.txt");   
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FI73_Predicated_042108_BigComp.txt"); 
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FI73InGene_102908_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FI73InGene_Remove_Hubs_110308_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FI73InGene_062008_NO_P_I_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FI73InGene_061008_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FI73InGeneUpperCase_111208_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FIsInGene_041709_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FIsInGene_Pathway_060109_BigComp.txt");
        //fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FIs_042109_BigComp.txt");
        fu.saveInteractions(interactionsInBiggest, RESULT_DIR + "FIsInGene_No_ZNF_042810_BigComp.txt");
    }

    /**
     * Calculate linked graph components from a set of FIs.
     * @param interactions
     * @return
     */
    public List<Set<String>> calculateGraphComponents(Set<String> interactions) {
        return JGraphTUtilities.calculateGraphComponents(interactions);
    }
    
    @Test
    public void checkSimpleStatistics() throws IOException {
        String intFileName =  RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = new FileUtility().loadInteractions(intFileName);
        System.out.println("Total interactions: " + fis.size());
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total ids: " + ids.size());
    }
    
    private Set<String> getHubs(int edgeNum) throws Exception{
    	FileUtility fu = new FileUtility();
    	Map<String,Set<String>> pairsMap = fu.loadSetMap("data/Reactome/ReactomeInteractions_GENEID.txt", 
    														"\t", true);
    	Set<String> keySet = pairsMap.keySet();
    	Set<String> hubSet = new HashSet<String>();
    	System.out.println("\nThe following hubs are found"+"(The edge number is "+edgeNum+"):");
    	for(String key : keySet){
    		Set<String> edgeSet = pairsMap.get(key);
    		if(edgeSet != null)
    			if(edgeSet.size() > edgeNum) {
    				//System.out.println(key);
    				hubSet.add(key);
    			}
    	}
    	fu.saveInteractions(hubSet, "data/GraphAnalysis/ReactomeHubs"+edgeNum+".txt");
    	return hubSet;
    }
    @Test
    public void getHubStatistics() throws Exception {
    	FileUtility fu = new FileUtility();
    	Map<String,String> omimMap = fu.loadInteractionPairs( "data/OMIM/mimGeneIDMap.txt", "\t",true);
    	int[] edgeNums = new int[] {
    		10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,
    		300,310,320,330,340,350,360,370,380,390,400
    	};
    	List<String> result = new ArrayList<String>();
    	FileUtility fu2 = new FileUtility();
    	fu2.setOutput( "data/OMIM/edgeNumVSChance.txt");
    	for(int edgeNum : edgeNums) {
    		Set<String> hubSet = getHubs(edgeNum);
    		int debugCounter = 0;
    		for(String hub : hubSet){
    			String map = omimMap.get(hub);
    			//if(map != null)
    				//System.out.println(hub);
    			if(map == null)
    				debugCounter++;
    		}
    		System.out.println("The edge number is: "+edgeNum);
    		System.out.println("The hub set size is: "+hubSet.size());
    		System.out.println("The disease percentage: "+(hubSet.size()-debugCounter)*1.0/hubSet.size());
    		fu2.printLine(edgeNum+"\t"+(hubSet.size()-debugCounter)*1.0/hubSet.size()+"\t"+hubSet.size());
    	}
    	fu2.close();
    }
    
    public Graph<String, DefaultEdge> createGraph(String intFileName) throws IOException {
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        SimpleGraph<String, DefaultEdge> graph = (SimpleGraph<String, DefaultEdge>) JGraphTUtilities.createGraph(ids, interactions);
        return graph;
    }
    
    protected Graph<String, DefaultEdge> createRandomGraph(Set<String> ids,
                                                           int totalEdges) {
        Graph<String, DefaultEdge> graph = new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);
        for (String id : ids) 
            graph.addVertex(id);
        // Need to generate enough edge numbers
        // want to control edge by this method 
        Set<String> edgeSet = new HashSet<String>();
        List<String> idList = new ArrayList<String>(ids);
        String key = null;
        while (edgeSet.size() < totalEdges) {
            int index1 = (int) (Math.random() * idList.size());
            int index2 = (int) (Math.random() * idList.size());
            if (index1 == index2)
                continue;
            if (index1 < index2)
                key = index1 + " " + index2;
            else
                key = index2 + " " + index1;
            if (edgeSet.contains(key))
                continue;
            edgeSet.add(key);
            graph.addEdge(idList.get(index1),
                          idList.get(index2));
        }
        return graph;
    }
}
