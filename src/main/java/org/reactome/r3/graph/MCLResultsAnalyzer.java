/*
 * Created on Oct 23, 2006
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.distribution.HypergeometricDistribution;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import junit.framework.TestCase;

public class MCLResultsAnalyzer extends TestCase {
    private final double P_VALUE_CUTOFF = 0.0001d;
    //private final String RESULT_DIR = "results/interaction/mclResults/";
    //private final String RESULT_DIR = "results/v2/mcl/";
    //private final String RESULT_DIR = "datasets/BreastCancer/";
    //private final String CLUSTERING_RESULT_FILE = "FIsWithWeightFromTTestMCLClustersGSE2034_I80_091509.txt";
    private final String RESULT_DIR = "results/v3/";
    private final String CLUSTERING_RESULT_FILE = "MCLCluster_FIsInGene_041709_I40.txt";
    //private final String CLUSTERING_RESULT_FILE = "FI73_042108.class_I24_c26";
    //private final String CLUSTERING_RESULT_FILE = "FIInteractions67.class_I24_c26";
    //private final String CLUSTERING_RESULT_FILE = "FIInteractionsWithComplexAndSet.cluster.I24s6c26";
    private int TOTAL_PROTEIN = 15199;
    private FileUtility fu;
    
    public MCLResultsAnalyzer() {
        fu = new FileUtility();
        try {
            Set<String> interactions = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
            TOTAL_PROTEIN = InteractionUtilities.grepIDsFromInteractions(interactions).size();
        }                                                
        catch(IOException e) {
            e.printStackTrace();
        }
    }
    
//    public void checkInteractionsFromRandomFile() throws Exception {
//        String dir = "/Users/wgm/Documents/gkteam/marcela/";
//        String fileName = dir + "random95c.txt";
//        Set<String> ids = fu.loadSet(fileName);
//        StringBuilder builder = new StringBuilder();
//        for (Iterator<String> it = ids.iterator(); it.hasNext();) {
//            builder.append(it.next());
//            if (it.hasNext())
//                builder.append(",");
//        }
//        Service serviceModel = new ObjectServiceFactory().create(InteractionService.class);
//        InteractionService service = (InteractionService) new XFireProxyFactory().create(serviceModel, 
//                                                                                         "http://localhost:8080/caBigR3WebApp/FIService");
//        
//        List<Interaction> interactions = service.queryForAnd(builder.toString());
//        System.out.println("Interactions: " + interactions.size());
//        // List interactions occur in the same clusters
//        List<ClusterInfo> clusters = loadClusters();
//        // Get rid of clusters having only one proteins
//        for (Iterator<ClusterInfo> it = clusters.iterator(); it.hasNext();) {
//            ClusterInfo c = it.next();
//            if (c.proteinNumber < 2)
//                it.remove();
//        }
//        List<Interaction> inClusterList = new ArrayList<Interaction>();
//        for (Iterator<Interaction> it = interactions.iterator(); it.hasNext();) {
//            Interaction i = it.next();
//            String id1 = i.getFirstProtein().getPrimaryAccession();
//            String id2 = i.getSecondProtein().getPrimaryAccession();
//            for (ClusterInfo c : clusters) {
//                if (c.proteins.contains(id1) && c.proteins.contains(id2)) {
//                    inClusterList.add(i);
//                    it.remove();
//                    break;
//                }
//            }
//        }
//        System.out.println("In Cluster Interactions: " + inClusterList.size());
//        // Output
//        fileName = fileName + ".result";
//        fu.setOutput(fileName);
//        fu.printLine("Interactions occuring in the same clusters:");
//        for (Interaction i : inClusterList) {
//            fu.printLine(i.getFirstProtein().getPrimaryAccession() + " " +
//                         i.getSecondProtein().getPrimaryAccession());
//        }
//        fu.printLine("");
//        fu.printLine("Interactions occuring not in the same clusters:");
//        //TODO: This cannot work. It seems all Interactions are the same!!!
//        //interactions.removeAll(inClusterList);
//        for (Interaction i : interactions) {
//            fu.printLine(i.getFirstProtein().getPrimaryAccession() + " " + 
//                         i.getSecondProtein().getPrimaryAccession());
//        }
//        fu.close();
//    }
//    
    /**
     * This method is used to generate p values for pathways in clusters.
     * @throws IOException
     */
    public void generatePValuesForPathwaysInClusters() throws IOException {
        // Load the map from pathways to protein numbers
        Map<String, Integer> pathwayToProtein = new HashMap<String, Integer>();
        String fileName = RESULT_DIR + "ProteinNumberInTopics052708.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // 1    426 Wnt signaling pathway(P)
            pathwayToProtein.put(tokens[2], Integer.parseInt(tokens[1]));
        }
        fu.close();
        // Load the map from cluster to protein numbers
        Map<Integer, Integer> clusterToProtein = new HashMap<Integer, Integer>();
        //fileName = RESULT_DIR + "ProteinsInClustersWithComplexAndSet.txt";
        fileName = RESULT_DIR + "ProteinsInClusters041208.txt";
        fu.setInput(fileName);
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // 1    298
            clusterToProtein.put(Integer.parseInt(tokens[0]),
                                 Integer.parseInt(tokens[1]));
        }
        fu.close();
        // Now calculation started
        int total = TOTAL_PROTEIN;
        //fileName = RESULT_DIR + "ClusterToTopicWithComplexAndSet.txt";
        fileName = RESULT_DIR + "ClusterToTopic041208.txt";
        fu.setInput(fileName);
        line = fu.readLine();
        // This is the file for outputting
        //String output = RESULT_DIR + "ClusterToTopicWithPValueWithComplexAndSet.txt";
        String output = RESULT_DIR + "ClusterToTopicWithPValue041208.txt";
        FileUtility outputFu = new FileUtility();
        outputFu.setOutput(output);
        outputFu.printLine("ClusterIndex\tProteinsInPathway\tP Value\tPathwayInCluster");
        int proteinInCluster;
        int proteinInPathway;
        int proteinInBoth;
        int cluster;
        double p;
        int clusterSizeCutoff = 9;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // 1    193 HIV Infection(R)
            cluster = Integer.parseInt(tokens[0]);
            //if (cluster > 60)
            //    break; // Only the first 60 clusters (protein numbers > 10)
            proteinInBoth = Integer.parseInt(tokens[1]);
            proteinInCluster = clusterToProtein.get(cluster);
            if (proteinInCluster < clusterSizeCutoff)
                break;
            proteinInPathway = pathwayToProtein.get(tokens[2]);
            p = calculatePValue(total, 
                                proteinInBoth, 
                                proteinInCluster,
                                proteinInPathway);
            outputFu.printLine(tokens[0] + "\t" + 
                               tokens[1] + "\t" + 
                               p + "\t" + 
                               tokens[2]);
        }
        outputFu.close();
    }
    
    private double calculatePValue(int total,
                                   int proteinInBoth,
                                   int proteinInCluster,
                                   int proteinInPathway) {
        double p = 0.0d;
        HypergeometricDistribution dist = new HypergeometricDistributionImpl(total, 
                                                                             proteinInPathway, 
                                                                             proteinInCluster);
//        HypergeometricDistribution dist = factory.createHypergeometricDistribution(total, 
//                                                                                   proteinInPathway, 
//                                                                                   proteinInCluster);
        double p1;
        for (int i = proteinInBoth; i < proteinInCluster + 1; i++) {
            p1 = dist.probability(i);
            //System.out.println("pdf: " + p1);
            p += p1;
        }
        return p;
    }
    
    public void generateCuratedTopicToClusterFile() throws IOException {
        FileUtility fu = new FileUtility();
        //String fileName = RESULT_DIR + "ClusterToTopicWithPValueSorted.txt";
        String fileName = RESULT_DIR + "ClusterToTopicWithPValue041208Sorted.txt";
        fu.setInput(fileName);
        // Get output for p value < 10-4 only.
        FileUtility outputFu = new FileUtility();
        String outFileName = RESULT_DIR + "ClusterToTopicCurated041208.txt";
        outputFu.setOutput(outFileName);
        String line = fu.readLine();
        outputFu.printLine(line);
        double pvalue;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            //1 193 2.773032056652674E-245  HIV Infection(R)
            pvalue = Double.parseDouble(tokens[2]);
            if (pvalue > P_VALUE_CUTOFF)
                continue;
            outputFu.printLine(line);
        }
        fu.close();
        outputFu.close();
    }
    
    /**
     * This method is used to generate a comprehensive cluster information file by
     * combining the most enrichment pathways, protein numbers and cliqueness.
     * @throws IOException
     */
    public void generateComprehensiveClusterInfoFile() throws IOException {
        List<ClusterInfo> allClusters = loadClusters();
        String pValueFileName = RESULT_DIR + "ClusterToTopicCurated041208.txt";
        fu.setInput(pValueFileName);
        String line = fu.readLine();
        int currentIndex = 0;
        List<TopicInClusterInfo> clusters = new ArrayList<TopicInClusterInfo>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            int index = Integer.parseInt(tokens[0]);
            if (index != currentIndex) {
                TopicInClusterInfo cluster = new TopicInClusterInfo();
                cluster.clusterIndex = index;
                cluster.pvalue = Double.parseDouble(tokens[2]);
                cluster.numbder = new Integer(tokens[1]);
                cluster.id = tokens[3];
                currentIndex = index;
                clusters.add(cluster);
            }
        }
        fu.close();
        String cliquenessFileName = RESULT_DIR + "ClusterToCliqueness041208.txt";
        fu.setInput(cliquenessFileName);
        line = fu.readLine();
        String outFileName = RESULT_DIR + "ClusterWithCompInfo041208.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        outFu.printLine("Cluster\tProteinInCluster\tFIsInCluster\tCliqueness\tProteinInPathway\tP_value\tPathway");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            int index = Integer.parseInt(tokens[0]);
            if (index > clusters.size() - 1)
                break;
            int number = Integer.parseInt(tokens[1]);
            double cliqueness = Double.parseDouble(tokens[2]);
            // Clusters are indexed from 1
            TopicInClusterInfo cluster = clusters.get(index - 1);
            ClusterInfo clusterInfo = allClusters.get(index - 1);
            outFu.printLine(index + "\t" + number + "\t" + clusterInfo.interactionNumber + "\t" + 
                            cliqueness + "\t" +
                            cluster.numbder + "\t" + cluster.pvalue + "\t" + cluster.id);
        }
        outFu.close();
        fu.close();
    }
    
    void assignPathwayAnnotations(List<ClusterInfo> clusters) throws IOException {
        FileUtility fu = new FileUtility();
        String fileName = RESULT_DIR + "ClusterToTopicWithPValueSorted.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        int index = 0;
        int preIndex = 0;
        //int proteinNumber;
        double pValue;
        String pathwayName;
        ClusterInfo info = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            //1   193 2.773032056652674E-245  HIV Infection(R)
            index = Integer.parseInt(tokens[0]);
            if (index != preIndex) {
                info = clusters.get(index - 1);
                preIndex = index;
            }
            //proteinNumber = Integer.parseInt(tokens[1]);
            // Choose protein numbers should be no less than 10% of total proteins in clusters
            //if (proteinNumber < info.proteinNumber * 0.1)
            //    continue;
            // P-value should not be higher than 0.0001
            pValue = Double.parseDouble(tokens[2]);
            if (pValue > P_VALUE_CUTOFF)
                continue;
            info.addTopic(tokens[3]);
        }
        fu.close();
    }
    
    public Map<String, Integer> getClusterForProteins(Set<String> ids) throws IOException {
        List<ClusterInfo> clusters = loadClusters();
        Map<String, Integer> idToCluster = new HashMap<String, Integer>();
        for (String id : ids) {
            for (ClusterInfo c : clusters) {
                if (c.proteins.contains(id))
                    idToCluster.put(id, c.index);
            }
        }
        return idToCluster;
    }
    
    public void getClusterForProtein() throws IOException {
        List<ClusterInfo> clusters = loadClusters();
        assignPathwayAnnotations(clusters);
        String[] ids = new String[] {
                "O43290",
                "P10909",
                "P21953",
                "P30042",
                "P40925",
                "P41222",
                "P78346",
                "Q14146",
                "Q8WXD2",
                "Q96E40",
                "Q96IT1",
                "Q99731",
                "Q9HAU5",
                "Q9Y5K5",
                "Q9Y6I8"
        };
        for (String id : ids) {
            for (ClusterInfo c : clusters) {
                if (c.proteins.contains(id)) {
                    System.out.println(id + " in " + c.index + " (" + c.proteinNumber + ")");
                    //print out annotation
                    if (c.topics != null) {
                        for (String t : c.topics) 
                            System.out.println("    " + t);
                    }
                    else
                        System.out.println("    no pathways assigned.");
                    break;
                }
            }
        }
    }
    
    public void mergePathwayInfoToClusterInfo() throws IOException {
        // Load Cluster Information
        FileUtility fu = new FileUtility();
        String fileName = RESULT_DIR + "I24C26ClusteringInfoSorted.txt";
        fu.setInput(fileName);
        Map<Integer, ClusterInfo> infoMap = new HashMap<Integer, ClusterInfo>();
        String line = fu.readLine(); // Escape the header line
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            ClusterInfo info = new ClusterInfo();
            info.index = Integer.parseInt(tokens[0]);
            info.proteinNumber = Integer.parseInt(tokens[1]);
            info.interactionNumber = Integer.parseInt(tokens[2]);
            info.interactionPerProtein = Double.parseDouble(tokens[3]);
            info.sortedIndex = Integer.parseInt(tokens[4]);
            infoMap.put(info.index, info);
        }
        fu.close();
        // Load cluster pathway information
        fileName = RESULT_DIR  + "ClusterToTopic.txt";
        fu.setInput(fileName);
        line = fu.readLine(); // Escape the first header line
        int index = 0;
        ClusterInfo info = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            index = Integer.parseInt(tokens[0]);
            if (info == null || info.index != index) 
                info = infoMap.get(index);
            // Want to merge number and pathway name together
            info.addTopic(tokens[2] + "(" + tokens[1] + ")");
        }
        fu.close();
        // This should be the total number
        int total = infoMap.size();
        fileName = RESULT_DIR + "I24C26ClusteringInfoWithPathway.txt";
        fu.setOutput(fileName);
        fu.printLine("Cluster\tProtein\tPathway\tPathway1\tPathway2\tPathway3\tInteraction\tInteraction/Protein\tSortedIndex");
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < total; i++) {
            info = infoMap.get(i + 1);
            if (info.proteinNumber < 2)
                continue; // Escape clusters containing one protein only
            builder.setLength(0);
            builder.append(info.index).append("\t");
            builder.append(info.proteinNumber).append("\t");
            if (info.topics == null)
                builder.append("\t\t\t\t");
            else {
                builder.append(info.topics.size()).append("\t");
                for (int j = 0; j < info.topics.size() && j < 3; j++)
                    builder.append(info.topics.get(j)).append("\t");
                // Add enough "\t"
                for (int j = info.topics.size(); j < 3; j++) {
                    builder.append("\t");
                }
            }
            builder.append(info.interactionNumber).append("\t");
            builder.append(info.interactionPerProtein).append("\t");
            builder.append(info.sortedIndex);
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    protected Map<String, Integer> assignClusterToProtein(List<ClusterInfo> clusters) {
        Map<String, Integer> proteinToCluster = new HashMap<String, Integer>();
        for (ClusterInfo info : clusters) {
            for (String id : info.proteins) {
                proteinToCluster.put(id, info.index);
            }
        }
        return proteinToCluster;
    }
    
    /**
     * Filter any annotations whose p-values > 0.0001 and the number of proteins
     * for that annotation is less than 10%.
     * @throws IOException
     */
    public void filterTopicsInClusters() throws IOException {
        List<ClusterInfo> clusters = loadClusters();
        // Remove all topic annotaions whose p-value > 0.0001
        FileUtility fu = new FileUtility();
        String inputFile = RESULT_DIR + "ClusterToTopicWithPValueWithComplexAndSetSorted.txt";
        fu.setInput(inputFile);
        FileUtility outFu = new FileUtility();
        String outputFile = RESULT_DIR + "ClusterToTopicFiltered.txt";
        outFu.setOutput(outputFile);
        String line = fu.readLine();
        outFu.printLine(line);
        int cluster = 0;
        int proteinCutOff = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            int tmp = Integer.parseInt(tokens[0]);
            if (tmp != cluster) {
                cluster = tmp;
                ClusterInfo info = clusters.get(cluster);
                proteinCutOff = (int) (info.proteinNumber * 0.1);
            }
            Double pValue = Double.valueOf(tokens[2]);
            int proteinInAnnot = Integer.parseInt(tokens[1]);
            if (pValue < 0.0001d && proteinInAnnot > proteinCutOff)
                outFu.printLine(line);
        }
        fu.close();
        outFu.close();
    }
    
    /**
     * This method is used to assign pathways to clusters.
     * @throws IOException
     */
    public void assignTopicsToClusters() throws Exception {
        // Load map first
        FileUtility fu = new FileUtility();
        //Map<String, Set<String>> idToTopic = fu.loadSetMap(R3Constants.PROTEIN_ID_TO_TOPIC);
        Map<String, Set<String>> idToTopic = fu.loadSetMap(R3Constants.GENE_TO_TOPIC);
        Map<String, Set<String>> topicToIds = InteractionUtilities.switchKeyValues(idToTopic);
        // Load the actual clusters
        String fileName = RESULT_DIR + CLUSTERING_RESULT_FILE;
        List<Set<String>> clusters = loadMCLClusters(fileName);
        // Output
        //fileName = RESULT_DIR + "ClusterToTopic.txt";
        //fileName = RESULT_DIR + "ClusterToTopicWithComplexAndSet.txt";
        //fileName = RESULT_DIR + "ClusterToTopic041208.txt";
        int total = TOTAL_PROTEIN;
        fileName = RESULT_DIR + "MCLClusters_FIsInGene_041709_I40ToTopic091609.txt";
        fu.setOutput(fileName);
        // Generate header
        fu.printLine("ClusterIndex\tProteins\tProteinsInPathway\tp-value\tPathwayInCluster");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            if (cluster.size() < 3)
                continue;
            List<String> lines = new ArrayList<String>();
            for (String t : topicToIds.keySet()) {
                // Want to calculate the p-value
                Set<String> idsInTopic = topicToIds.get(t);
                Set<String> shared = new HashSet<String>(idsInTopic);
                shared.retainAll(cluster);
                double pvalue = MathUtilities.calculateHypergeometricPValue(total, 
                                                                            cluster.size(), 
                                                                            idsInTopic.size(), 
                                                                            shared.size());
                if (pvalue > 1.0E-2)
                    continue;
                lines.add(i + "\t" + 
                          cluster.size() + "\t" +
                          shared.size() + "\t" + 
                          pvalue + "\t" +
                          t);
            }
            Collections.sort(lines, new Comparator<String>() {
                public int compare(String line1, String line2) {
                    String[] tokens = line1.split("\t");
                    Double pvalue1 = new Double(tokens[3]);
                    tokens = line2.split("\t");
                    Double pvalue2 = new Double(tokens[3]);
                    return pvalue1.compareTo(pvalue2);
                }
            });
            for (String line : lines) {
                fu.printLine(line);
            }
        }
        fu.close();
    }
    
    private Map<String, Integer> countTopics(List<String> topics) {
        Map<String, Integer> topicCount = new HashMap<String, Integer>();
        for (String t : topics) {
            Integer i = topicCount.get(t);
            if (i == null) {
                topicCount.put(t, 1);
            }
            else {
                topicCount.put(t, ++i);
            }
        }
        return topicCount;
    }
    
//    public void generateProteinToGOMap() throws IOException {
//        // There should be two mappings
////        Map<String, Set<String>> protein2GOF = new HashMap<String, Set<String>>();
////        Map<String, Set<String>> protein2GOP = new HashMap<String, Set<String>>();
//        GODataAnalyzer goAnalyzer = new GODataAnalyzer();
////        goAnalyzer.generateProtein2GOMap(protein2GOF, protein2GOP);
//        FileUtility fu = new FileUtility();
////        fu.saveSetMap(protein2GOF, RESULT_DIR + "ProteinIdToGOF.txt");
////        fu.saveSetMap(protein2GOP, RESULT_DIR + "ProteinIdToGOP.txt");
//        Map<String, Set<String>> protein2GOC = new HashMap<String, Set<String>>();
//        goAnalyzer.generateProtein2GOCCMap(protein2GOC);
//        fu.saveSetMap(protein2GOC, RESULT_DIR + "ProteinIdToGOC.txt");
//    }
//    
//    public void sortAnnotationsInClusters() throws IOException {
//        String[] fileNames = new String[] {
//                //"ClusterToGOBP.txt",        
//                //"ClusterToGOCC.txt",               
//                //"ClusterToTopicWithPValueWithComplexAndSet.txt",
//                //"ClusterToGOMF.txt"
//                "ClusterToTopicWithPValue041208.txt"
//        };
//        Map<String, String> goIdToTerm = new GODataAnalyzer().loadGOIdToTermMap();
//        FileUtility fu = new FileUtility();
//        FileUtility outFu = new FileUtility(); // This is for output
//        for (String fileName : fileNames) {
//            fileName = RESULT_DIR + fileName;
//            fu.setInput(fileName);
//            String line = fu.readLine();
//            int currentIndex = 1;
//            List<TopicInClusterInfo> topics = new ArrayList<TopicInClusterInfo>();
//            String outFileName = fileName.substring(0, fileName.length() - 4) + "Sorted.txt";
//            outFu.setOutput(outFileName);
//            boolean isForGo = fileName.contains("GO");
//            if (isForGo)
//                outFu.printLine("ClusterIndex\tProteinsInGO\tP-value\tGOInCluster\tGOTerm");
//            else // For pathway
//                outFu.printLine("ClusterIndex\tProteinsInPathway\tP-value\tPathwayInCluster");
//            int cluster;
//            Comparator<TopicInClusterInfo> sorter = new Comparator<TopicInClusterInfo>() {
//                public int compare(TopicInClusterInfo t1,
//                                   TopicInClusterInfo t2) {
//                    Double p1 = t1.pvalue;
//                    Double p2 = t2.pvalue;
//                    return p1.compareTo(p2);
//                }
//            };
//            while ((line = fu.readLine()) != null) {
//                String[] tokens = line.split("\t");
//                // 1    52  2.8874891006951427E-35  GO:0003723
//                cluster = Integer.parseInt(tokens[0]);
//                if (cluster != currentIndex) {
//                    // Need to output
//                    Collections.sort(topics, sorter);
//                    for (TopicInClusterInfo t : topics) {
//                        if (isForGo) {
//                            outFu.printLine(t.clusterIndex + "\t" +
//                                            t.numbder + "\t" +
//                                            t.pvalue + "\t" + 
//                                            t.id + "\t" + 
//                                            goIdToTerm.get(t.id));
//                        }
//                        else 
//                            outFu.printLine(t.clusterIndex + "\t" +
//                                            t.numbder + "\t" +
//                                            t.pvalue + "\t" + 
//                                            t.id);
//                    }
//                    topics.clear();
//                }
//                TopicInClusterInfo topic = new TopicInClusterInfo();
//                topic.clusterIndex = cluster;
//                currentIndex = cluster;
//                topic.numbder = Integer.parseInt(tokens[1]);
//                topic.pvalue = Double.parseDouble(tokens[2]);
//                topic.id = tokens[3];
//                topics.add(topic);
//            }
//            fu.close();
//            outFu.close();
//        }
//    }
//    
    
    /**
     * This method is used to generate protein ids to go ids mappping.
     * @throws IOException
     */
    public void assignGOTermsToClustersInfo() throws IOException {
        FileUtility fu = new FileUtility();
        String goDir = "/Users/wgm/Documents/EclipseWorkspace/caBigR3/results/interaction/mclResults/";
        // There should be two mappings
        // The following is for GO MF
//        Map<String, Set<String>> protein2GOF = fu.loadSetMap(goDir + "ProteinIdToGOF.txt");
//        // Need to retain proteins occur in clusters only so that the total
//        // proteins in clusters are correct
//        Set<String> totalIds = new HashSet<String>();
//        List<ClusterInfo> clusters = loadClusters();
//        for (ClusterInfo c : clusters) {
//            totalIds.addAll(c.proteins);
//        }
//        for (Iterator<String> it = protein2GOF.keySet().iterator(); it.hasNext();) {
//            String id = it.next();
//            if (!totalIds.contains(id))
//                it.remove();
//        }
//        Map<String, Integer> proteinsInGOF = countGOUsage(protein2GOF);
////        exportProteinNumberMap(fu, 
////                               RESULT_DIR + "ProteinNumberInGOCC.txt",
////                               proteinsInGOF);
//        assignGOTermsToClusters(clusters, 
//                                protein2GOF,
//                                proteinsInGOF, 
//                                RESULT_DIR + "ClusterToGOMF.txt");
//        // The following is for GO BP
//        Map<String, Set<String>> protein2GOP = fu.loadSetMap(goDir + "ProteinIdToGOP.txt");
//        // Need to retain proteins occur in clusters only so that the total
//        // proteins in clusters are correct
//        Set<String> totalIds = new HashSet<String>();
//        List<ClusterInfo> clusters = loadClusters();
//        for (ClusterInfo c : clusters) {
//            totalIds.addAll(c.proteins);
//        }
//        for (Iterator<String> it = protein2GOP.keySet().iterator(); it.hasNext();) {
//            String id = it.next();
//            if (!totalIds.contains(id))
//                it.remove();
//        }
//        Map<String, Integer> proteinsInGOP = countGOUsage(protein2GOP);
////        exportProteinNumberMap(fu, 
////                               RESULT_DIR + "ProteinNumberInGOBP.txt",
////                               proteinsInGOF);
//        assignGOTermsToClusters(clusters, 
//                                protein2GOP,
//                                proteinsInGOP, 
//                                RESULT_DIR + "ClusterToGOBP.txt");
//         The following is for GO CC
        Map<String, Set<String>> protein2GOC = fu.loadSetMap(goDir + "ProteinIdToGOC.txt");
        Set<String> totalIds = new HashSet<String>();
        List<ClusterInfo> clusters = loadClusters();
        for (ClusterInfo c : clusters) {
            totalIds.addAll(c.proteins);
        }
        for (Iterator<String> it = protein2GOC.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            if (!totalIds.contains(id))
                it.remove();
        }
        Map<String, Integer> proteinsInGOC = countGOUsage(protein2GOC);
//        exportProteinNumberMap(fu, 
//                               RESULT_DIR + "ProteinNumberInGOCC.txt",
//                               proteinsInGOC);
        assignGOTermsToClusters(clusters, 
                                protein2GOC, 
                                proteinsInGOC, 
                                RESULT_DIR + "ClusterToGOCC.txt");
    }

    private void assignGOTermsToClusters(List<ClusterInfo> clusters, 
                                         Map<String, Set<String>> protein2GO,
                                         Map<String, Integer> proteinsInGO,
                                         String output) throws IOException {
        // Count how many proteins using term for each cluster
        // First for MF
        for (ClusterInfo c : clusters) {
            List<String> goList = new ArrayList<String>();
            c.topics = goList;
            for (String p : c.proteins) {
                Set<String> goIds = protein2GO.get(p);
                // Some proteins may have no GO annotations
                if (goIds != null)
                    goList.addAll(goIds);
            }
            // Want to convert the list to map
            c.proteinsInTopics = countTopics(goList);
        }
        // Now p-values can be calculated
        // This is the file for outputting
        int total = TOTAL_PROTEIN;
        FileUtility outputFu = new FileUtility();
        outputFu.setOutput(output);
        outputFu.printLine("ClusterIndex\tProteinsInGO\tP-value\tGOInCluster");
        int proteinInCluster;
        int proteinInGO;
        int proteinInBoth;
        int cluster;
        double p;
        int index = 0;
        for (ClusterInfo c : clusters) {
            index ++;
            final Map<String, Integer> map = c.proteinsInTopics;
            List<String> keys = new ArrayList<String>(map.keySet());
            Collections.sort(keys,
                             new Comparator<String>() {
                public int compare(String s1, String s2) {
                    int c1 = map.get(s1);
                    int c2 = map.get(s2);
                    return c2 - c1;
                }
            });
            // Want to sort based on number
            for (String goId : keys) {
                proteinInBoth = map.get(goId);
                proteinInCluster = c.proteinNumber;
                proteinInGO = proteinsInGO.get(goId);
                p = calculatePValue(total, 
                                    proteinInBoth, 
                                    proteinInCluster,
                                    proteinInGO);
                outputFu.printLine(index + "\t" + 
                                   proteinInBoth + "\t" + 
                                   p + "\t" + 
                                   goId);
            }
        }
        outputFu.close();
    }
    
    private Map<String, Integer> countGOUsage(Map<String, Set<String>> protein2GO) {
        Map<String, Integer> proteinsInGO = new HashMap<String, Integer>();
        for (Iterator<String> it = protein2GO.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            Set<String> goIds = protein2GO.get(id);
            for (String goId : goIds) {   
                Integer c = proteinsInGO.get(goId);
                if (c == null)
                    proteinsInGO.put(goId, 1);
                else
                    proteinsInGO.put(goId, ++c);
            }
        }
        return proteinsInGO;
    }
    
    public void countProteinsInClusters() throws IOException {
        String inputFileName = RESULT_DIR + CLUSTERING_RESULT_FILE;
        FileUtility fu = new FileUtility();
        fu.setInput(inputFileName);
        FileUtility outputFu = new FileUtility();
        //String outputFileName = RESULT_DIR + "ProteinsInClustersI24C26.txt";
        //String outputFileName = RESULT_DIR + "ProteinsInClustersWithComplexAndSet.txt";
        String outputFileName = RESULT_DIR + "ProteinsInMCLClusters.txt";
        outputFu.setOutput(outputFileName);
        // To use FileUtility
        String line = null;
        String[] tokens = null;
        int c = 1;
        outputFu.printLine("Cluster\tProtein");
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            outputFu.printLine(c + "\t" + tokens.length);
            c ++;
        }
        fu.close();
        outputFu.close();
    }
    
    public void generateInteractionsForCluster() throws IOException {
        // Load all interactions
        Set<String> interactions = fu.loadInteractions(RESULT_DIR + "FIInteractions67.txt");
        List<ClusterInfo> clusters = loadClusters();
        // Just want to try the first cluster
        ClusterInfo firstCluster = clusters.get(0);
        List<String> clusterProteins = firstCluster.proteins;
        System.out.println("Total Proteins in Cluster 1: " + clusterProteins.size());
        FileUtility outputFu = new FileUtility();
        String outFileName = RESULT_DIR + "InteractionsInCluster1.txt";
        outputFu.setOutput(outFileName);
        for (String interaction : interactions) {
            String[] tokens = interaction.split(" ");
            if (clusterProteins.contains(tokens[0]) &&
                clusterProteins.contains(tokens[1]))
                outputFu.printLine(interaction);
        }
        outputFu.close();
    }
    
    /**
     * This method is used to calculate cluqueness to cluster mapping information.
     * @throws IOException
     */
    public void generateCliquenessToCluster() throws IOException {
        Set<String> allFIs = fu.loadInteractions(R3Constants.INTERACTION_FILE_NAME);
        List<ClusterInfo> clusters = loadClusters();
        CliquenessAnalyzer cliquenessAnalyzer = new CliquenessAnalyzer();
        int index = 0;
        String outFileName = RESULT_DIR + "ClusterToCliqueness041208.txt";
        fu.setOutput(outFileName);
        fu.printLine("Cluster\tProteinsInCluster\tCliqueness");
        for (ClusterInfo cluster : clusters) {
            Set<String> proteins = new HashSet<String>(cluster.proteins);
            if (proteins.size() < 6)
                break; // Do a cutoff
            // Search for interactions
            Set<String> interactions = new HashSet<String>();
            for (String interaction : allFIs) {
                index = interaction.indexOf(" ");
                String id1 = interaction.substring(0, index);
                String id2 = interaction.substring(index + 1);
                if (proteins.contains(id1) &&
                    proteins.contains(id2))
                    interactions.add(interaction);
            }
            // Do a quick check
            Set<String> idsInInteractions = InteractionUtilities.grepIDsFromInteractions(interactions);
            if (idsInInteractions.size() != proteins.size()) {
                proteins.removeAll(idsInInteractions);
                //throw new IllegalStateException(cluster.index + " has proteins not in FIs: " + proteins);
                System.out.println(cluster.index + " has proteins not in FIs: " + proteins);
            }
            // Note: There are some proteins in the clusters that are not connected to others. This
            // is very weird.
            double cliqueness = cliquenessAnalyzer.calculateAverageCliqueness(interactions, 
                                                                              idsInInteractions);
            fu.printLine(cluster.index + "\t" + cluster.proteinNumber + "\t" + cliqueness);
        }
        fu.close();
    }
    
    /**
     * Get the sif files for Cytoscape for the first ten concentrated clusters.
     * @throws IOException
     */
    public void generateSifFiles() throws IOException {
        // Load the first ten cluster index
        List<Integer> clusterIndex = new ArrayList<Integer>();
        FileUtility fu = new FileUtility();
        String fileName = RESULT_DIR + "I24C26ClusteringInfoSorted.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        int c = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            clusterIndex.add(Integer.parseInt(tokens[0]));
            c ++;
            if (c == 10)
                break;
        }
        fu.close();
        fileName = RESULT_DIR + "I24C26Clustering.txt";
        fu.setInput(fileName);
        List<String> clusters = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            clusters.add(line);
        }
        fu.close();
        Set<String> interactions = fu.loadInteractions(RESULT_DIR + "FIInteractions06.txt");
        c = 1;
        int index = 0;
        String protein1, protein2;
        for (Integer i : clusterIndex) {
            System.out.println("Handling cluster " + i);
            fileName = RESULT_DIR + "InteractionsInCluster" + i + ".sif";
            fu.setOutput(fileName);
            line = clusters.get(i);
            String[] proteins = line.split("\t");
            List<String> queryList = Arrays.asList(proteins);
            // Check interactions
            for (String inter : interactions) {
                index = inter.indexOf(" ");
                protein1 = inter.substring(0, index);
                protein2 = inter.substring(index + 2);
                if (queryList.contains(protein1) && 
                    queryList.contains(protein2)) {
                    fu.printLine(protein1 + " pp " + protein2);
                }
            }
            fu.close();
        }
    }
    
    public List<ClusterInfo> loadClusters() throws IOException {
        FileUtility fu = new FileUtility();
//        Set<String> interactions = fu.loadInteractions(RESULT_DIR + "FIInteractions06.txt");
//        String fiFile = RESULT_DIR + "I24C26Clustering.txt";
        Set<String> interactions = fu.loadInteractions(R3Constants.INTERACTION_FILE_NAME);
        //Set<String> interactions = fu.loadInteractions("results/v2/FIInteractionsWithComplexAndSet.txt");
        String fiFile = RESULT_DIR + CLUSTERING_RESULT_FILE;
        fu.setInput(fiFile);
        List<ClusterInfo> clusterInfos = new ArrayList<ClusterInfo>();
        String line = null;
        int index = 1;
        String protein1, protein2;
        int c = 0;
        while ((line = fu.readLine()) != null) {
            ClusterInfo info = new ClusterInfo();
            clusterInfos.add(info);
            info.index = index;
            index ++;
            String[] proteins = line.split("\t");
            info.proteins = Arrays.asList(proteins);
            if (proteins.length == 1) {
                info.interactionNumber = 0;
                info.proteinNumber = 1;
                info.interactionPerProtein = 0.0d;
                continue;
            }
            info.proteinNumber = proteins.length;
            Arrays.sort(proteins);
            // Check interactions
            c = 0;
            for (int i = 0; i < proteins.length - 1; i++) {
                for (int j = i + 1; j < proteins.length; j++) {
                    String pair = proteins[i] + " " + proteins[j];
                    if (interactions.contains(pair))
                        c ++;
                }
            }
            info.interactionNumber = c;
            info.interactionPerProtein = (double) c / proteins.length;
        }
        fu.close();
        return clusterInfos;
    }
    
    /**
     * Load MCL clustering results from the specified file name.
     * @param fileName
     * @return
     * @throws IOException
     */
    public List<Set<String>> loadMCLClusters(String fileName) throws IOException {
        List<Set<String>> clusters = new ArrayList<Set<String>>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            Set<String> set = new HashSet<String>();
            String[] tokens = line.split("\t");
            for (String gene : tokens)
                set.add(gene);
            clusters.add(set);
        }
        fu.close();
        return clusters;
    }
    
    public void generateClusterInfo() throws IOException {
        FileUtility fu = new FileUtility();
        List<ClusterInfo> clusterInfos = loadClusters();
        // Sorting based on the values of interactionPerProtein
        Collections.sort(clusterInfos, new Comparator<ClusterInfo>() {
            public int compare(ClusterInfo info1, ClusterInfo info2) {
                Double value1 = info1.interactionPerProtein;
                Double value2 = info2.interactionPerProtein;
                return value2.compareTo(value1); // In the reverse order
            }
        });
        fu.setOutput(RESULT_DIR + "I24C26ClusteringInfoSorted.txt");
        int sortedIndex = 1;
        // Header
        fu.printLine("Cluster\tProtein\tInteraction\tInterPerProtein\tSortedIndex");
        for (ClusterInfo info : clusterInfos) {
            fu.printLine(info.index + "\t" + info.proteinNumber + "\t" + 
                         info.interactionNumber + "\t" + info.interactionPerProtein + 
                         "\t" + sortedIndex);
            sortedIndex ++;
        }
        fu.close();
    }
    
    public void generateTopicAndClusterContingencyTable() throws IOException {
        String fileName = RESULT_DIR + "ClusterToTopic.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine(); // Escape the header
        // Used to hold list of ClusterInfo for each pathway
        Map<String, List<ClusterInfo>> topicToCluster = new HashMap<String, List<ClusterInfo>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Check the topic first
            List<ClusterInfo> infos = topicToCluster.get(tokens[2]);
            if (infos == null) {
                infos = new ArrayList<ClusterInfo>();
                topicToCluster.put(tokens[2], infos);
            }
            ClusterInfo info = new ClusterInfo();
            info.index = Integer.parseInt(tokens[0]);
            info.proteinNumber = Integer.parseInt(tokens[1]);
            infos.add(info);
        }
        fu.close();
        int size = 84; // there are 84 pathways having no less than 50 proteins
        // Get the list of pathways to be printed out
        List<String> pathways = new ArrayList<String>();
        fileName = RESULT_DIR + "ProteinNumberInTopics.txt";
        fu.setInput(fileName);
        line = fu.readLine();
        int c = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            pathways.add(tokens[2]);
            c ++;
            if (c == size)
                break;
        }
        fu.close();
//      Print out a table: columns for clusters, rows for pathways
        fileName = RESULT_DIR + "TopicClusterContingencyTable84.txt";
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        int[] numberInClusers = new int[size];
        for (String t : pathways) {
            System.out.println(t);
            List<ClusterInfo> infos = topicToCluster.get(t);
            // infos should be sorted based on the input file, ClusterToTopic.txt.
            for (ClusterInfo info : infos) {
                if (info.index > size)
                    continue; // Not in the first 84 clusters
                numberInClusers[info.index - 1] = info.proteinNumber;
            }
            // output the array
            for (int i : numberInClusers) 
                builder.append(i + "\t");
            fu.printLine(builder.toString());
            // reset the array
            for (int i = 0; i < size; i++)
                numberInClusers[i] = 0;
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * This internal private class is used to hold cluster information
     */
    static class ClusterInfo {
        int index;
        int proteinNumber;
        int interactionNumber;
        double interactionPerProtein;
        int sortedIndex;
        List<String> proteins;
        List<String> topics;
        Map<String, Integer> proteinsInTopics;
        
        public ClusterInfo() {
        }
        
        public void addTopic(String topic) {
            if (topics == null)
                topics = new ArrayList<String>();
            topics.add(topic);
        }
        
    }
    
    /**
     * This internal class is used to describe a topic (pathway or go) in a cluster
     */
    private class TopicInClusterInfo {
        private int clusterIndex;
        private String id;
        private double pvalue;
        private int numbder;
        
        public TopicInClusterInfo() {
            
        }
    }
}
