/*
 * Created on Apr 23, 2007
 *
 */
package org.reactome.r3.cluster;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.imageio.ImageIO;

import org.reactome.r3.util.FileUtility;

/**
 * A simple hierarchical clustering to merge pathways based on FI or protein ids.
 * @author guanming
 *
 */
public class HierarchicalCluster {
    /**
     * The methods used to hierarchical clustering. Right now, only three methods are used.
     */
    public enum ClusterDistanceMethod {
        SINGLE,
        COMPLETE,
        AVERAGE
    }
    
    private DistanceCalculator distanceCalculator;
    // Default using average
    private ClusterDistanceMethod method = ClusterDistanceMethod.AVERAGE;
    
    private final boolean debug = true;
    
    public HierarchicalCluster() {
    }
    
    public void setDistanceCalculator(DistanceCalculator calculator) {
        this.distanceCalculator = calculator;
    }
    
    public void setMethod(ClusterDistanceMethod method) {
        this.method = method;
    }
    
    public ClusterDistanceMethod getMethod() {
        return this.method;
    }
    
    public void grepAllClusters(HierarchicalClusterNode newFirstNode,
                                Set<HierarchicalClusterNode> allNodes) {
        Set<HierarchicalClusterNode> current = new HashSet<HierarchicalClusterNode>();
        current.add(newFirstNode);
        Set<HierarchicalClusterNode> next = new HashSet<HierarchicalClusterNode>();
        while (current.size() > 0) {
            for (HierarchicalClusterNode tmp : current) {
                allNodes.add(tmp);
                if (tmp.childNode1 == null)
                    continue;
                next.add(tmp.childNode1);
                next.add(tmp.childNode2);
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    /**
     * To run a hierarchical clustering, the client to this class should call this method 
     * by passing a collection of ids
     * @param ids an initial collection of ids to be clustered.
     * @return the top-most cluster node.
     */
    public HierarchicalClusterNode cluster(Collection<String> ids) {
        // Make sure DistanceCalculator should be assigned
        if (distanceCalculator == null)
            throw new IllegalStateException("DistanceCalculator should not be null!");
        // Create a list of cluster nodes from ids
        List<HierarchicalClusterNode> nodes = new ArrayList<HierarchicalClusterNode>();
        for (String id : ids) {
            HierarchicalClusterNode node = new HierarchicalClusterNode();
            node.ids = new HashSet<String>();
            node.ids.add(id);
            nodes.add(node);
        }

        double min = Double.MAX_VALUE;
        while (nodes.size() > 1) {
            // Search the pair with the minimum shortest path
            HierarchicalClusterNode minNode1 = null;
            HierarchicalClusterNode minNode2 = null;
            for (int i = 0; i < nodes.size() - 1; i++) {
                HierarchicalClusterNode node1 = nodes.get(i);
                for (int j = i + 1; j < nodes.size(); j++) {
                    HierarchicalClusterNode node2 = nodes.get(j);
                    double dist = calculateDistance(node1, node2);
                    if (dist < min) {
                        minNode1 = node1;
                        minNode2 = node2;
                        min = dist;
                    }
                    else if (dist == min && minNode1 != null && minNode2 != null) { // Try to have as many ids in the same clusters
                        // Get the highest possible combination: this selection will
                        // generate more leaf nodes.
                        if (node1.ids.size() > minNode1.ids.size() ||
                            node2.ids.size() > minNode2.ids.size() ||
                            node1.ids.size() > minNode2.ids.size() ||
                            node2.ids.size() > minNode1.ids.size()) {
                            minNode1 = node1;
                            minNode2 = node2;
                        }
                    }
                }
            }
            // Find the pair
            HierarchicalClusterNode newCluster = new HierarchicalClusterNode();
            newCluster.ids = new HashSet<String>();
            newCluster.ids.addAll(minNode1.ids);
            newCluster.ids.addAll(minNode2.ids);
            newCluster.setChildNode1(minNode1);
            newCluster.setChildNode2(minNode2);
            newCluster.pathDistance = min;
            nodes.remove(minNode1);
            nodes.remove(minNode2);
            nodes.add(newCluster);
            min = Double.MAX_VALUE;
        }
        return nodes.get(0); // There should be only one node in the list.
    }
    
    private double calculateDistance(HierarchicalClusterNode node1,
                                     HierarchicalClusterNode node2) {
        switch (method) {
            case AVERAGE :
                return calculateAverageDistance(node1, node2);
            case SINGLE :
                return calculateSingleDistance(node1, node2);
            case COMPLETE :
                return calculateCompleteDistance(node1, node2);
        }
        throw new IllegalStateException("Unknow cluster distance method: " + method);
    }
    
    private double calculateAverageDistance(HierarchicalClusterNode node1,
                                            HierarchicalClusterNode node2) {
        double total = 0.0;
        double count = 0;
        for (String id1 : node1.ids) {
            for (String id2 : node2.ids) {
                total += distanceCalculator.calculateDistance(id1, id2);
                count ++;
            }
        }
        return total / count;
    }
    
    private double calculateSingleDistance(HierarchicalClusterNode node1,
                                           HierarchicalClusterNode node2) {
        double rtn = Double.MAX_VALUE;
        for (String id1 : node1.ids) {
            for (String id2 : node2.ids) {
                double dist = distanceCalculator.calculateDistance(id1, id2);
                if (dist < rtn)
                    rtn = dist;
            }
        }
        return rtn;
    }
    
    private double calculateCompleteDistance(HierarchicalClusterNode node1, 
                                             HierarchicalClusterNode node2) {
        double rtn = Double.MIN_VALUE;
        for (String id1 : node1.ids) {
            for (String id2 : node2.ids) {
                double dist = distanceCalculator.calculateDistance(id1, id2);
                if (dist > rtn)
                    rtn = dist;
            }
        }
        return rtn;
    }
    
    public void drawDendrogram(List<HierarchicalClusterNode> clusters,
                              int width,
                              String fileName) throws Exception {
        drawDendrogram(clusters, width, null, null, fileName);
    }
    
    public void drawDendrogram(List<HierarchicalClusterNode> clusters, 
                               int width,
                               double[] lines,
                               List<String> highlights,
                               String fileName) throws Exception {
        int padding = 2;
        // List all single nodes first
        List<HierarchicalClusterNode> leaves = new ArrayList<HierarchicalClusterNode>();
        for (HierarchicalClusterNode node : clusters) {
            if (node.ids.size() == 1) {
                leaves.add(node);
            }
        }
        System.out.println("total singles: " + leaves.size());
        // Used to get graphic context
        BufferedImage image = new BufferedImage(5, 5, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D g2 = (Graphics2D) image.getGraphics();
        Font font = new Font("Monospaced", 
                             Font.PLAIN, 
                             12);
        g2.setFont(font);
        FontMetrics metrics = g2.getFontMetrics();
        String sample = leaves.get(0).ids.iterator().next();
        // Need to get the widest string
        Rectangle2D stringBounds = null;
        for (HierarchicalClusterNode node : leaves) {
            String id = node.ids.iterator().next();
            Rectangle2D tmp = metrics.getStringBounds(id, g2);
            if (stringBounds == null ||
                tmp.getWidth() > stringBounds.getWidth())
                stringBounds = tmp;
        }
        // Calculate height
        int height = (int) (stringBounds.getHeight() + 2) * leaves.size();
        image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
        g2 = (Graphics2D) image.getGraphics();
        g2.setFont(font);
        g2.setPaint(Color.WHITE);
        g2.fillRect(0, 0, image.getWidth(), image.getHeight());
        int x = (int) (width - stringBounds.getWidth() - padding);
        int y = (int) (stringBounds.getHeight());
        g2.setPaint(Color.BLACK);
        //Map<String, Integer> sampleToTimeSpane = new NatureGBMAnalyzer().loadSampleToTimeSpan();
        //int c1 = 0; int c2 = 0;
//        Map<String, Integer> sampleToSurvival = new NatureGBMAnalyzer().loadSampleToSurvivalRate();
//        Map<String, String> sampleToRecurrent = new NatureGBMAnalyzer().loadSampleToRecurrent();
//       // Map<String, List<Integer>> sampleToNetworkClusters = new NatureGBMAnalyzer().loadSampleToNetworkClusters();
//        // Min and max
//        int min = Integer.MAX_VALUE;
//        int max = Integer.MIN_VALUE;
//        for (String s : sampleToSurvival.keySet()) {
//            Integer r = sampleToSurvival.get(s);
//            if (r > max)
//                max = r;
//            if (r < min)
//                min = r;
//        }
//        max = 1600;
        if (highlights == null)
            highlights = new ArrayList<String>(); // For easy checking
        for (HierarchicalClusterNode single : leaves) {
            String id = single.ids.iterator().next();
//            List<Integer> networkClusters = sampleToNetworkClusters.get(id);
//            String output = null;
//            if (networkClusters.toString().equals("[0, 1]"))
//                output = "1";
//            else
//                output = "2";
//            System.out.println(id + "\t" + networkClusters + "\t" + output);
//            String recurrent = sampleToRecurrent.get(id);
//            if (recurrent.equals("No"))
//                g2.setPaint(Color.green);
//            else
//                g2.setPaint(Color.red);
//            Integer span = sampleToTimeSpane.get(id);
//            Integer survial = sampleToSurvival.get(id);
//            if (survial > max)
//                survial = max;
//            // Need to adjust shaded for green
//            int c1 = (int)((double)(survial - min) / (max - min) * 255);
//            System.out.println(c1);
//            Color c = new Color(0, c1, 0);
//            g2.setPaint(c);
//            if (span < 2) {
//                g2.setPaint(Color.red);
//                c1 ++;
//            }
//            else {
//                g2.setPaint(Color.green);
//                c2 ++;
//            }
//            g2.fillRect(x, (int) (y - stringBounds.getHeight()), 
//                        (int)stringBounds.getWidth(), (int)stringBounds.getHeight() + 2);
            if (highlights.contains(id))
                g2.setPaint(Color.BLUE);
            else
                g2.setPaint(Color.BLACK);
            g2.drawString(id, x, y);
            single.x = x;
            single.y = (int)(y - stringBounds.getHeight() / 2);
            y += (stringBounds.getHeight() + 2);
            //System.out.println(id + "\t" + sampleToTimeSpane.get(id));
        }
        HierarchicalClusterNode root = clusters.get(0);
        // Need to map distance for pixel
        double ratio = (width - stringBounds.getWidth() - 2 * padding) / root.pathDistance;
        drawClusterNodes(root, 
                         ratio, 
                         x,
                         g2);
        // Draw some lines
        if (lines != null) {
            g2.setPaint(Color.BLUE);
            for (double cutoff : lines) {
                int x1 = (int) (x - ratio * cutoff);
                g2.drawLine(x1, 
                            0, 
                            x1,
                            image.getHeight());
            }
        }
        ImageIO.write(image, "png", new File(fileName));
    }
    
    private void drawClusterNodes(HierarchicalClusterNode cluster,
                                  double ratio,
                                  int x0,
                                  Graphics2D g2) {
        if (cluster.ids.size() == 1)
            return ; // Don't need to draw
        drawClusterNodes(cluster.childNode1, ratio, x0, g2);
        drawClusterNodes(cluster.childNode2, ratio, x0, g2);
        // If both children have been draw, draw this node
        // Otherwise go to children
        if (cluster.childNode1.x > 0 && cluster.childNode2.x > 0) {
            // Draw this cluster
            cluster.x = (int) (x0 - ratio * cluster.pathDistance);
            cluster.y = (cluster.childNode1.y + cluster.childNode2.y) / 2;
            g2.drawLine(cluster.childNode1.x, cluster.childNode1.y, 
                        cluster.x, cluster.childNode1.y);
            g2.drawLine(cluster.childNode2.x, cluster.childNode2.y, 
                        cluster.x, cluster.childNode2.y);
            g2.drawLine(cluster.x, cluster.childNode1.y,
                        cluster.x, cluster.childNode2.y);
        }
    }
    
    public void outputCluster(String space,
                              HierarchicalClusterNode cluster) {
        if (cluster.ids.size() == 1)
            System.out.println(space + cluster.ids);
        else
            System.out.println(space + cluster.pathDistance + ": " + cluster.ids + "(" + cluster.ids.size() + ")");
        if (cluster.childNode1 != null)
            outputCluster(space + "  ", cluster.childNode1);
        if (cluster.childNode2 != null)
            outputCluster(space + "  ", cluster.childNode2);
    }
    
    /**
     * Get all nodes from the root, and returned in a sorted list.
     * @param root
     * @return
     */
    public List<HierarchicalClusterNode> grepAllNodesInDFS(HierarchicalClusterNode root) {
        List<HierarchicalClusterNode> rtn = new ArrayList<HierarchicalClusterNode>();
        grepAllNodesInDFS(root, rtn);
        return rtn;
    }
    
    private void grepAllNodesInDFS(HierarchicalClusterNode node, List<HierarchicalClusterNode> list) {
        list.add(node);
        if (node.childNode1 != null)
            grepAllNodesInDFS(node.childNode1, list);
        if (node.childNode2 != null)
            grepAllNodesInDFS(node.childNode2, list);
    }
    
//    @Test
//    public void runDrawDendrogram() throws Exception {
//        //String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "SampleClusteringFromMutationsShortestPath.txt";
//        //String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "SampleClustersAvgShortestPath070809.txt";
//        //String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "SampleClusteringFromHomoCNVAndMutation_063009.txt";
//        //String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "TCGASampleClustersFromNetworkClusters_Weighted.txt";
//        //String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "HierarchicalClusterResultsForParsons100909.txt";
//        String fileName = CancerResequenceDataSetAnalyzer.OVARIAN_DIR_NAME + "HierarchicalClusterOfPlatCorGenes.txt";
//        List<HierarchicalClusterNode> clusters = loadHierarchicalClusters(fileName);
//        //outputCluster("", clusters.get(0));
//        drawDendrogram(clusters, 
//                       500,
//                       //NatureGBMAnalyzer.TCGA_GBM_DIR + "TCGASampleClustersFromNetworkClusters_Recurrent_Weighted_Survival.png");
//                       //NatureGBMAnalyzer.TCGA_GBM_DIR + "HierarchicalClusterResultsForParsons100909.png");
//                       CancerResequenceDataSetAnalyzer.OVARIAN_DIR_NAME + "HierarchicalClusterOfPlatCorGenes.png");
//    }
    
    /**
     * Draw a diagram for sample information so that it can be displayed with heat map
     * from R.
     * @param sampleToInfo
     * @param infoList
     * @param colors
     * @throws IOException
     */
    public void drawSampleInformation(Map<String, String> sampleToInfo,
                                      List<String> samples,
                                      List<String> infoList,
                                      List<Color> colors,
                                      int length,
                                      int width,
                                      String outFileName) throws IOException {
        double leg = (double) length / samples.size();
        BufferedImage image = new BufferedImage(width, 
                                                length, 
                                                BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D g2 = (Graphics2D) image.getGraphics();
         
        Rectangle2D rect = new Rectangle2D.Double();
        int index = 0;
        // Use default light gray
        g2.setPaint(Color.lightGray);
        for (String sample : samples) {
            double y = leg * index;
            Color oldColor = g2.getColor();
            rect.setRect(0.0d,
                         y,
                         width,
                         leg);
            String info = sampleToInfo.get(sample);
            int colorIndex = infoList.indexOf(info);
            Color color = null;
            if (colorIndex >= 0)
                color = colors.get(colorIndex);
            if (color != null)
                g2.setPaint(color);
            g2.fill(rect);
            g2.draw(rect);
            g2.setPaint(oldColor);
            index ++;
        }
        ImageIO.write(image, 
                      "png",
                      new File(outFileName));
    }
    
    public List<HierarchicalClusterNode> loadHierarchicalClusters(String fileName, double cutoff) throws IOException {
        List<HierarchicalClusterNode> clusters = loadHierarchicalClusters(fileName);
        List<HierarchicalClusterNode> rtn = new ArrayList<HierarchicalClusterNode>();
        readClusters(clusters.get(0), cutoff, rtn);
        return rtn;
    }
    
    public void sortClustersBasedOnSizes(List<HierarchicalClusterNode> clusters) {
        Collections.sort(clusters, new Comparator<HierarchicalClusterNode>() {
            public int compare(HierarchicalClusterNode cluster1, HierarchicalClusterNode cluster2) {
                return cluster2.getIds().size() - cluster1.getIds().size();
            }
        });
    }
    
    public void readClusters(HierarchicalClusterNode cluster, 
                             double cutoff, 
                             List<HierarchicalClusterNode> list) {
        if (cluster.pathDistance <= cutoff)
            list.add(cluster);
        else {
            readClusters(cluster.childNode1, cutoff, list);
            readClusters(cluster.childNode2, cutoff, list);
        }
    }
    
    public List<HierarchicalClusterNode> loadHierarchicalClusters(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        List<HierarchicalClusterNode> clusterNodes = new ArrayList<HierarchicalClusterNode>();
        boolean inData = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("hierarchical layouting...")) {
                inData = true;
                continue;
            }
            if (!inData)
                continue;
            // Start the parsing
            line = line.trim();
            HierarchicalClusterNode cluster = new HierarchicalClusterNode();
            Set<String> ids = extractIDs(line);
            cluster.ids = ids;
            if (!line.startsWith("[")) {
                // Get distance
                int index = line.indexOf(":");
                cluster.pathDistance = new Double(line.substring(0, index));
            }
            // attach to the parent
            if (clusterNodes.size() > 0) {
                for (int i = clusterNodes.size() - 1; i >= 0; i--) {
                    HierarchicalClusterNode parentNode = clusterNodes.get(i);
                    // Check if it is a parent of the current node
                    if (isParentNode(parentNode, cluster)) {
                        if (parentNode.childNode1 == null) {
                            parentNode.childNode1 = cluster;
                            break;
                        }
                        else if (parentNode.childNode2 == null) {
                            parentNode.childNode2 = cluster;
                            break;
                        }
                    }
                }
            }
            clusterNodes.add(cluster);
        }
        fu.close();
        return clusterNodes;
    }
    
    private boolean isParentNode(HierarchicalClusterNode parent, 
                                 HierarchicalClusterNode current) {
        Set<String> copy = new HashSet<String>(current.ids);
        copy.removeAll(parent.ids);
        if (copy.size() == 0)
            return true;
        return false;
    }
    
    private Set<String> extractIDs(String line) {
        int index1 = line.indexOf("[");
        int index2 = line.lastIndexOf("]");
        line = line.substring(index1 + 1, index2);
        String[] tokens = line.split(", ");
        Set<String> rtn = new HashSet<String>();
        for (String token : tokens)
            rtn.add(token);
        return rtn;
    }

}
