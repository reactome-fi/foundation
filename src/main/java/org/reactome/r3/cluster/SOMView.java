/*
 * Created on Jul 13, 2007
 *
 */
package org.reactome.r3.cluster;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.Point;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.ToolTipManager;

import org.gk.util.GKApplicationUtilities;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to display SOM results in a 2D rectangle.
 * @author guanming
 *
 */
public class SOMView extends JFrame {
    
    public SOMView() {
    }
    
    public void displaySOM(SOMWeightNodes nodes) {
        // Just don't want the tooltip dismissed
        ToolTipManager.sharedInstance().setDismissDelay(100000);
        List<ReferenceNode> nodeList = nodes.getReferenceNodes();
        // Get the maxX and maxY
        int maxX = 0;
        int maxY = 0;
        for (ReferenceNode node : nodeList) {
            Point p = node.getLocation();
            if (p.x > maxX)
                maxX = p.x;
            if (p.y > maxY)
                maxY = p.y;
        }
        // Use JLabel for the time being
        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(maxY, maxX));
        // This map is used to color labels
        Map<JLabel, Double> labelToDist = new HashMap<JLabel, Double>();
        for (int y = 0; y < maxY; y ++) {
            for (int x = 0; x < maxX; x++) {
                ReferenceNode node = nodes.getNodeAt(x, y);
                double distWithNeighbors = nodes.calculateDistWithNeighbors(node);
                JLabel label = createReferenceLabel();
                labelToDist.put(label, distWithNeighbors);
                panel.add(label);
                if (node.getLabels() != null) {
                    Set<String> pathwayNames = node.getLabels();
                    StringBuilder builder = new StringBuilder();
                    builder.append("<html>");
                    for (Iterator<String> it = pathwayNames.iterator();
                         it.hasNext();) {
                        builder.append(it.next());
                        if (it.hasNext())
                            builder.append("<br>");
                    }
                    builder.append("</html>");
                    label.setText(builder.toString());
                    label.setToolTipText(label.getText());
                }
            }
        }
        showColors(labelToDist);
        getContentPane().add(panel, BorderLayout.CENTER);
        validate();
    }
    
    private void showColors(Map<JLabel, Double> labelToDist) {
        // The the range of distance
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (Iterator<JLabel> it = labelToDist.keySet().iterator(); it.hasNext();) {
            JLabel label = it.next();
            double dist = labelToDist.get(label);
            if (dist < min)
                min = dist;
            if (dist > max)
                max = dist;
        }
        // Use blue to indicate the color
        double diff = max - min;
        for (Iterator<JLabel> it = labelToDist.keySet().iterator(); it.hasNext();) {
            JLabel label = it.next();
            double dist = labelToDist.get(label);
            float colorValue = (float) ((dist - min) / diff);
            Color c = new Color(colorValue, colorValue, colorValue);
            label.setOpaque(true);
            label.setBackground(c);
        }
    }
    
    private JLabel createReferenceLabel() {
        JLabel label = new JLabel();
        label.setBorder(BorderFactory.createEtchedBorder());
        return label;
    }
    
    public SOMWeightNodes loadSOM(String soFile) throws IOException {
        // Load SOM results into a SOMWeightNodes that contain
        // a list of ReferenceNode
        // Here SOMWeightNodes are used as a model actually
        SOMWeightNodes somNodes = new SOMWeightNodes();
        FileUtility fu = new FileUtility();
        fu.setInput(soFile);
        String line = null;
        boolean isInResult = false;
        ReferenceNode node = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("SOM results:")) {
                isInResult = true;
            }
            else if (isInResult) {
//                location:14, 23
//                index:0.0
//                label:Allantoin degradation(P)
//                vector:0.0,0.0,0.0,
                if (line.startsWith("location:")) {
                    // Keep the old one
                    if (node != null)
                        somNodes.addNode(node);
                    node = new ReferenceNode();
                    Point location = extractLocation(line);
                    node.setLocation(location);
                }
                else if (line.startsWith("index:")) {
                    node.setIndex(extractIndex(line));
                }
                else if (line.startsWith("label:")) {
                    extractLabels(line, node);
                }
                else if (line.startsWith("vector:")) {
                    extractVector(line, node);
                }
            }
        }
        if (node != null) // Last node
            somNodes.addNode(node);
        fu.close();
        return somNodes;
    }
    
    private void extractVector(String line,
                               ReferenceNode node) {
        int index = line.indexOf(":");
        String sub = line.substring(index + 1,
                                    line.length() - 1); // Remove the last ","
        String[] tokens = sub.split(",");
        double[] vector = new double[tokens.length];
        for (int i = 0; i < tokens.length; i++) {
            vector[i] = Double.parseDouble(tokens[i]);
        }
        node.setReferenceVector(vector);
    }
    
    private void extractLabels(String line, 
                               ReferenceNode node) {
        int index = line.indexOf(":");
        String sub = line.substring(index + 1);
        if (sub.length() == 0)
            return;
        String[] tokens = sub.split("(\\), )");
        for (String token : tokens) {
            if (token.endsWith(")"))
                node.addLabel(token);
            else
                node.addLabel(token + ")");
        }
    }
    
    private Point extractLocation(String line) {
        int index = line.indexOf(":");
        int index1 = line.indexOf(",");
        int x = Integer.parseInt(line.substring(index + 1, index1));
        int y = Integer.parseInt(line.substring(index1 + 1).trim());
        return new Point(x, y);
    }
    
    private double extractIndex(String line) {
        int index = line.indexOf(":");
        return Double.parseDouble(line.substring(index + 1));
    }
    
    public static void main(String[] args) {
        SOMView view = new SOMView();
        view.setSize(600, 600);
        GKApplicationUtilities.center(view);
        view.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        view.setVisible(true);
        String fileName = R3Constants.RESULT_DIR + "SOM_24_24_200.txt";
        //String fileName = R3Constants.RESULT_DIR + "SOM_40_30_100_ORI.txt";
        try {
            SOMWeightNodes nodes = view.loadSOM(fileName);
            view.displaySOM(nodes);
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }
    
}
