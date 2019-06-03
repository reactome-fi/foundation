/*
 * Created on Jan 11, 2008
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to calculate cliqueness of the functional interaction network.
 * @author guanming
 *
 */
public class CliquenessAnalyzer {
    
    public CliquenessAnalyzer() {
    }
    
    @Test
    public void calculateAverageCliquenessFromRandomGraph() throws IOException {
        FileUtility fu = new FileUtility();
        //String intFileName = R3Constants.INTERACTION_FILE_NAME;
        //String intFileName = R3Constants.RESULT_DIR + "FIInteractions60_121707_BigComp.txt";
        String intFileName = R3Constants.RESULT_DIR + ShortestPathAnalyzer.BIG_COMP_INT_FILE;
        Set<String> interactions = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        for (int i = 0; i < 10; i++) {
            System.out.println("Randome graph sampling " + i + "...");
            Graph<String, DefaultEdge> graph = graphAnalyzer.createRandomGraph(ids, interactions.size());
            System.out.println("Graph: " + graph.vertexSet().size() + " vertices, " + 
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
            double averageCliqueness = calculateAverageCliqueness(newInteractions, 
                                                                  ids);
            System.out.println("Average cliqueness: " + averageCliqueness);
        }
    }
    
    public double calculateAverageCliqueness(Set<String> interactions,
                                              Set<String> ids) {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(interactions);
        double totalCliqueness = 0.0;
        Double cliqueness;
        int totalCountable = 0;
        for (String id : ids) {
            cliqueness = bfs.calculateCliqueness(id, idToPartners);
            if (cliqueness == null)
                continue;
            //System.out.println("Cliquness of " + id + ": " + cliqueness);
            totalCliqueness += cliqueness;
            totalCountable ++;
        }
        return totalCliqueness / totalCountable;
    }
    
    @Test
    public void calculateAverageCliqueness() throws IOException {
        String intFileName = R3Constants.RESULT_DIR + "FI73_042108_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + ShortestPathAnalyzer.BIG_COMP_INT_FILE;
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        double averageCliqueness = calculateAverageCliqueness(interactions, ids);
        System.out.println("Average Cliqueness: " + averageCliqueness);
    }
    
}
