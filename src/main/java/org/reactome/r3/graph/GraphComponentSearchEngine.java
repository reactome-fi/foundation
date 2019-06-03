/*
 * Created on Jun 9, 2008
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.r3.graph.GraphComponent.ScoreCalculator;
import org.reactome.r3.util.R3Constants;

/**
 * This class implements a greedy-based search algorithms to find a graph component.
 * The implementation of this class is based on some methods in class 
 * csplugins.jActiveModules.GreedySearchThread.
 * @author wgm
 *
 */
public class GraphComponentSearchEngine {
    // Search results should be here
    private Map<String, GraphComponent> nodeToComponents;
    // predefined depth
    private int maxDepth;
    // predefined search length
    private int searchDepth;
    // Cached values to avoid a long passing parameters
    private Graph<String, DefaultEdge> graph;
    // Map from a vertex to its neighbors
    private Map<String, Set<String>> vertexToNeighbors = new HashMap<String, Set<String>>();
    private double bestScore;
    // Used to calculate score for component
    private ScoreCalculator scoreCalculator;
    
    public GraphComponentSearchEngine() {
    }
    
    /**
     * The maximum distance from the seed to its neighbor that can be used for search.
     * @param depth
     */
    public void setMaxDepth(int depth) {
        this.maxDepth = depth;
    }
    
    /**
     * The length that should be used for searching.
     * @param depth
     */
    public void setSearchDepth(int depth) {
        this.searchDepth = depth;
    }
    
    public ScoreCalculator getScoreCalculator() {
        if (scoreCalculator == null) {
            scoreCalculator = new MutualInformationScoreCalculator();
        }
        return scoreCalculator;
    }
    
    public void setScoreCalculator(ScoreCalculator calculator) {
        this.scoreCalculator = calculator;
    }
    
    /**
     * Test method.
     * @throws IOException
     */
    @Test
    public void testSearch() throws IOException {
        String fileName = R3Constants.INTERACTION_FILE_NAME;
        setMaxDepth(2);
        setSearchDepth(2);
        search(fileName);
    }
    
    /**
     * This method is used to start the algorithm from a passed interaction file.
     * @param intFileName
     * @throws IOException
     */
    public void search(String intFileName) throws IOException {        
        GraphAnalyzer analyzer = new GraphAnalyzer();
        Graph<String, DefaultEdge> graph = analyzer.createGraph(intFileName);
        search(graph);
    }
    
    public void search(Graph<String, DefaultEdge> graph) {
        if (nodeToComponents == null)
            nodeToComponents = new HashMap<String, GraphComponent>();
        else
            nodeToComponents.clear();
        this.graph = graph;
        Set<String> vertices = graph.vertexSet();
        int c = 0;
        for (String vertex : vertices) {
            searchForSeed(vertex);
            // Just a quick check
            GraphComponent comp = nodeToComponents.get(vertex);
            if (comp == null)
                continue; // Cannot start the search
            System.out.println(c + "\t" + vertex + "\t" + comp.getScore() + "\t" + comp.getAllNodes().size());
            c++;
        }
    }
    
    private void searchForSeed(String vertex) {
        // Get all linked nodes in a predefined distance
        Set<String> possibleNodes = new HashSet<String>();
        getMaxNeighbors(vertex, possibleNodes);
        GraphComponent component = new GraphComponent();
        component.setScoreCalculator(getScoreCalculator());
        component.addNode(vertex);
        bestScore = Double.NEGATIVE_INFINITY;
        // Start actual searching
        searchForSeedRecursively(vertex, 
                                 component,
                                 searchDepth,
                                 possibleNodes);
        searchForRemovable(component);
        // Do a sorting
        Set<String> nodes = component.getAllNodes();
        for (String node : nodes) {
            GraphComponent comp = nodeToComponents.get(node);
            if (comp == null ||
                comp.getScore() < component.getScore()) {
                nodeToComponents.put(node, component);
            }
        }
    }
    
    /**
     * This method is used to do a reverse searching by removing the leaf nodes in the
     * GraphComonent.
     * @param component
     */
    private void searchForRemovable(GraphComponent component) {
        Set<String> removable = component.getRemovableNodes();
        Set<String> checked = new HashSet<String>();
        boolean isChanged = false;
        while (removable.size() > 0) {
            isChanged = false;
            for (String v : removable) {
                checked.add(v);
                // Try it first
                String source = component.getSourceForRemovableNode(v);
                component.removeNode(v);
                if (component.getScore() > bestScore) {
                    bestScore = component.getScore();
                    isChanged = true;
                }
                else { // Should add the removed node back
                    component.addNode(v);
                    if (source != null)
                        component.addEdge(source, v);
                }
            }
            // Need to do another round of testing
            if (!isChanged)
                break;
            removable = component.getRemovableNodes();
            removable.removeAll(checked);
        }
    }
    
    /**
     * The main method for a recursive search for a specified node.
     * @param v
     * @param component
     * @param depth
     * @param possibleNodes
     * @return
     */
    private boolean searchForSeedRecursively(String v,
                                             GraphComponent component,
                                             int depth,
                                             Set<String> possibleNodes) {
        boolean improved = false;
        // Don't use >= which will absorb many genes don't have expression values
        if (component.getScore() > bestScore) {
            bestScore = component.getScore();
            depth = searchDepth; // Restart the search
            improved = true;
        }
        if (depth > 0) {
            boolean anyCallImproved = false;
            int dependentCount = 0;
            Set<String> neighbors = getNeighbors(v);
            for (String neighbor : neighbors) {
                if (possibleNodes.contains(neighbor) &&
                    !component.containsNode(neighbor)) {
                    component.addNode(neighbor);
                    component.addEdge(v, neighbor);
                    boolean thisCallImproved = searchForSeedRecursively(neighbor,
                                                                        component,
                                                                        depth - 1,
                                                                        possibleNodes);
                    if (!thisCallImproved) {
                        component.removeNode(neighbor);
                    } // end of if ()
                    else {
                        dependentCount += 1;
                        anyCallImproved = true;
                    } 
                } 
            }
            improved |= anyCallImproved;
        }
        return improved;
    }
    
    private void getMaxNeighbors(String vertex,
                                 Set<String> neighbors) {
        int step = 0;
        Set<String> current = new HashSet<String>();
        current.add(vertex);
        Set<String> next = new HashSet<String>();
        while (step <= maxDepth && current.size() > 0) {
            for (String v : current) {
                neighbors.add(v);
                Set<String> closeNeighbors = getNeighbors(v);
                next.addAll(closeNeighbors);
            }
            next.removeAll(neighbors);
            current.clear();
            current.addAll(next);
            next.clear();
            step ++;
        }
        neighbors.remove(vertex); // Remove itself
    }
    
    /**
     * A helper method to get neighbors for a String node.
     * @param v
     * @return
     */
    private Set<String> getNeighbors(String v) {
        Set<String> neighbors = vertexToNeighbors.get(v);
        if (neighbors != null)
            return neighbors;
        Set<DefaultEdge> edges = graph.edgesOf(v);
        // grep the next layer
        neighbors = new HashSet<String>();
        vertexToNeighbors.put(v, neighbors);
        for (DefaultEdge e : edges) {
            String s = graph.getEdgeSource(e);
            String t = graph.getEdgeTarget(e);
            if (s == v)
                neighbors.add(t);
            else
                neighbors.add(s);
        }
        return neighbors;
    }
    
    /**
     * Return the found components. GraphComponents in the returned list are sorted based scores
     * in the ascending order.
     * @return
     */
    public List<GraphComponent> getFoundComponents() {
        Set<GraphComponent> componentSet = new HashSet<GraphComponent>(nodeToComponents.values());
        postProcessComponents(componentSet);
        // Want to do a sorting based on scores
        List<GraphComponent> componentList = new ArrayList<GraphComponent>(componentSet);
        Collections.sort(componentList, new Comparator<GraphComponent>() {
            public int compare(GraphComponent comp1, 
                               GraphComponent comp2) {
                Double score1 = comp1.getScore();
                Double score2 = comp2.getScore();
                return score2.compareTo(score1); // In descending order
            }
        });
        return componentList;
    }
    
    private void postProcessComponents(Set<GraphComponent> comps) {
        // Filter all components that have no scores
        for (Iterator<GraphComponent> it = comps.iterator(); it.hasNext();) {
            GraphComponent comp = it.next();
            if (comp.getScore() == Double.NEGATIVE_INFINITY)
                it.remove();
        }
    }
    
}
