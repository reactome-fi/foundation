/*
 * Created on Jan 10, 2008
 *
 */
package org.reactome.r3.graph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.reactome.r3.util.InteractionUtilities;

/**
 * This is a specific implementation of FloydWarshall used to calculate shortest path
 * among all pairs of a provided FI network.
 * Note: this is not a generic implementation of FloydWarshall algorithm.
 * @author guanming
 *
 */
public class FloydWarshall {
    // The following two dimensional integer arrays are used
    // to hold distance information.
    private int[][] distMinusOne;
    private int[][] dist;
    
    public FloydWarshall() {
    }
    
    /**
     * Protein ids in the passed interactions should be sorted: the first id should be
     * less than the second.
     * @param interactions
     */
    private void setInteractions(Set<String> interactions) {
        // Get all protein ids first
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        List<String> idList = new ArrayList<String>(ids);
        Collections.sort(idList);
        // Generate adjMatrix
        int size = idList.size();
        distMinusOne = new int[size][];
        String pair = null;
        long time1 = System.currentTimeMillis();
        for (int i = 0; i < size; i++)
            distMinusOne[i] = new int[size];
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                pair = idList.get(i) + " " + idList.get(j);
                if (interactions.contains(pair))
                    distMinusOne[i][j] = 1;
                else
                    distMinusOne[i][j] = Integer.MAX_VALUE;
            }
        }
        // Do a copy
        for (int i = 0; i < size - 1; i++)
            for (int j = i + 1; i < size; i++)
                distMinusOne[j][i] = distMinusOne[i][j];
        long time2 = System.currentTimeMillis();
        System.out.println("Time for setting: " + (time2 - time1));
    }
    
    public void calculateAllPairs(Set<String> interactions) {
        // k = -1
        setInteractions(interactions);
        int size = distMinusOne.length; 
        // Iterate over k. At each step, the previously computed matrix
        // should be copied
        System.out.println("Total size: " + size);
        // Initialize the distance matrix
        dist = new int[size][];
        for (int i = 0; i < size; i++)
            dist[i] = new int[size];
        for (int k = 0; k < size; k++) {
            System.out.println("k: " + k);
            // comput dist[i][j]
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (i == j)
                        continue;
                    dist[i][j] = Math.min(distMinusOne[i][j], 
                                          distMinusOne[i][k] + distMinusOne[k][j]);
                }
            }
            // Need to copy dist to distMinusOne
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    distMinusOne[i][j] = dist[i][j];
        }
    }
    
    /**
     * Return the calculated results. The results are returned in
     * a two-dimension integer array. The original passed proteins
     * have been sorted alphabetically. So the client to this class
     * should figure out what is what by itself.
     * @return
     */
    public int[][] getAllPairDistances() {
        return dist;
    }
    
}
