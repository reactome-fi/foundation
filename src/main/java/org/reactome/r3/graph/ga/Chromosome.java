/*
 * Created on Apr 2, 2007
 * Modified from the above source.
 */
package org.reactome.r3.graph.ga;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.random.JDKRandomGenerator;
import org.apache.commons.math.random.RandomDataImpl;

public class Chromosome {
    
    protected List<Gene> geneList;
    private Map<String, Integer> distMap;
    
    private int cost;
    
    /**
     * The constructor
     *
     * @param cities The order that this chromosome would
     * visit the
     * cities. These cities can be thought of as the
     * genes of this chromosome.
     */
    public Chromosome(List<Gene> genes,
                      Map<String, Integer> distMap) {
        this.distMap = distMap;
        this.geneList = new ArrayList<Gene>();
        cost = 0;
        int size = genes.size();
        // Need to initialize
        for (int i = 0; i < size; i++)
            geneList.add(null);
        int[] random = new RandomDataImpl(new JDKRandomGenerator()).nextPermutation(size, size);
        for (int i = 0; i < size; i++)
            geneList.set(i, genes.get(random[i]));
        calculateCost();
        //System.out.println("Gene List: " + geneList);
    }
    
    public Chromosome(Map<String, Integer> distMap) {
        this.distMap = distMap;
        this.geneList = new ArrayList<Gene>();
    }
    
    /**
     * Calculate the cost of of the specified list of
     * cities.
     *
     * @param cities A list of cities.
     */
    public void calculateCost() {
        cost = 0;
        Gene gene1, gene2;
        for (int i = 0; i < geneList.size() - 1; i++) {
            gene1 = geneList.get(i);
            gene2 = geneList.get(i + 1);
            cost += GAUtilities.distance(gene1, gene2, distMap);
        }
        gene1 = geneList.get(geneList.size() - 1);
        gene2 = geneList.get(0);
        cost += GAUtilities.distance(gene1, gene2, distMap);
    }
    
    /**
     * Get the cost for this chromosome. This is the
     * amount of distance that must be traveled
     * for this chromosome's plan for moving through
     * the cities. The cost of a chromosome determines
     * its rank in terms of being able to mate and not
     * being killed off.
     */
    public int getCost() {
        return cost;
    }
    
    /**
     * Get the ith city in this chromosome. For example the
     * value 3 would return the fourth city this chromosome
     * would visit. Similarly 0 would return the first city
     * that this chromosome would visit.
     *
     * @param i The city you want.
     * @return The ith city.
     */
    public Gene getGene(int i) {
        return geneList.get(i);
    }
    
    public List<Gene> getGenes() {
        return this.geneList;
    }
    
    /**
     * Set the order of cities that this chromosome
     * would visit.
     *
     * @param list A list of cities.
     */
    public void setGenes(List<Gene> genes) {
        geneList.clear();
        geneList.addAll(genes);
        calculateCost();
    }
    
    /**
     * Set the index'th city in the city list.
     *
     * @param index The city index to change
     * @param value The city number to place into the index.
     */
    public void getGene(int index, Gene gene) {
        geneList.set(index, gene);
    }
}
