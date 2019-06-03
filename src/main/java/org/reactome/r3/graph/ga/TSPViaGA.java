/*
 * Created on Apr 2, 2007
 *
 */
package org.reactome.r3.graph.ga;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * This class is used to calculate TSP by using genetic algorithm.
 * @author guanming
 *
 */
public class TSPViaGA {
    /**
     * How many chromosomes to use.
     */
    private static final int POPULATION_SIZE = 1000;

    /**
     * What percent of new-borns to mutate.
     */
    private static final double MUTATION_PERCENT = 0.10;

    /**
     * The part of the population eligable for mateing.
     */
    private int matingPopulationSize = POPULATION_SIZE / 2;

    /**
     * The part of the population favored for mating.
     */
    private int favoredPopulationSize = matingPopulationSize / 2;
    
    /**
     * How much genetic material to take during a mating.
     */
    protected int cutLength;
    
    private List<Chromosome> chromosomes;
        
    public TSPViaGA(List<String> ids,
                    Map<String, Integer> idPairToDist) {
        List<Gene> genes = createGenesFromIds(ids);
        init(genes, idPairToDist);
    }
    
    private void init(List<Gene> genes,
                      Map<String, Integer> idPairToDist) {
        cutLength = genes.size() / 5;
        chromosomes = new ArrayList<Chromosome>();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Chromosome c = new Chromosome(genes,
                                          idPairToDist);
            chromosomes.add(c);
        }
    }
    
    private List<Gene> createGenesFromIds(List<String> ids) {
        List<Gene> genes = new ArrayList<Gene>(ids.size());
        for (String id : ids) {
            Gene g = new Gene(id);
            genes.add(g);
        }
        return genes;
    }
    
    public Chromosome calculateOptimal() {
        double thisCost = 500.0;
        double oldCost = 0.0;
        double dcost = 500.0;
        int countSame = 0;
        int generation = 0;
        
        while(countSame < 100) {

          generation++;
          int ioffset = matingPopulationSize;
          int mutated = 0;

          // Mate the chromosomes in the favoured population
          // with all in the mating population
          for ( int i = 0; i < favoredPopulationSize; i++ ) {
            Chromosome cmother = chromosomes.get(i);
            // Select partner from the mating population
            int father = (int) (Math.random() * matingPopulationSize);
            Chromosome cfather = chromosomes.get(father);

            mate(cmother,
                 cfather,
                 chromosomes.get(ioffset),
                 chromosomes.get(ioffset+1));
            ioffset += 2;
          }

          // Now sort the new mating population
          sortChromosomes(chromosomes);

          double cost = chromosomes.get(0).getCost();
          System.out.println("Optimal cost: " + cost);
          dcost = Math.abs(cost-thisCost);
          thisCost = cost;
          double mutationRate = 100.0 * (double) mutated / matingPopulationSize;

          if ( (int)thisCost == (int)oldCost ) {
            countSame++;
          } 
          else {
            countSame = 0;
            oldCost = thisCost;
          }
        }
        return chromosomes.get(0);
    }
    
    /**
     * Sort the chromosomes by their cost.
     *
     * @param chromosomes An array of chromosomes to sort.
     * @param num How much of the chromosome list to sort.
     */
    private void sortChromosomes(List<Chromosome> chromosomes) {
        Collections.sort(chromosomes, new Comparator<Chromosome>() {
            public int compare(Chromosome c1, Chromosome c2) {
                double d1 = c1.getCost();
                double d2 = c2.getCost();
                if (d1 < d2)
                    return -1;
                if (d1 > d2)
                    return 1;
                return 0;
            }
        });
    }
    
    /**
     * Assuming this chromosome is the "mother" mate with
     * the passed in "father".
     *
     * @param father The father.
     * @param offspring1 Returns the first offspring
     * @param offspring2 Returns the second offspring.
     * @return The amount of mutation that was applied.
     */
    public void mate(Chromosome mother,
                     Chromosome father, 
                     Chromosome offspring1, 
                     Chromosome offspring2) {
        int size = mother.getGenes().size();
        int cutpoint1 = (int) (Math.random() * (size - cutLength));
        int cutpoint2 = cutpoint1 + cutLength;
        
        boolean taken1 [] = new boolean[size];
        boolean taken2 [] = new boolean[size];
        // For mother (this)
        List<Gene> off1 = new ArrayList<Gene>(size);
//      For father (father from passed)
        List<Gene> off2 = new ArrayList<Gene>(size);
        for (int i = 0; i < size; i++) {
            off1.add(null);
            off2.add(null);
        }
        
        // Cross
        for ( int i=0; i<size; i++) {
            if (i >= cutpoint1 && i < cutpoint2) {
                off1.set(i, father.getGene(i));
                off2.set(i, mother.getGene(i)); // This is mother
            }
        }
        
        // Popup the top part
        for ( int i=0; i<cutpoint1; i++ ) {
            // off1
            for (int j = 0; j < size; j++) {
                Gene g = mother.getGene(j);
                if (!off1.contains(g)) {
                    off1.set(i, g);
                    break;
                }
            }
            // off2
            for (int j = 0; j < size; j++) {
                Gene g = father.getGene(j);
                if (!off2.contains(g)) {
                    off2.set(i, g);
                    break;
                }
            }
        }
        // Another part
        for ( int i = size - 1; i >= cutpoint2; i-- ) {
            for (int j = size - 1; j >= 0; j--) {
                Gene g = mother.getGene(j);
                if (!off1.contains(g)) {
                    off1.set(i, g);
                    break;
                }
            }
            for (int j = size - 1; j >= 0; j--) {
                Gene g = father.getGene(j);
                if (!off2.contains(g)) {
                    off2.set(i, g);
                    break;
                }
            }
        }
        
        // Mutation
        if (Math.random() < MUTATION_PERCENT) {
            int iswap1 = (int) (Math.random() * size);
            int iswap2 = (int) (Math.random() * size);
            Gene g1 = off1.get(iswap1);
            Gene g2 = off1.get(iswap2);
            off1.set(iswap1, g2);
            off1.set(iswap2, g1);
        }
        if (Math.random() < MUTATION_PERCENT) {
            int iswap1 = (int) (Math.random() * size);
            int iswap2 = (int) (Math.random() * size);
            Gene g1 = off2.get(iswap1);
            Gene g2 = off2.get(iswap2);
            off2.set(iswap1, g2);
            off2.set(iswap2, g1);
        }
        offspring1.setGenes(off1);
        offspring2.setGenes(off2);
    }
    
}
