/*
 * This class is based on the paper by Sengoku, H. and Yoshihara, I titled by
 * A Fast TSP Solver Using GA on Java.
 * Created on Apr 3, 2007
 *
 */
package org.reactome.r3.graph.ga;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.random.JDKRandomGenerator;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;

public class FastTSPSolverViaGA {
    
    // Total size of population
    private final int POPULATION_SIZE = 100;
    // Percent of population should be selected away
    private final int SELECTION_SIZE = (int) (POPULATION_SIZE * 0.3);
    // Percent of population should go mutating
    private final int MUTATION_SIZE = (int) (POPULATION_SIZE * 0.2);
    // Used to choose neighbor species for selecting away
    private final int DIFF_THRESHOLD = 0;
    // The whole population
    private List<Chromosome> chromosomes;
    private Map<String, Integer> idPairToDist;
    
    public FastTSPSolverViaGA(List<String> ids,
                              Map<String, Integer> idPairToDist) {
        List<Gene> genes = createGenesFromIds(ids);
        init(genes, idPairToDist);
        this.idPairToDist = idPairToDist;
    }
    
    private void init(List<Gene> genes,
                      Map<String, Integer> idPairToDist) {
        chromosomes = new ArrayList<Chromosome>();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Chromosome c = new Chromosome(genes,
                                          idPairToDist);
            chromosomes.add(c);
        }
    }
    
    public Chromosome calculateOptimal() {
        double thisCost = 500.0;
        double oldCost = 0.0;
        int countSame = 0;
        int generation = 0;
        Chromosome elite = null;
        while (countSame < 50) {
            generation++;
            // Three evolve operations
            long time1 = System.currentTimeMillis();
            select();
            long time2 = System.currentTimeMillis();
            System.out.println("Select: " + (time2 - time1));
            multiplicate();
            long time3 = System.currentTimeMillis();
            System.out.println("Multiplicate: " + (time3 - time2));
            mutate();
            long time4 = System.currentTimeMillis();
            System.out.println("Mutate: " + (time4 - time3));
            GAUtilities.sortChromosomes(chromosomes);
            elite = chromosomes.get(0);
            thisCost = elite.getCost();
            System.out.println("Generation " + generation + ": " + thisCost);
            if (thisCost == oldCost) {
                countSame ++;
            }
            else {
                oldCost = thisCost;
                countSame = 0;
            }
        }
        return chromosomes.get(0);
    }
    
    private List<Gene> createGenesFromIds(List<String> ids) {
        List<Gene> genes = new ArrayList<Gene>(ids.size());
        for (String id : ids) {
            Gene g = new Gene(id);
            genes.add(g);
        }
        return genes;
    }
    
    /**
     * Eliminate chromosomes based on similarities and the cost functions.
     */
    private void select() {
        // Sort chromosomes first
        GAUtilities.sortChromosomes(chromosomes);
        // track the total eliminate number
        List<Chromosome> eliminate = new ArrayList<Chromosome>();
        // Using DIFF_THRESHOLD = 0 means don't take any similarities
        for (int i = 0; i < chromosomes.size() - 1; i++) {
            Chromosome c1 = chromosomes.get(i);
            for (int j = i + 1; j < chromosomes.size(); j++) {
                Chromosome c2 = chromosomes.get(j);
                int diff = c2.getCost() - c1.getCost();
                if (diff < DIFF_THRESHOLD) {
                    eliminate.add(c2);
                }
            }
        }
        if (eliminate.size() <= SELECTION_SIZE) {
            chromosomes.removeAll(eliminate); 
            int more = SELECTION_SIZE - eliminate.size();
            for (int i = 0; i < more; i++)
                chromosomes.remove(chromosomes.size() - 1);
        }
        else {
            int extra = eliminate.size() - SELECTION_SIZE;
            // Don't do this: an expensive step if eliminate is big
            //for (int i = 0; i < extra; i++)
            //    eliminate.remove(0);
            for (int i = eliminate.size() - 1; i > eliminate.size() - 1 - SELECTION_SIZE; i--) 
                chromosomes.remove(eliminate.get(i));
         }
    }
    
    private void multiplicate() {
        // Need to choose 2 * SELECTION_SIZE chromosomes for crossover
        // Use a no-seed RandomGenerator to avoid radom sequences are the same
        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
        Object[] crossoverPairs = randomizer.nextSample(chromosomes, 2 * SELECTION_SIZE);
        for (int i = 0; i < crossoverPairs.length - 1; i++) {
            Chromosome mother = (Chromosome) crossoverPairs[i];
            Chromosome father = (Chromosome) crossoverPairs[i + 1];
            Chromosome child = crossover(mother, father);
            chromosomes.add(child);
        }
    }
    
    private Chromosome crossover(Chromosome mother,
                                 Chromosome father) {
        Chromosome child = new Chromosome(idPairToDist);
        // Need to calculate geneList for this child
        List<Gene> listA = mother.getGenes();
        List<Gene> listB = father.getGenes();
        int totalGenes = listA.size();
        boolean fA = true;
        boolean fB = true;
        // Get a random
        Gene gene = getRandomGene(listA);
        int indexA = listA.indexOf(gene);
        int indexB = listB.indexOf(gene);
        List<Gene> childList = new ArrayList<Gene>();
        childList.add(gene);
        do {
            indexA --;
            if (indexA < 0)
                indexA += totalGenes; // Go to the other end
            indexB ++;
            if (indexB > totalGenes - 1)
                indexB -= totalGenes; // Start from 0
            if (fA) {
                Gene g = listA.get(indexA);
                if (childList.contains(g))
                    fA = false;
                else
                    childList.add(0, g); // Add to the head
            }
            if (fB) {
                Gene g = listB.get(indexB);
                if (childList.contains(g)) 
                    fB = false;
                else
                    childList.add(g); // Add to the bottom
            }
        }
        while (fA || fB);
        if (childList.size() < listA.size()) {
            List<Gene> copy = new ArrayList<Gene>(listA);
            copy.removeAll(childList);
            randomList(copy);
            for (Gene g : copy)
                childList.add(g);
        }
        child.setGenes(childList);
        return child;
    }
    
    private void randomList(List<Gene> genes) {
        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
        int size = genes.size();
        int[] order = randomizer.nextPermutation(size, size);
        List<Gene> copy = new ArrayList<Gene>(genes);
        for (int i = 0; i < size; i++)
            genes.set(i, copy.get(order[i]));
    }
    
    private Gene getRandomGene(List<Gene> genes) {
        int index = (int) (Math.random() * genes.size());
        return genes.get(index);
    }
    
    private void mutate() {
        GAUtilities.sortChromosomes(chromosomes);
        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
        int[] permutate = randomizer.nextPermutation(POPULATION_SIZE, MUTATION_SIZE);
        boolean isElitedChosen = false;
        for (int i : permutate) {
            if (i == 0)
                isElitedChosen = true;
            Chromosome c = chromosomes.get(i);
            mutate(c);
        }
        // Make sure elite is always chosen
        if (!isElitedChosen) {
            Chromosome c = chromosomes.get(0);
            mutate(c);
        }
    }
    
    private void mutate(Chromosome c) {
        List<Gene> geneList = c.getGenes();
        int size = geneList.size();
        boolean isSwapped = true;
        while (isSwapped) {
            isSwapped = false;
            // Whenever there is a swap the whole whole should be recompute
            for (int i = 0; i < size - 3; i ++) {
                // These two genes are for the first edge
                Gene gene1 = geneList.get(i);
                Gene gene2 = geneList.get(i + 1);
                int dist1 = GAUtilities.distance(gene1, gene2, idPairToDist);
                for (int j = i + 2; j < size; j++) {
                    Gene gene3 = geneList.get(j);
                    Gene gene4 = null;
                    if (j == size - 1)
                        gene4 = geneList.get(0);
                    else
                        gene4 = geneList.get(j + 1);
                    int dist2 = GAUtilities.distance(gene3, gene4, idPairToDist);
                    // Check possible distance
                    int dist3 = GAUtilities.distance(gene1, gene3, idPairToDist);
                    int dist4 = GAUtilities.distance(gene2, gene4, idPairToDist);
                    if (dist1 + dist2 > dist3 + dist4) {
                        // Need to swap gene2 and gene3
                        geneList.set(i + 1, gene3);
                        geneList.set(j, gene2);
                        gene2 = gene3;
                        isSwapped = true;
                        dist1 = GAUtilities.distance(gene1, gene2, idPairToDist);
                    }
                }
            }
        }
        c.calculateCost();
    }
    
    public static void main(String[] args) {
        // A test for nodes in a circle
        List<String> ids = new ArrayList<String>();
        for (int i = 0; i < 100; i++) {
            ids.add(i + "");
        }
        // Generate map
        Map<String, Integer> distMap = new HashMap<String, Integer>();
        for (int i = 0; i < 99; i++) {
            String id1 = i + "";
            for (int j = i + 1; j < 100; j++) {
                String id2 = j + "";
                int dist = j - i;
                if (i == 0 && j == 99)
                    dist = 1;
                int compare = id1.compareTo(id2);
                if (compare < 0)
                    distMap.put(id1 + "->" + id2, dist);
                else
                    distMap.put(id2 + "->" + id1, dist);
            }
        }
        FastTSPSolverViaGA tspSolver = new FastTSPSolverViaGA(ids,
                                                              distMap);
        Chromosome c = tspSolver.calculateOptimal();
        System.out.println("Solution: " + c.getCost());
        for (Gene gene : c.getGenes())
            System.out.println(gene);
    }
}
