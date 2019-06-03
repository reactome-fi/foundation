/*
 * Created on Apr 4, 2007
 *
 */
package org.reactome.r3.graph.ga;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

public class GAUtilities {

    /**
     * Sort the chromosomes by their cost.
     *
     * @param chromosomes An array of chromosomes to sort.
     * @param num How much of the chromosome list to sort.
     */
    public static void sortChromosomes(List<Chromosome> chromosomes) {
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
    
    public static int distance(Gene gene1,
                                  Gene gene2,
                                  Map<String, Integer> distMap) {
        // create a key for search
        String id1 = gene1.getId();
        String id2 = gene2.getId();
        String key = null;
        int compare = id1.compareTo(id2);
        if (compare < 0)
            key = id1 + "->" + id2;
        else
            key = id2 + "->" + id1;
        return distMap.get(key);
    }
    
}
