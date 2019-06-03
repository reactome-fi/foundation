/*
 * Created on Apr 27, 2010
 *
 */
package org.reactome.r3.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.util.InteractionUtilities;

/**
 * This is a simple clique search. Should not be applied for generic clique search.
 * @author wgm
 *
 */
public class SimpleCliqueSearcher {
  
    public SimpleCliqueSearcher() {
    }
    
    /**
     * This method is used to check if the found clique is correct by checking neighbors FI partners.
     * @param clique
     * @param totalFIs
     * @return true for new ids have been added
     */
    public boolean validateClique(Collection<String> clique,
                                  Set<String> totalFIs) {
        Map<String, Set<String>> idToPartners = new BreadthFirstSearch().generateIdToPartnersMap(totalFIs);
        int originalSize = clique.size();
        int preSize = originalSize;
        while (true) {
            Set<String> neighbors = new HashSet<String>();
            for (String id : clique) {
                Set<String> partners = idToPartners.get(id);
                neighbors.addAll(partners);
            }
            neighbors.removeAll(clique);
            // Check if any neighbor can be a member of clique: in other words, any neighbor
            // can interact with all clique
            for (String n : neighbors) {
                Set<String> partners = idToPartners.get(n);
                if (partners.containsAll(clique)) {
                    //System.out.println(n + " should be in clique!");
                    clique.add(n);
                    break;
                }
            }
            if (preSize == clique.size())
                break;
            preSize = clique.size();
        }
        return ! (originalSize < clique.size());
    }
    
    /**
     * This simple method is used to find an embedded clique from the passed FIs.
     * @param fis
     * @return
     */
    public Collection<String> searchClique(Set<String> fis) {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> idList = new ArrayList<String>(ids);
        Collections.sort(idList);
        int count = 0;
        Set<String> currentFIs = new HashSet<String>(fis);
        while (true) {
            //System.out.println(count + ": " + idList.size());
            count ++;
            String toBeRemoved = null;
            Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(currentFIs);
            for (Iterator<String> it = idList.iterator(); it.hasNext();) {
                String id = it.next();
                Set<String> partners = idToPartners.get(id);
                if (partners == null)
                    it.remove();
            }
            for (int i = 0; i < idList.size() - 1; i++) {
                String id1 = idList.get(i);
                Set<String> partners1 = idToPartners.get(id1);
                for (int j = i + 1; j < idList.size(); j++) {
                    String id2 = idList.get(j);
                    if (!(currentFIs.contains(id1 + "\t" + id2))) {
                        Set<String> partners2 = idToPartners.get(id2);
                        if (partners1.size() < partners2.size()) 
                            toBeRemoved = id1;
                        else
                            toBeRemoved = id2;
                        break;
                    }
                }
                if (toBeRemoved != null)
                    break;
            }
            if (toBeRemoved == null)
                break;
            idList.remove(toBeRemoved);
            currentFIs = InteractionUtilities.getFIs(idList, fis);
            if (currentFIs.size() == 0)
                break;
        }
        return idList;
    }
    
}
