/*
 * Created on Apr 2, 2007
 *
 */
package org.reactome.r3.graph.ga;

public class Gene {
    private String id;
    
    public Gene() {
    }
    
    public Gene(String id) {
        this();
        this.id = id;
    }
    
    public void setId(String id) {
        this.id = id;
    }
    
    public String getId() {
        return this.id;
    }
    
    public String toString() {
        return id;
    }
}
