/*
 * Created on Jul 17, 2007
 *
 */
package org.reactome.r3.cluster;

import java.io.FileOutputStream;

import org.junit.Test;
import org.reactome.r3.util.R3Constants;

/**
 * For test running
 * @author guanming
 *
 */
public class SOMTest {
    
    @Test
    public void testOriginalSOM() throws Exception {
        OriginalSOM som = new OriginalSOM();
        // Want to get out the initial SOM
        String iniFileName = R3Constants.RESULT_DIR + "SOM_24_24_100_INI.txt";
        FileOutputStream fos = new FileOutputStream(iniFileName);
        som.output(fos);
        som.learn();
        som.label();
        String resultFileName = R3Constants.RESULT_DIR + "SOM_24_24_100.txt";
        fos = new FileOutputStream(resultFileName);
        som.output(fos);
    }
    
}
