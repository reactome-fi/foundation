/*
 * Created on Jun 1, 2015
 *
 */
package org.reactome.pagerank;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;


/**
 * This class is used to generate files needed in order to run HotNet2.
 * @author gwu
 *
 */
public class HotNet2FilesGenerator {
    private final String HOTNET2_DIR = "/Users/gwu/ProgramFiles/hotnet2-1.0.0/";
    private FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public HotNet2FilesGenerator() {
    }
    
    @Test
    public void parseResults() throws IOException {
        String dirName = HOTNET2_DIR + "hotnet_output/delta_0.00127124561251/";
        String srcFileName = dirName + "components.txt";
        fu.setInput(srcFileName);
        FileUtility outFu = new FileUtility();
        int count = 0;
        String line = null;
        int sizeCutoff = 4;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length >= sizeCutoff) {
                String outFileName = dirName + "component" + count + ".txt";
                outFu.setOutput(outFileName);
                for (String token : tokens)
                    outFu.printLine(token);
                outFu.close();
            }
            count ++;
        }
        fu.close();
    }
    
    /**
     * Generate a gene score file.
     * @throws IOException
     */
    @Test
    public void generateGeneScoreFile() throws IOException {
        String srcFileName = "results/ICGC_PanCancer/ICGCPanCancerScores_NoMutation_052815.txt";
        String outFileName = "results/ICGC_PanCancer/ICGCPanCancerScores_NoMutation_052815_Gene_Score.txt";
        fu.setInput(srcFileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fu.printLine(tokens[0] + "\t" + tokens[tokens.length - 1]);
        }
        fu.close();
    }
    
    /**
     * Generate random network by keeping the original degrees of genes in the FI network.
     * The Python code is too slow to be used for the ReactomeFI network.
     * @throws IOException
     */
    @Test
    public void permuteNetworks() throws IOException {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        // Generate a map from genes to indices
        Map<String, Integer> geneToIndex = new HashMap<String, Integer>();
        for (int i = 0; i < geneList.size(); i++)
            geneToIndex.put(geneList.get(i), i + 1);
        Map<String, Set<String>> geneToParnters = InteractionUtilities.generateProteinToPartners(fis);
        // Convert the above into map from gene index to gene's degree
        Map<Integer, Integer> indexToDegree = new HashMap<Integer, Integer>();
        for (String gene : geneToParnters.keySet()) {
            Integer index = geneToIndex.get(gene);
            Set<String> partners = geneToParnters.get(gene);
            indexToDegree.put(index, partners.size());
        }
        int numberOfPerm = 100;
        RandomData randomizer = new RandomDataImpl();
        generateRandomNetwork(geneList, 
                              geneToIndex,
                              indexToDegree,
                              numberOfPerm, 
                              randomizer);
    }

    private void generateRandomNetwork(List<String> geneList,
                                       Map<String, Integer> geneToIndex,
                                       Map<Integer, Integer> indexToDegree,
                                       int numberOfPerm, RandomData randomizer) throws IOException {
        // Directory to hold the randomized network
        String dirName = "/Users/gwu/ProgramFiles/hotnet2-1.0.0/influence_matrices/reactomefi/permuted/";
        for (int i = 0; i < numberOfPerm; i++) {
            String subDir = dirName + (i + 1) + "/";
            new File(subDir).mkdir();
            // Output gene index
            String geneIndexFile = subDir + "reactomefi_index_genes";
            fu.setOutput(geneIndexFile);
            for (int j = 0; j < geneList.size(); j++)
                fu.printLine((j + 1) + "\t" + geneList.get(j));
            fu.close();
            long time1 = System.currentTimeMillis();
            // Output edge list
            String edgeListFile = subDir + "reactomefi_edge_list";
            fu.setOutput(edgeListFile);
            for (Integer index : indexToDegree.keySet()) {
                Integer degree = indexToDegree.get(index);
                Set<Integer> randomIndices = MathUtilities.randomSampling(geneToIndex.values(),
                                                                          degree, 
                                                                          randomizer);
                randomIndices.remove(index); // Avoid itself
                for (Integer randomIndex : randomIndices) {
                    fu.printLine(index + "\t" + randomIndex + "\t1"); // 1 is assumed as edge weight
                }
            }
            fu.close();
            calculateInfluenceMatrix(subDir);
            long time2 = System.currentTimeMillis();
            System.out.println("Time for " + i + ": " + (time2 - time1));
        }
    }
    
    /**
     * Use ProcessBuilder to call a Python script to calculate the influence matrix.
     */
    private void calculateInfluenceMatrix(String dirName) throws IOException {
        String workingDirName = "/Users/gwu/ProgramFiles/hotnet2-1.0.0/";
        StringBuilder parameters = new StringBuilder();
        parameters.append("-e ").append(dirName).append("reactomefi_edge_list ");
        parameters.append("-i ").append(dirName).append("reactomefi_index_genes ");
        parameters.append("-p reactomefi ");
        parameters.append("-a 0.5 ");
        parameters.append("-o ").append(dirName);
//        System.out.println(parameters.toString());
        ProcessBuilder builder = new ProcessBuilder("bash", // have to go through bash. Otherwise, the process cannot be run.
                                                    "createPPRMat.sh",
                                                    parameters.toString());
        builder.directory(new File(workingDirName));
        builder.redirectErrorStream(true);
        Map<String, String> env = builder.environment();
        env.put("PYTHONPATH", ".:/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages");
        Process process = builder.start();
        InputStream is = process.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line = null;
        while ((line = br.readLine()) != null)
            System.out.println(line);
        br.close();
        isr.close();
        is.close();
    }
    
    @Test
    public void createGenesAndEdgesFiles() throws IOException {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        String output = HOTNET2_DIR + "influence_matrices/reactomefi/reactomefi_index_genes";
        fu.setOutput(output);
        Map<String, Integer> geneToIndex = new HashMap<String, Integer>();
        for (int i = 0; i < geneList.size(); i++) {
            fu.printLine((i + 1) + "\t" + geneList.get(i));
            geneToIndex.put(geneList.get(i), i + 1);
        }
        fu.close();
        output = HOTNET2_DIR + "influence_matrices/reactomefi/reactomefi_edge_list";
        fu.setOutput(output);
        for (String fi : fis) {
            String[] tokens = fi.split("\t");
            Integer index1 = geneToIndex.get(tokens[0]);
            Integer index2 = geneToIndex.get(tokens[1]);
            fu.printLine(index1 + "\t" + index2 + "\t1"); // 1 is added in the example file: assuming it is the weight.
        }
        fu.close();
    }
   
}
