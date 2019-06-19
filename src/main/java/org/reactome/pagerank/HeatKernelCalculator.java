/*
 * Created on Feb 14, 2013
 *
 */
package org.reactome.pagerank;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.ProcessRunner;
import org.reactome.r3.util.R3Constants;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * This class is used to calculate heat kernel based on this paper by Fan Chung:
 * A local graph partitioning algorithm using heat kernel pagerank.
 * @author gwu
 *
 */
public class HeatKernelCalculator {
    private final String R_MATRIX_SCRIPT = R3Constants.R_SRC_DIR + "MatrixExponential.R";
    
    public HeatKernelCalculator() {
    }
    
    /**
     * Use this method to calculate heat kernel
     * @param nodes
     * @param nodeToNeighbor
     * @param temp temperature
     * @return
     */
    public void calculateHeatKernel(List<String> nodes,
                                    Map<String, Set<String>> nodeToNeighbor,
                                    double temp) throws IOException {
        // Create the transition probability matrix W
        // Actually there is no need to create w. What we need is matrix L = I - W
        DoubleMatrix2D l = new DenseDoubleMatrix2D(nodes.size(), nodes.size());
        for (int i = 0; i < nodes.size(); i++) {
            String gene = nodes.get(i);
            Set<String> neighbor = nodeToNeighbor.get(gene);
            for (int j = 0; j < nodes.size(); j++) {
                String gene1 = nodes.get(j);
                if (neighbor.contains(gene1)) {
                    l.set(i, j, -1.0d / neighbor.size());
                }
                else if (i == j)
                    l.set(i, j, 1.0d);
            }
        }
        scalarMultiple(l, -temp);
        String inMatrixFileName = "tmp/I_matrix.txt";
        outputMatrix(l, inMatrixFileName);
        String outMatrixFileName = "tmp/Kernel.txt";
        outMatrixFileName = "tmp/HeatKernel_time_01_021913.txt";
        runMatrixExp(inMatrixFileName,
                     outMatrixFileName);
    }
    
    /**
     * This is a HotNet way to calculate heat kernel. For details, see
     * DISCOVERY OF MUTATED SUBNETWORKS ASSOCIATED WITH CLINICAL DATA IN CANCER.
     * @param nodes
     * @param nodeToNeighbor
     * @param temp
     * @throws IOException
     */
    public void calculateHotNetHeatKernel(List<String> nodes,
                                          Map<String, Set<String>> nodeToNeighbor,
                                          double temp) throws IOException {
        Collections.sort(nodes); // Sort all nodes so that we can retrieve them later on.
        // Create the transition probability matrix W
        // Actually there is no need to create w. What we need is matrix L = I - W
        DoubleMatrix2D l = new DenseDoubleMatrix2D(nodes.size(), nodes.size());
        for (int i = 0; i < nodes.size(); i++) {
            String gene = nodes.get(i);
            Set<String> neighbor = nodeToNeighbor.get(gene);
            for (int j = 0; j < nodes.size(); j++) {
                String gene1 = nodes.get(j);
                if (neighbor.contains(gene1)) {
                    l.set(i, j, -1.0d);
                }
                else if (i == j)
                    l.set(i, j, neighbor.size());
            }
        }
        scalarMultiple(l, -temp);
        String inMatrixFileName = "tmp/HotNet_L_matrix.txt";
        outputMatrix(l, inMatrixFileName);
//        String outMatrixFileName = "tmp/Kernel.txt";
//        outMatrixFileName = "tmp/HeatKernel_HotNet_time_01_021913.txt";
//        runMatrixExp(inMatrixFileName,
//                     outMatrixFileName);
    }
    
    private void scalarMultiple(DoubleMatrix2D matrix, final double scalar) {
        matrix.assign(new DoubleFunction() {
            @Override
            public double apply(double argument) {
                return argument * scalar;
            }
        });
    }
    
//    private void matrixAdd(DoubleMatrix2D matrix1, DoubleMatrix2D matrix2) {
//        matrix1.assign(matrix2, new DoubleDoubleFunction() {
//            
//            @Override
//            public double apply(double x, double y) {
//                return x + y;
//            }
//        });
//    }
    
    private void outputMatrix(DoubleMatrix2D matrix,
                              String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        // Use an double array for output is much faster than
        // matrix.getQuick(int, int)
        double[][] values = matrix.toArray();
        for (int i = 0; i < matrix.rows(); i++) {
            for (int j = 0; j < matrix.columns(); j++) {
                double value = values[i][j];
//                builder.append(String.format("%.5f", value));
                builder.append(value);
                builder.append("\t");
            }
            // Remove the last tab delimit
            builder.deleteCharAt(builder.length() - 1);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private void runMatrixExp(String inFileName,
                              String outFileName) throws IOException {
        ProcessRunner runner = new ProcessRunner();
        String[] parameters = new String[]{
                R_MATRIX_SCRIPT,
                inFileName,
                outFileName
        };
        String[] outputs = runner.runRScript(parameters);
        for (String output : outputs) {
            System.out.println(output);
        }
    }
    
    @Test
    public void testCalculateHeatKernel() throws IOException {
        FileUtility fu = new FileUtility();
//        String fiBigCompFile = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
        String fiBigCompFile = "results/v3/FIsInGene_041709_BigComp.txt";
        Set<String> fis = fu.loadInteractions(fiBigCompFile);
        Map<String, Set<String>> nodeToNeighbor = InteractionUtilities.generateProteinToPartners(fis);
        List<String> genes = new ArrayList<String>(nodeToNeighbor.keySet());
        double temp = 0.1d;
        long time1 = System.currentTimeMillis();
//        calculateHeatKernel(genes, nodeToNeighbor, temp);
        calculateHotNetHeatKernel(genes, nodeToNeighbor, temp);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
        // Output the matrix
//        outputMatrix(kernel, "tmp/HeatKernel_time_01_021913.txt");
    }
    
    @Test
    public void testRMatrixExponential() throws IOException {
        String inFileName = "tmp/example.txt";
        String outFileName = "tmp/example_exp.txt";
        ProcessRunner runner = new ProcessRunner();
        String[] parameters = new String[]{
                R_MATRIX_SCRIPT,
                inFileName,
                outFileName
        };
        runner.runRScript(parameters);
    }
    
    /**
     * Use this method to load a pre-generated heat kernel matrix, which is
     * generated by using R, and has the first row as its header.
     */
    public DoubleMatrix2D loadHeatKernelMatrix(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // The output from R has a header, which cannot be stripped out.
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        int size = tokens.length;
        int i = 0;
        int j = 0;
//        long time1 = System.currentTimeMillis();
        DoubleMatrix2D matrix = new DenseDoubleMatrix2D(size, size);
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (size != tokens.length)
                throw new IllegalStateException("Wrong tokens in a line: it should have " + size + " but has " + tokens.length);
            j = 0;
            for (String token : tokens) {
                matrix.set(i, j, Double.parseDouble(token)); // Using parseDouble is a little faster than new Double().
                j++;
            }
            i ++;
        }
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for parsing: " + (time2 - time1));
        fu.close();
        return matrix;
    }
    
    /**
     * Use this method to load a pre-generated, serialized heat kernel matrix. Loading a serialized heat kernel matrix
     * is much faster than loading a text file and parse it into a matrix.
     * @param fileName
     * @return
     * @throws IOException
     * @throws ClassNotFoundException
     */
    public DoubleMatrix2D loadSerializedHeatKernelMatrix(String fileName) throws IOException, ClassNotFoundException {
        FileInputStream fis = new FileInputStream(fileName);
        ObjectInputStream ois = new ObjectInputStream(fis);
        DoubleMatrix2D matrix = (DoubleMatrix2D) ois.readObject();
        ois.close();
        fis.close();
        return matrix;
    }
    
    @Test
    public void testLoadHeatKernel() throws Exception {
//        String heatKernelFileName = "results/pagerank/HeatKernel_HotNet_time_01_2012_021913.txt";
        String heatKernelFileName = "results/pagerank/HeatKernel_HotNet_time_01_2009_041713.txt";
        long time1 = System.currentTimeMillis();
        DoubleMatrix2D matrix = loadHeatKernelMatrix(heatKernelFileName);
        System.out.println("Matrix: " + matrix.rows() + ", " + matrix.columns());
        long time2 = System.currentTimeMillis();
        System.out.println("Total time for loading a matrix from a file: " + (time2 - time1));
        // Want to use object serialization
//        String fileName = "results/pagerank/HeatKernel_HotNet_time_01_2012_021913.ser";
        String fileName = "results/pagerank/HeatKernel_HotNet_time_01_2009_041713.ser";
        FileOutputStream fos = new FileOutputStream(fileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(matrix);
        oos.flush();
        oos.close();
        fos.close();
        long time3 = System.currentTimeMillis();
        System.out.println("Total time for saving matrix object: " + (time3 - time2));
        FileInputStream fis = new FileInputStream(fileName);
        ObjectInputStream ois = new ObjectInputStream(fis);
        DoubleMatrix2D loadedMatrix = (DoubleMatrix2D) ois.readObject();
        ois.close();
        fis.close();
        long time4 = System.currentTimeMillis();
        System.out.println("Total time for loading matrix object: " + (time4 - time3));
        // Do a random test
        Random random = new Random();
        int size = matrix.rows(); // Should have same columns and rows.
        for (int c = 0; c < 10; c ++) {
            int i = random.nextInt(size);
            int j = random.nextInt(size);
            System.out.printf("(%d, %d): %f\n", i, j, matrix.get(i, j));
            System.out.printf("(%d, %d): %f\n", j, i, matrix.get(j, i));
            System.out.printf("(%d, %d): %f\n", i, j, loadedMatrix.get(i, j));
            System.out.printf("(%d, %d): %f\n", j, i, loadedMatrix.get(j, i));
            System.out.println();
        }
    }
    
    @Test
    public void checkHeatKernelMatrix() throws IOException {
        long time1 = System.currentTimeMillis();
        String fileName = "results/pagerank/HeatKernel_HotNet_time_01_021913.txt";
        DoubleMatrix2D matrix = loadHeatKernelMatrix(fileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
        // Do a random test
        Random random = new Random();
        int size = matrix.rows(); // Should have same columns and rows.
        for (int c = 0; c < 10; c ++) {
            int i = random.nextInt(size);
            int j = random.nextInt(size);
            System.out.printf("(%d, %d): %f\n", i, j, matrix.get(i, j));
            System.out.printf("(%d, %d): %f\n", j, i, matrix.get(j, i));
            System.out.println();
        }
    }
    
}
