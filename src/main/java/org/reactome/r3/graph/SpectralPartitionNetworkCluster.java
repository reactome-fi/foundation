/*
 * Created on Mar 22, 2010
 *
 */
package org.reactome.r3.graph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

/**
 * The following implementation is based on Newman, MEJ. PNAS 103(23): 8577-8582 (2006).
 * @author wgm
 *
 */
public class SpectralPartitionNetworkCluster {
    private final boolean DEBUG = false;
    private final int CLUSTER_SIZE_CUTOFF = 4;
    // Used to check if zero. This number may need to be changed in different platform. 
    // Under MacOS, this should be sufficient.
    private final double ZERO_CUTOFF = 1.0e-9;
    // Cache this map to avoid run converting of FIs multiple times
    private Map<String, Set<String>> idToPartners = null;
    
    public SpectralPartitionNetworkCluster() {
    }
    
    /**
     * Cluster a set of FIs based on a method described in Newman, MEJ. PNAS (2006) 103: 8577-8582.
     * There are several steps involved in this method:
     * 1). Construct a network based on the passed interactions.
     * 2). Check the network to see how many graph components are.
     * 3). Run clustering on each graph component. If a graph component has a size <= CLUSTER_SIZE_CUTOFF (pre-defined),
     * just return
     * 4). Generate a matrix as described in the original paper for a graph component.
     * 5). Do an EVD for the matrix
     * 6). Find the biggest positive eigenvalue and its corresponding eigenvector
     * 7). Split the network into two clusters based on the eigenvector from 6)
     * 8). Run fine tuning method by Kernighan-Lin algorithm (a kind of).
     * Note: a passed set of FIs may form several connected graph components. The algorithm
     * runs on each graph component to calculate total edges. However, the degrees for each
     * vertex has not changed. So idToPartnerss needs to be calculated once, and is cached in 
     * this implementation.
     * More note: the method used to check if the modularity can be increase after flipping a node from one cluster
     * to another cluster is based on the following: based on eq 1 in the original paper, in order to calculate
     * the difference of modularity after a node flip, we only need to focus on elements related to k, the index of
     * node has been flipped (or to be flipped). All other elements should vanish during calcullation since nothing
     * has been changed. Based on this, we can write DELTA_Q = -4 * Sk * [SIGMA(i in G, i != k)(Aki - Kk * Ki / 2m) * Si].
     * In sum (SIGMA), k should be excluded since k is changed. By caching this sum for each node, we can get an 
     * implementation in O[(m + n)n].
     * @param fis
     * @return a list of set of String. Using set so that it can be used by other classes in other packages (esp the cancer package).
     */
    public List<Set<String>> cluster(Set<String> fis) {
        idToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
        final List<Set<String>> clusters = new ArrayList<Set<String>>();
        List<Set<String>> graphComponents = JGraphTUtilities.calculateGraphComponents(fis);
        for (Set<String> component : graphComponents) {
            Set<String> fisInComponent = InteractionUtilities.getFIs(component, fis);
            List<String> idList = new ArrayList<String>(component);
            Collections.sort(idList);
            cluster(idList,
                    fisInComponent,
                    false,
                    clusters);
        }
        Collections.sort(clusters, new Comparator<Set<String>>() {
            public int compare(Set<String> cluster1, Set<String> cluster2) {
                return cluster2.size() - cluster1.size();
            }
        });
        return clusters;
    }
    
    public List<Set<String>> clusterGenes(Set<String> genes) throws IOException {
        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
        return cluster(fisInGenes);
    }
    
    private double calculateModularity(EigenvalueDecomposition evd,
                                       int[] clusterIndex) {
        double modularity = 0.0;
        DoubleMatrix1D eigenvalues = evd.getRealEigenvalues();
        DoubleMatrix2D eigenmatrix = evd.getV();
        for (int i = 0; i < eigenvalues.size(); i++) {
            double sub = 0.0;
            for (int j = 0; j < eigenmatrix.rows(); j++) {
                sub += eigenmatrix.get(j, i) * clusterIndex[j];
//                if (clusterIndex[j] == 1)
//                    sub += eigenmatrix.get(j, i);
//                else
//                    sub -= eigenmatrix.get(j, i);
            }
            //System.out.println("sub value: " + i + ": " + sub + ", eigenvalue " + eigenvalues.get(i));
            modularity += sub * sub * eigenvalues.get(i);
        }
        return modularity;
    }
    
    /**
     * Use double arrays instead of calling DoubleMatrix to increase a performance much, much better.
     * @param evdMatrix
     * @param eigenvalues
     * @param clusterIndex
     * @param flipIndex
     * @return
     */
    private double calculateModularity(double[][] evdMatrix,
                                       double[] eigenvalues,
                                       int[] clusterIndex,
                                       int flipIndex) {
        double modularity = 0.0;
        for (int i = 0; i < eigenvalues.length; i++) {
            double sub = 0.0;
            double[] eigenvector = evdMatrix[i];
            for (int j = 0; j < eigenvector.length; j++) {
                if (j == flipIndex) {
                    sub += eigenvector[j] * (-clusterIndex[j]);
                }
                else {
                    sub += eigenvector[j] * clusterIndex[j];
                    // Comparion is slower than the above multiplication
//                    if (clusterIndex[j] == 1)
//                        sub += eigenvector[j];
//                    else
//                        sub -= eigenvector[j];
                }
            }
            modularity += sub * sub * eigenvalues[i];
        }
        return modularity;
    }
    
    private void cluster(List<String> idsInComponent,
                         Set<String> fisInComponent,
                         boolean isForSubCluster,
                         List<Set<String>> clusters) {
        // If a graph component has size not greater than CLUSTER_SIZE_CUTOFF,
        // Treat it as a cluster
        if (idsInComponent.size() <= CLUSTER_SIZE_CUTOFF) {
            clusters.add(new HashSet<String>(idsInComponent));
            return; 
        }
        // Convert to a matrix
        DoubleMatrix2D matrix = convertFIsToMatrix(idsInComponent, 
                                                   fisInComponent,
                                                   isForSubCluster);
        EigenvalueDecomposition evd = new EigenvalueDecomposition(matrix);
        int index = searchForPositiveBiggestEigenvalue(evd);
        if (index < 0) {
            // This is done: no more split can be done
            clusters.add(new HashSet<String>(idsInComponent));
            return;
        }
        if (DEBUG)
            System.out.println("\nLargest eigenvalue index: " + index);
        double[] eigenvector = getEigenvector(evd, index);
        int[] clusterIndex = generateClusterVector(eigenvector);
        double modularity = fineTuneOnKernighanLin(clusterIndex, 
                                                   matrix, 
                                                   evd);
//        double modularity = fineTuneOnKernighanLin(clusterIndex, 
//                                                   eigenvector,
//                                                   evd);
        if (modularity <= 0) {
            clusters.add(new HashSet<String>(idsInComponent));
            return;
        }
        List<String> cluster1 = new ArrayList<String>();
        List<String> cluster2 = new ArrayList<String>();
        split(idsInComponent, 
              cluster1,
              cluster2, 
              clusterIndex);
        if (cluster1.size() == 0 || cluster2.size() == 0) {
            // No more clustering needed
            if (cluster1.size() > 0)
                clusters.add(new HashSet<String>(cluster1));
            else
                clusters.add(new HashSet<String>(cluster2));
            return;
        }
        // Work on sub-clusters
        cluster(cluster1, 
                fisInComponent, 
                true,
                clusters);
        cluster(cluster2, 
                fisInComponent, 
                true,
                clusters);
    }
    
    private void split(List<String> idList,
                       List<String> cluster1,
                       List<String> cluster2,
                       int[] clusterIndex) {
        for (int i = 0; i < idList.size(); i++) {
            if (clusterIndex[i] == 1)
                cluster1.add(idList.get(i));
            else
                cluster2.add(idList.get(i));
        }
    }
    
    /**
     * The returned array is used to compare the change of modularity after a node is switched
     * from one cluster to another. It is not needed to calculate the actual values to find the 
     * node that can bring in largest modularity increase or least modularity decrease. What we
     * need is to used the return array and the original cluster index to check the change, as
     * used to find the best clustering using the eigenvector for the biggest positive eigenvalue.
     * @param clusterIndex
     * @param matrix
     * @return
     */
    private double[] calculateLeftProduct(int[] clusterIndex,
                                          double[][] matrix) {
        double[] rtn = new double[clusterIndex.length];
        for (int i = 0; i < rtn.length; i++) {
            double[] vector = matrix[i];
            double tmp = 0.0;
            for (int j = 0; j < vector.length; j++) {
                if (j == i) // Exclude itself. This element should not be used for modularity increase calculation.
                    continue;
                tmp += clusterIndex[j] * vector[j];
            }
            rtn[i] = tmp;
        }
        return rtn;
    }
    
    private void updateLeftProduct(double[] leftProduct,
                                   int[] clusterIndex, // has been updated
                                   double[][] matrix,
                                   int updateIndex) {
        for (int i = 0; i < clusterIndex.length; i++) {
            if (i == updateIndex)
                continue; // This element is not used
            leftProduct[i] += (2 * clusterIndex[updateIndex] * matrix[i][updateIndex]);
        }
    }
    
    private double fineTuneOnKernighanLin(int[] clusterIndex,
                                          DoubleMatrix2D matrix,
                                          EigenvalueDecomposition evd) {
        // Original modularity
        double originalModularity = calculateModularity(evd, clusterIndex);
        // Use array for much faster access of values
        double[][] matrixArray = matrix.toArray();
        while (true) {
            // For changed index
            List<Integer> changedIndex = new ArrayList<Integer>();
            List<Double> modularityIncreases = new ArrayList<Double>();
            double biggestIncrease;
            int maxIndex;
            double testIncrease;
            // To hold test clusterIndex
            int[] testClusterIndex = new int[clusterIndex.length];
            System.arraycopy(clusterIndex, 0,
                             testClusterIndex, 0, 
                             clusterIndex.length);
            List<Integer> toBeChanged = new ArrayList<Integer>();
            for (int i = 0; i < clusterIndex.length; i++)
                toBeChanged.add(i);
            // The double array is used to calculate the difference of modularity.
            double[] leftProduct = calculateLeftProduct(testClusterIndex, matrixArray);
            // The following loop is used to search for a best set of flips
            while (toBeChanged.size() > 0) {
                biggestIncrease = Double.NEGATIVE_INFINITY;
                maxIndex = -1;
                // The double array is used to calculate the difference of modularity.
                //double[] leftProduct = calculateLeftProduct(testClusterIndex, matrixArray);
                // The following loop is used to find the index that should be flip
                for (Integer index : toBeChanged) {
                    testIncrease = -testClusterIndex[index] * leftProduct[index];
                    // Compare the difference to avoid random effects from double float problem.
//                    if (testIncrease - biggestIncrease > ZERO_CUTOFF) {
                    if (testIncrease > biggestIncrease) {
                        biggestIncrease = testIncrease;
                        maxIndex = index;
                    }
                }
                // Make a change
                toBeChanged.remove(new Integer(maxIndex)); // Need to use Integer instead of int.
                testClusterIndex[maxIndex] = -testClusterIndex[maxIndex];
                // Use update instead of re-calculate to increase the speed.
                updateLeftProduct(leftProduct,
                                  testClusterIndex,
                                  matrixArray, 
                                  maxIndex);
                changedIndex.add(maxIndex);
                modularityIncreases.add(biggestIncrease);
//                double tmp = calculateModularity(evd, testClusterIndex);
//                System.out.println("Index " + c);
//                System.out.println("Original: " + originalModularity);
//                System.out.println("modularity after flip: " + tmp);
//                System.out.println("Increased to: " + (originalModularity + 4 * addAll(modularityIncreases)));
//                System.out.println("Increases: " + modularityIncreases);
//                c++;
            }
            // Search the highest modualrity
            biggestIncrease = Double.NEGATIVE_INFINITY;
            maxIndex = -1;
            double incrementalIncrease = 0.0;
            for (int i = 0; i < modularityIncreases.size(); i++) {
                incrementalIncrease += modularityIncreases.get(i);
//                if (incrementalIncrease - biggestIncrease > ZERO_CUTOFF) {
                if (incrementalIncrease > biggestIncrease) {
                    biggestIncrease = incrementalIncrease;
                    maxIndex = i;
                }
            }
            if (biggestIncrease < ZERO_CUTOFF)
                break;
            if (DEBUG) {
                System.out.println("Max index: " + maxIndex);
                System.out.println("Fine tuning clustering...");
                System.out.println("Change modualarity from " + originalModularity + " to " + 
                                   (originalModularity + 4 * biggestIncrease));
                System.out.println("changed index: " + changedIndex);
            }
            // Make changes based on cached changedIndex
            for (int i = 0; i < maxIndex + 1; i++) {
                Integer index = changedIndex.get(i);
                clusterIndex[index] = -clusterIndex[index];
            }
            originalModularity += 4 * biggestIncrease;
            //originalModularity = calculateModularity(evd, clusterIndex);
            //break;
        }
        return originalModularity;
    }
    
    /**
     * Fine tuning the cluster index vector based on Kernighan-Lin algorithm. This implementation 
     * should be deleted. Leave it here for test purpose only.
     * @param clusterIndex
     * @param eigenvector
     * @param evd
     */
    private double fineTuneOnKernighanLin(int[] clusterIndex,
                                          double[] eigenvector,
                                          EigenvalueDecomposition evd) {
        // Original modularity
        double originalModularity = calculateModularity(evd, clusterIndex);
        double[][] eigenmatrix = new Algebra().transpose(evd.getV()).toArray();
        double[] eigenvalues = evd.getRealEigenvalues().toArray();
        while (true) {
            // For changed index
            List<Integer> changedIndex = new ArrayList<Integer>();
            List<Double> changedModularities = new ArrayList<Double>();
            double maxModularity;
            int maxIndex;
            double testModularity;
            // To hold test clusterIndex
            int[] testClusterIndex = new int[clusterIndex.length];
            System.arraycopy(clusterIndex, 0,
                             testClusterIndex, 0, 
                             clusterIndex.length);
            List<Integer> toBeChanged = new ArrayList<Integer>();
            for (int i = 0; i < clusterIndex.length; i++)
                toBeChanged.add(i);
            // The following loop is used to search for a best set of flips
            while (toBeChanged.size() > 0) {
                maxModularity = Double.NEGATIVE_INFINITY;
                maxIndex = -1;
                // The following loop is used to find the index that should be flip
                for (Integer index : toBeChanged) {
                    testModularity = calculateModularity(eigenmatrix,
                                                         eigenvalues, 
                                                         testClusterIndex,
                                                         index);
                    if (testModularity > maxModularity) {
                        maxModularity = testModularity;
                        maxIndex = index;
                    }
                }
                // Make a change
                toBeChanged.remove(new Integer(maxIndex)); // Need to use Integer instead of int.
                testClusterIndex[maxIndex] = -testClusterIndex[maxIndex];
                changedIndex.add(maxIndex);
                changedModularities.add(maxModularity);
            }
            // Search the highest modualrity
            maxModularity = Double.NEGATIVE_INFINITY;
            maxIndex = -1;
            for (int i = 0; i < changedModularities.size(); i++) {
                if (changedModularities.get(i) > maxModularity) {
                    maxModularity = changedModularities.get(i);
                    maxIndex = i;
                }
            }
            if (originalModularity >= maxModularity)
                break;
            if (DEBUG) {
                System.out.println("Max index: " + maxIndex);
                System.out.println("Fine tuning clustering...");
                System.out.println("Change modualarity from " + originalModularity + " to " + maxModularity);
                System.out.println("changed index: " + changedIndex);
            }
            // Make changes based on cached changedIndex
            for (int i = 0; i < maxIndex + 1; i++) {
                Integer index = changedIndex.get(i);
                clusterIndex[index] = -clusterIndex[index];
            }
            originalModularity = maxModularity;
        }
        return originalModularity;
    }
    
    /**
     * Based on the passed eigenvector, create a cluster index vector: 1 for first cluster
     * (element >= 0, note: equal here. This is rather arbitrary), -1 for second cluster (element < 0)
     * @param eigenvector
     * @return
     */
    private int[] generateClusterVector(double[] eigenvector) {
        int[] clusterIndex = new int[eigenvector.length];
        for (int i = 0; i < eigenvector.length; i++) {
            if (eigenvector[i] >= 0)
                clusterIndex[i] = 1;
            else
                clusterIndex[i] = -1;
        }
        return clusterIndex;
    }
    
    /**
     * This method is used to search for the biggest. Based on the test results, the
     * largest eigenvalue should be at the end of eigenvalues list. It seems that 
     * there is no need to check it (Note added by Guanming on Oct 25, 2015).
     * @param evd
     * @return
     */
    private int searchForPositiveBiggestEigenvalue(EigenvalueDecomposition evd) {
        // Largest should be at the end of eigen value array
        DoubleMatrix1D eigenvalues = evd.getRealEigenvalues();
        double value = eigenvalues.get(eigenvalues.size() - 1);
        if (value <= ZERO_CUTOFF)
            return -1;
        return eigenvalues.size() - 1;
    }
    
    private double[] getEigenvector(EigenvalueDecomposition evd, int colIndex) {
        DoubleMatrix2D matrix = evd.getV();
        double[] vector = new double[matrix.rows()];
        for (int i = 0; i < vector.length; i++) {
            vector[i] = matrix.get(i, colIndex);
        }
        return vector;
    }
    
    /**
     * This method is used to convert a set of FIs into an adjacent matrix.
     * It is assumed that all genes in the gene list should be linked together
     * to form a connected graph component. Has not tested the results if these
     * genes are not in the same graph component.
     * @param fis
     * @return
     */
    private DoubleMatrix2D convertFIsToMatrix(List<String> geneList,
                                              Set<String> fis,
                                              boolean isForSubCluster) {
        DoubleMatrix2D matrix = new SparseDoubleMatrix2D(geneList.size(),
                                                        geneList.size());
        for (int i = 0; i < geneList.size(); i++) {
            String gene1 = geneList.get(i);
            for (int j = i; j < geneList.size(); j++) {
                String gene2 = geneList.get(j);
                // Check if gene1 and gene2 has interaction
                String fi = gene1 + "\t" + gene2;
                int value = 0;
                if (fis.contains(fi)) 
                    value = 1;
                matrix.set(i, j, value);
                matrix.set(j, i, value);
            }
        }
        DoubleMatrix2D secondMatrix = convertFIsToSecondMatrix(geneList, fis);
        subtractMatrix(matrix, secondMatrix);
        if (isForSubCluster)
            subtractMatrix(matrix, generateMatrixForSubCluster(geneList, matrix));
        return matrix;
    }
    
    private DoubleMatrix2D generateMatrixForSubCluster(List<String> geneList,
                                                       DoubleMatrix2D original) {
        DoubleMatrix2D matrix = new SparseDoubleMatrix2D(original.rows(), 
                                                         original.columns());
        for (int i = 0; i < geneList.size(); i++) {
            double total = 0.0;
            for (int k = 0; k < original.columns(); k++)
                total += original.get(i, k);
            matrix.set(i, i, total);
        }
        return matrix;
    }
     
    private void subtractMatrix(DoubleMatrix2D original, 
                                DoubleMatrix2D substract) {
        for (int i = 0; i < original.rows(); i++) {
            for (int j = 0; j < original.columns(); j++) {
                double value = original.get(i, j) - substract.get(i, j);
                original.set(i, j, value);
            }
        }
    }
    
    private DoubleMatrix2D convertFIsToSecondMatrix(List<String> geneList,
                                                 Set<String> fis) {
        int m = fis.size();
        DoubleMatrix2D matrix = new SparseDoubleMatrix2D(geneList.size(), 
                                                             geneList.size()); 
        for (int i = 0; i < geneList.size(); i++) {
            String gene1 = geneList.get(i);
            int k1 = idToPartners.get(gene1).size();
            for (int j = i; j < geneList.size(); j++) {
                String gene2 = geneList.get(j);
                int k2 = idToPartners.get(gene2).size();
                double value = (double) k1 * k2 / (2.0 * m);
                matrix.set(i, j, value);
                // Do a symmetric set
                matrix.set(j, i, value);
            }
        }
        return matrix;
    }
    
    public double calculateModualarity(Collection<Set<String>> clusters,
                                       Set<String> fis) {
        return new NetworkModularityCalculator().calculateModularity(clusters, fis);
//        double score = 0.0;
//        int insideEdge = 0;
//        int allEdge = 0;
//        int index = 0;
//        Map<String, Set<String>> idToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
//        for (Set<String> cluster : clusters) {
//            // Get edges in this cluster
//            insideEdge = 0;
//            for (String fi : fis) {
//                index = fi.indexOf("\t");
//                String gene1 = fi.substring(0, index);
//                String gene2 = fi.substring(index + 1);
//                if (cluster.contains(gene1) && cluster.contains(gene2)) {
//                    insideEdge ++;
//                }
//            }
//            int sub = 0;
//            for (String gene1 : cluster) {
//                int k = idToPartners.get(gene1).size();
//                sub += k;
//            }
//            score += (2.0 * insideEdge - sub * sub / (2.0 * fis.size())) / (2.0 * fis.size());
//        }
//        return score;
    }
//    
//    /**
//     * This method is used to test the clustering results from this class.
//     */
//    @Test
//    public void testClustering() throws Exception {
//        // Load interaction
////        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
////        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
////        // Pick up 250 genes randomly
////        Set<String> selectedGenes = MathUtilities.randomSampling(fiGenes, 1000);
////        Set<String> selectedFIs = InteractionUtilities.getFIs(selectedGenes, fis);
//        
////        List<Set<String>> graphComponents = new GraphAnalyzer().calculateGraphComponents(selectedFIs);
////        Set<String> firstGraphComponent = graphComponents.get(0);
////        System.out.println("Total genes in the graph component: " + firstGraphComponent.size());
////        selectedGenes = firstGraphComponent;
////        selectedFIs = InteractionUtilities.getFIs(selectedGenes, fis);
////        
//        FileUtility fu = new FileUtility();
//        //fu.saveInteractions(selectedFIs, "tmp/SpectralTestFIs.txt");
////        Set<String> selectedFIs = fu.loadInteractions("tmp/BigSpectralTestFIs.txt");
//        Set<String> selectedFIs = fu.loadInteractions("tmp/SpectralTestFIs.txt");
//        Set<String> selectedGenes = InteractionUtilities.grepIDsFromInteractions(selectedFIs);
//        System.out.println("Total genes and FIs: " + selectedGenes.size() + ", " + selectedFIs.size());
//        
//        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
//        
//        long time1 = System.currentTimeMillis();
//        List<Set<String>> clusters = clusterAnalyzer.cluster(selectedGenes);
//        long time2 = System.currentTimeMillis();
//        System.out.println("\nTime for edge: " + (time2 - time1));
//        double modualarity = clusterAnalyzer.calculateModularity(clusters, selectedFIs);
//        System.out.println("Modularity from edge betweenness: " + modualarity);
//        for (int i = 0; i < clusters.size(); i++) {
//            System.out.println(i + "\t" + clusters.get(i).size());
//        }
//        
//        long time11 = System.currentTimeMillis();
//        List<Set<String>> newClusters = cluster(selectedFIs);
//        long time21 = System.currentTimeMillis();
//        System.out.println("\nTime for spectral: " + (time21 - time11));
//        double newModualarity = clusterAnalyzer.calculateModularity(newClusters, selectedFIs);
//        System.out.println("Modularity from new implementation: " + newModualarity);
//        for (int i = 0; i < newClusters.size(); i++) {
//            System.out.println(i + "\t" + newClusters.get(i).size());
//        }
//    }
//    
    
    @Test
    public void testCluster() throws Exception {
        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // Randomly pick up 1000 genes
        int randomNumber = 2500;
        Set<String> randomGenes = MathUtilities.randomSampling(fiGenes, randomNumber);
        Set<String> randomFIs = InteractionUtilities.getFIs(randomGenes, fis);
        System.out.println("Total FIs: " + randomFIs.size());
        long time1 = System.currentTimeMillis();
        List<Set<String>> clusters = cluster(randomFIs);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1) + " ms");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            System.out.println(i + "\t" + cluster.size());
        }
    }
    
    /**
     * This method is used to test EVD in a Java matrix lib, MTL.
     */
    @Test
    public void testEVDInMTL() throws Exception {
        // Load interaction
        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // Pick up 250 genes randomly
        Set<String> selectedGenes = MathUtilities.randomSampling(fiGenes, 250);
        Set<String> selectedFIs = InteractionUtilities.getFIs(selectedGenes, fis);
        idToPartners = new BreadthFirstSearch().generateIdToPartnersMap(selectedFIs);
        long mem = Runtime.getRuntime().totalMemory();
        System.out.println("memory before loading (M): " + mem/(1024 * 1024));
        List<String> geneList = new ArrayList<String>(idToPartners.keySet());
        Collections.sort(geneList);
        long time2 = System.currentTimeMillis();
        // Try colt for EVD
        DoubleMatrix2D coltMatrix = convertFIsToMatrix(geneList, 
                                                       selectedFIs, 
                                                       true);
        EigenvalueDecomposition coltEVD = new EigenvalueDecomposition(coltMatrix);
        DoubleMatrix2D v = coltEVD.getV();
        System.out.println("\nEVD from colt:");
        for (int i = 0; i < v.columns(); i++) {
            for (int j = i; j < v.columns(); j++) {
                double total = 0.0;
                for (int k = 0; k < v.rows(); k++)
                    total += v.get(i, k) * v.get(j, k);
                System.out.println(i + ", " + j + ": " + total);
            }
        }
    }
    
//    public static void main(String[] args) {
//        try {
//            new SpectralPartitionNetworkCluster().testClustering();
//        }
//        catch(Exception e) {
//            e.printStackTrace();
//        }
//    }
}
