/*
 * Created on May 22, 2017
 *
 */
package org.reactome.r3.util;

import java.awt.Dimension;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.BoxAndWhiskerToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.ui.ApplicationFrame;

/**
 * Plotter based on JFreeChart API
 * @author gwu
 *
 */
public class Plotter extends ApplicationFrame {
    
//    static {
//       ChartFactory.setChartTheme(StandardChartTheme.createJFreeTheme());
//    }
    
    /**
     * Default constructor.
     */
    public Plotter(String title) {
        super(title);
    }
    
    /**
     * Make sure the types and values match each other.
     * @param values
     * @param types
     * @param mainTitle
     * @param xLabel
     * @param yLabel
     */
    public void plotBoxAndWhisker(List<List<Double>> values,
                                  List<String> types,
                                  String mainTitle,
                                  String xLabel,
                                  String yLabel) {
        if (types.size() != values.size())
            throw new IllegalArgumentException("Sizes of values and types don't match!");
        DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();
        for (int i = 0; i < types.size(); i++) {
            List<Double> list = values.get(i);
            String type = types.get(i);
            dataset.add(list, mainTitle, type);
        }
        // The following code is based on http://www.java2s.com/Code/Java/Chart/JFreeChartBoxAndWhiskerDemo.htm
        CategoryAxis xAxis = new CategoryAxis(xLabel);
        NumberAxis yAxis = new NumberAxis(yLabel);
        BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
        renderer.setToolTipGenerator(new BoxAndWhiskerToolTipGenerator());
        CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis, renderer);
        JFreeChart chart = new JFreeChart(mainTitle, plot);
        showChart(chart);
    }
    
    /**
     * Plot two lists of doubles as an XY line.
     * @param xValues
     * @param yValues
     * @param mainTitle
     * @param xLabel
     * @param yLabel
     */
    public void plotScatterPlot(List<Double> xValues,
                                List<Double> yValues,
                                String mainTitle,
                                String xLabel,
                                String yLabel) {
        DefaultXYDataset dataset = new DefaultXYDataset();
        // Convert two lists into a double array
        double[][] data = new double[2][];
        // First array is for x-value
        data[0] = new double[xValues.size()];
        for (int i = 0; i < xValues.size(); i++)
            data[0][i] = xValues.get(i);
        data[1] = new double[yValues.size()];
        for (int i = 0; i < yValues.size(); i++)
            data[1][i] = yValues.get(i);
        dataset.addSeries(mainTitle, data);
        JFreeChart chart = ChartFactory.createScatterPlot(mainTitle,
                xLabel, 
                yLabel, 
                dataset, 
                PlotOrientation.HORIZONTAL,
                true,
                true, 
                false);
        showChart(chart);
    }
    
    /**
     * Plot multiple datasets as a histogram.
     * @param dataSetToValues
     * @param bins
     * @param mainTitle
     * @param xLabel
     * @param yLabel
     */
    public void plotHistograpm(Map<String, List<Double>> dataSetToValues,
                               int bins,
                               String mainTitle,
                               String xLabel,
                               String yLabel) {
        // Get the maximum and minimum values
        double[] minMax = {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        dataSetToValues.values()
        .stream()
        .flatMap(List::stream)
        .distinct()
        .forEach(value -> {
            minMax[0] = Math.min(minMax[0], value);
            minMax[1] = Math.max(minMax[1], value);
        });
        HistogramDataset dataset = new HistogramDataset();
        dataSetToValues.keySet().stream()
                                .sorted((set1, set2) -> {
                                    List<Double> values1 = dataSetToValues.get(set1);
                                    List<Double> values2 = dataSetToValues.get(set2);
                                    return values1.size() - values2.size(); 
                                 })
                                .forEach(name -> {
                                    List<Double> values = dataSetToValues.get(name);
                                    double[] aValues = values.stream().mapToDouble(v -> v.doubleValue()).toArray();
                                    dataset.addSeries(name,
                                                      aValues,
                                                      bins,
                                                      minMax[0],
                                                      minMax[1]);
                                });
        // create the chart...
        JFreeChart chart = ChartFactory.createHistogram(mainTitle,
                                                        xLabel, 
                                                        yLabel, 
                                                        dataset,
                                                        PlotOrientation.VERTICAL, 
                                                        true, 
                                                        true, 
                                                        true);
        showChart(chart);
    }

    private void showChart(JFreeChart chart) {
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setFillZoomRectangle(true);
//        chartPanel.setMouseWheelEnabled(true);
        chartPanel.setPreferredSize(new Dimension(750, 540));
        setContentPane(chartPanel);
        pack();
        setVisible(true);
    }
    
}
