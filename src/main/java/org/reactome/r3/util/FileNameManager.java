/*
 * Created on Jul 24, 2007
 *
 */
package org.reactome.r3.util;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;

/**
 * This singleton class is used to manage file name used in this r3 project. All file names
 * should be defined in an external file, resources/FileNames.xml.
 * @author guanming
 *
 */
public class FileNameManager {
    private final String MAPPING_FILE_NAME = "resources" + File.separator + "FileNames.xml";
    private Map<String, String> keyToName;
    private static FileNameManager manager;
    
    private FileNameManager() {
        loadFileNames();
    }
    
    public static FileNameManager getManager() {
        if (manager == null)
            manager = new FileNameManager();
        return manager;
    }
    
    
    /**
     * Load the key to file names. Key should be 
     * className + "." + methodName + "." fileName.
     */
    private void loadFileNames() {
//        <fileNames>
//            <dir name="result" value="/Users/wgm/Documents/EclipseWorkspace/caBigR3/results/v2" />
//            <class name="TopicAnalyzer">
//                <method name="deleteMe">
//                    <file dir="result" name="TopicIDNumber.txt" />
//                </method>
//            </class>
//        </fileNames> 
        // Load an XML file using JDOM
        try {
            SAXBuilder builder = new SAXBuilder();
            Document document = builder.build(MAPPING_FILE_NAME);
            Element root = document.getRootElement();
            Map<String, String> dirMap = new HashMap<String, String>();
            keyToName = new HashMap<String, String>();
            List children = root.getChildren();
            for (Iterator it = children.iterator(); it.hasNext();) {
                Element child = (Element) it.next();
                if (child.getName().equals("dir")) {
                    String name = child.getAttributeValue("name");
                    String value = child.getAttributeValue("value");
                    dirMap.put(name, value);
                }
                else if (child.getName().equals("class")) {
                    // Get the class name
                    parseClassElement(child, dirMap);
                }
            }
        }
        catch(Exception e) {
            System.err.println("Cannot load mapping file!");
            throw new RuntimeException("Mapping file is not loaded!");
        }
    }
    
    private void parseClassElement(Element clsElm, Map<String, String> dirMap) {
        String clsName = clsElm.getAttributeValue("name");
        // If a file element is not wrapped around a method, this file name
        // will be a member variable. For example:
//        <file dir="${uniprotDir}" name="uniprot_sprot_human.dat" />
//        <file dir="${uniprotDir}" name="uniprot_trembl_human.dat" />
//        <method name="loadSwissProtIds">
//            <file dir="${uniprotDir}" name="SwissProtIDs.txt" />
//        </method>
        List fileElements = clsElm.getChildren("file");
        parseFileElements(dirMap,
                          clsName,
                          null,
                          fileElements);                             
        List methodElms = clsElm.getChildren("method");
//      <method name="deleteMe">
//          <file dir="${result}" name="TopicIDNumber.txt" />
//      </method>
        for (Iterator it = methodElms.iterator(); it.hasNext();) {
            Element methodElm = (Element) it.next();
            String methodName = methodElm.getAttributeValue("name");
            List fileElms = methodElm.getChildren("file");
            parseFileElements(dirMap, clsName, methodName, fileElms);
        }
    }

    private void parseFileElements(Map<String, String> dirMap, 
                                   String clsName, 
                                   String methodName, 
                                   List fileElms) {
        if (fileElms == null || fileElms.size() == 0)
            return;
        for (Iterator it1 = fileElms.iterator(); it1.hasNext();) {
            Element fileElm = (Element) it1.next();
            String dir = fileElm.getAttributeValue("dir");
            if (dir.startsWith("${")) { // a variable
                dir = dir.substring(2, dir.length() - 1); // Get rid of ${ and }
                dir = dirMap.get(dir);
            }
            String name = fileElm.getAttributeValue("name");
            String fileName = dir + File.separator + name;
            if (methodName == null)
                keyToName.put(clsName + "." + name,
                              fileName);
            else
                keyToName.put(clsName + "." + methodName + "." + name,
                              fileName);
                          
        }
    }
    
    public String getFileName(String fileName) {
        // Need to get the calling method and class
        StackTraceElement[] elements = new Throwable().fillInStackTrace().getStackTrace();
        // The first element is the current method and class. Take the second one
        String clsName = elements[1].getClassName();
        String methodName = elements[1].getMethodName();
        if (methodName.equals("<init>")) // From member variables
            methodName = null;
        String key = null;
        if (methodName == null)
            key = clsName + "." + fileName;
        else
            key = clsName + "." + methodName + "." + fileName;
        String mappedFile = keyToName.get(key);
        if (mappedFile == null)
            throw new IllegalArgumentException("FileNameManager.getFileName(): " + 
                                               fileName + " is not registered.");
        return mappedFile;
    }
}
