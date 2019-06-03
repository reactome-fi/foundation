/*
 * Created on Jun 29, 2006
 *
 */
package org.reactome.r3.util;

import java.lang.reflect.Method;

/**
 * Class is used to launch in public method in an application environment. This test is 
 * different from UnitTest, which can be used to hold GUI windows..
 * @author guanming
 *
 */
public class Launcher {
    
    public static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: java org.reactome.r3.Launcher clsName methodName");
            System.exit(1);
        }
        try {
            Class cls = Class.forName(args[0]);
            Object obj = cls.newInstance();
            Method method = cls.getMethod(args[1], new Class[]{});
            method.invoke(obj, new Object[]{});
        }
        catch(Exception e) {   
            e.printStackTrace();
        }
    }

}
