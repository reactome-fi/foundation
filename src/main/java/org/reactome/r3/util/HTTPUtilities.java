/*
 * Created on Jun 6, 2017
 *
 */
package org.reactome.r3.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.RequestEntity;
import org.apache.commons.httpclient.methods.StringRequestEntity;

/**
 * @author gwu
 *
 */
public class HTTPUtilities {
    public static final String REQUEST_JSON = "json";
    public static final String REQUEST_XML = "xml";
    public static final String HTTP_POST = "Post";
    public static final String HTTP_GET = "Get";
    
    private static HttpClient initializeHTTPClient(PostMethod post, String query) throws UnsupportedEncodingException {
        RequestEntity entity = new StringRequestEntity(query, "text/plain", "UTF-8");
//        RequestEntity entity = new StringRequestEntity(query, "application/XML", "UTF-8");
        post.setRequestEntity(entity);
//        post.setRequestHeader("Accept", "application/JSON, application/XML, text/plain");
              post.setRequestHeader("Accept", "application/json");
        HttpClient client = new HttpClient();
        return client;
    }
    
    public static String callHttp(String url,
                                  String type,
                                  String query,
                                  String requestType) throws IOException {
        HttpMethod method = null;
        HttpClient client = null;
        if (type.equals(HTTP_POST)) {
            method = new PostMethod(url);
            client = initializeHTTPClient((PostMethod) method, query);
        } else {
            method = new GetMethod(url); // Default
            client = new HttpClient();
        }
        if (requestType.equals(REQUEST_JSON))
            method.setRequestHeader("Accept", "application/json");
        else { // Use XML as the default
            method.setRequestHeader("Accept", "text/plain, application/xml");
        }
        int responseCode = client.executeMethod(method);
        if (responseCode == HttpStatus.SC_OK) {
            InputStream is = method.getResponseBodyAsStream();
            return readMethodReturn(is);
        } else {
            System.err.println("Error from server: " + method.getResponseBodyAsString());
            System.out.println("Response code: " + responseCode);
            throw new IllegalStateException(method.getResponseBodyAsString());
        }
    }
    
    private static String readMethodReturn(InputStream is) throws IOException {
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader reader = new BufferedReader(isr);
        StringBuilder builder = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null)
            builder.append(line).append("\n");
        reader.close();
        isr.close();
        is.close();
        // Remove the last new line
        String rtn = builder.toString();
        // Just in case an empty string is returned
        if (rtn.length() == 0)
            return rtn;
        return rtn.substring(0, rtn.length() - 1);
    }
    
}
