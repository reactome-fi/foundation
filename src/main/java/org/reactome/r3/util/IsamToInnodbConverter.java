/*
 * Created on Feb 23, 2012
 *
 */
package org.reactome.r3.util;

import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import org.junit.Test;

/**
 * This class is used to convert MyISAM tables into InnoDB tables. This class basically is converted
 * from Imre's MyISAMToInnodb with a little modification.
 * @author gwu
 *
 */
public class IsamToInnodbConverter {
    
    public IsamToInnodbConverter() {
    }
    
    private void installDriver() throws Exception {
        String dbDriver = "com.mysql.jdbc.Driver";
        Class.forName(dbDriver).newInstance();
    }
    
    private Connection connect(String dbHost,
                               String dbName,
                               String dbUser,
                               String dbPwd) throws SQLException {
        String connectionStr = "jdbc:mysql://" + dbHost + ":3306/" + dbName;
             //+ "?autoReconnect=true";
        Properties prop = new Properties();
        prop.setProperty("user", dbUser);
        prop.setProperty("password", dbPwd);
        prop.setProperty("useOldUTF8Behavior", "true");
        prop.setProperty("zeroDateTimeBehavior", "convertToNull");
        prop.setProperty("dontTrackOpenResources", "true");
        return DriverManager.getConnection(connectionStr, prop);
    }
    
    /**
     * Convert MyISAM to InnoDB
     * @throws Exception
     */
    @Test
    public void convertMyISAMToInnoDB() throws Exception {
        installDriver();
        Connection connection = connect("localhost", 
                                        "reactome_39_plus_i",
                                        "root", 
                                        "macmysql01");
        DatabaseMetaData metaData = connection.getMetaData();
        ResultSet results = metaData.getTables(null, null, null, new String[]{"TABLE"});
        int count = 0;
        List<String> tables = new ArrayList<String>();
        while (results.next()) {
            System.out.println(results.getString("TABLE_NAME") + "\t" + 
                               results.getString("TABLE_TYPE"));
            count ++;
            String table = results.getString("TABLE_NAME");
            tables.add(table);
        }
        results.close();
        java.sql.Statement stat = connection.createStatement();
        for (String tableName : tables) {
            // Need to drop all fulltext index
            dropFullTextIndices(tableName, stat, connection);
            stat.execute("ALTER TABLE " + tableName + " ENGINE=InnoDB");
        }
        stat.close();
        connection.close();
        System.out.println("Total tables: " + count);
    }
    
    private void dropFullTextIndices(String tableName, 
                                     java.sql.Statement stat,
                                     Connection conn) throws Exception {
        ResultSet results = stat.executeQuery("SHOW INDEX FROM " + tableName);
        Statement dropIndexStat = conn.createStatement();
        while (results.next()) {
            String indexType = results.getString("Index_type");
            if (indexType.equals("FULLTEXT")) {
                String indexName = results.getString("Key_name");
                dropIndexStat.executeUpdate("ALTER TABLE " + tableName + " DROP INDEX " + indexName);
            }
        }
        dropIndexStat.close();
    }
    
}
