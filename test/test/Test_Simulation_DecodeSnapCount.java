package test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Properties;
import sim.Simulation_MSM_Population;

public class Test_Simulation_DecodeSnapCount {

    public static void main(String[] arg) throws IOException, FileNotFoundException, ClassNotFoundException {

        File singleSimFolder = new File(
                                      //"C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Results_Analyse\\Import_1_A_R_B2_10000");                                      
                                      //"C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test\\Import_1_A_R_Pop_100000");
                                      "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\NG Vacc\\Vacc_C10_S000_T050_D720");

        Simulation_MSM_Population sim = new Simulation_MSM_Population();

        Path propFile = new File(singleSimFolder, Simulation_MSM_Population.FILENAME_PROP).toPath();
        Properties prop;
        prop = new Properties();
        try (InputStream inStr = java.nio.file.Files.newInputStream(propFile)) {
            prop.loadFromXML(inStr);
        }

        sim.setBaseDir(propFile.toFile().getParentFile());
        sim.loadProperties(prop);       
        //sim.decodeSnapCountFile();
        Simulation_MSM_Population.decodeExportedPopulationFiles(singleSimFolder, null);
    }

}
