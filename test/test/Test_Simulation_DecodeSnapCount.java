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

        File singleSimFolder = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\1000_Runs\\Import_1_BaseModel");//"C:\\Users\\Bhui\\Desktop\\TestDir");

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
        Simulation_MSM_Population.decodeExportedPopulationFiles(singleSimFolder);
    }

}
