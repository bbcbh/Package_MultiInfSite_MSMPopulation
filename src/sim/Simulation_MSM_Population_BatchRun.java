package sim;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;

/**
 *
 * @author Ben Hui
 * @version 20180601
 *
 * <pre>
 * History
 * 20180601
 *     - Simplfy directory setting by removing one layer of folder
 *
 * </pre>
 */
public class Simulation_MSM_Population_BatchRun {

    public static void main(String[] arg) throws IOException, InterruptedException, FileNotFoundException, ClassNotFoundException {

        Simulation_MSM_Population sim;
        File baseDir = new File(arg[0]);
        String[] folderNames;

        System.out.println("BaseDir = " + baseDir.getAbsolutePath());

        if (arg.length > 1) {
            folderNames = Arrays.copyOfRange(arg, 1, arg.length);
        } else {
            File[] dirNames = baseDir.listFiles(new FileFilter() {
                @Override
                public boolean accept(File file) {
                    return file.isDirectory();
                }
            });

            folderNames = new String[dirNames.length];
            for (int i = 0; i < folderNames.length; i++) {
                folderNames[i] = dirNames[i].getName();
            }
        }

        for (String singleSimFolder : folderNames) {
            sim = new Simulation_MSM_Population();
            File parentDir = new File(baseDir, singleSimFolder);
            Path propFile = new File(parentDir, Simulation_MSM_Population.FILENAME_PROP).toPath();
            Properties prop;
            prop = new Properties();
            try (InputStream inStr = java.nio.file.Files.newInputStream(propFile)) {
                prop.loadFromXML(inStr);
            }

            System.out.println("Generate results set for " + propFile.toFile().getParentFile().getAbsolutePath());
            sim.setBaseDir(propFile.toFile().getParentFile());
            sim.loadProperties(prop);
            sim.generateOneResultSet();
            sim.decodeSnapCountFile();
            //sim.decodeExportedPopulationFiles();

        }

    }

}
