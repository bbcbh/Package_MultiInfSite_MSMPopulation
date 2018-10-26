package sim;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;
import opt.Abstract_Optimisation_MSM;
import opt.Optimisation_MSM_TranProb_NumInfFit_GA;
import opt.Optimisation_MSM_TransProb_NumInfFit_NM;

/**
 *
 * @author Ben Hui
 * @version 20181011
 *
 * <pre>
 * History
 * 20180601
 *     - Simplfy directory setting by removing one layer of folder
 * 20181011
 *     - Combine optimisation and simulation into a single simulation interface
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

            sim.setBaseDir(propFile.toFile().getParentFile());
            sim.loadProperties(prop);

            Integer type = (Integer) sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_SIM_TYPE];

            if (type != null && (type == 1 || type == 2) ) {
                // Optimsiation 
                System.out.println("Performing optimisation from " + propFile.toFile().getParentFile().getAbsolutePath() +  " of type = " + type);

                String[] rArg = new String[5];
                // ImportPath
                rArg[0] = propFile.toFile().getParentFile().getAbsolutePath();
                // Target
                rArg[1] = sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_NUM_INF_TARGET]
                        != null ? sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_NUM_INF_TARGET].toString() : "860,830,20";
                // Number of sim 
                rArg[2] = sim.getPropVal()[Simulation_MSM_Population.PROP_NUM_SIM_PER_SET]
                        != null ? sim.getPropVal()[Simulation_MSM_Population.PROP_NUM_SIM_PER_SET].toString() : "";
                // Base seed
                rArg[3] = sim.getPropVal()[Simulation_MSM_Population.PROP_BASESEED]
                        != null ? sim.getPropVal()[Simulation_MSM_Population.PROP_BASESEED].toString() : "";
                // Number of threads to use 
                rArg[4] = sim.getPropVal()[Simulation_MSM_Population.PROP_USE_PARALLEL]
                        != null ? sim.getPropVal()[Simulation_MSM_Population.PROP_USE_PARALLEL].toString() : "";

                Abstract_Optimisation_MSM run;
                if (type == 1) {
                    run = new Optimisation_MSM_TransProb_NumInfFit_NM(rArg);
                } else {
                    run = new Optimisation_MSM_TranProb_NumInfFit_GA(rArg);
                }
                run.runOptimisation();

            } else {
                System.out.println("Generate results set for " + propFile.toFile().getParentFile().getAbsolutePath());
                sim.generateOneResultSet();
                sim.decodeSnapCountFile();
            }

        }

    }

}
