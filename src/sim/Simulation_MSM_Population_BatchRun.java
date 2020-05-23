package sim;

import opt.Optimisation_MSM_TranProb_NumInfRangeFit_GA;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;
import opt.Abstract_Optimisation_MSM;
import opt.Optimisation_MSM_TranProb_NumInfFit_GA;
import opt.Optimisation_MSM_TransProb_NumInfFit_NM;
import static sim.Simulation_MSM_Population.FLAG_NO_REMOVAL;

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
        String extraFlag = "";

        System.out.println("BaseDir = " + baseDir.getAbsolutePath());

        if (arg.length < 1) {
            folderNames = new String[]{baseDir.getName()};
            baseDir = baseDir.getParentFile();
        } else {
            ArrayList<File> matchFiles = new ArrayList<>();

            for (int i = 1; i < arg.length; i++) {

                if (!arg[i].startsWith("-")) {

                    final Pattern regEx = Pattern.compile(arg[i]);

                    File[] dirNames = baseDir.listFiles(new FileFilter() {
                        @Override
                        public boolean accept(File file) {
                            return file.isDirectory() && regEx.matcher(file.getName()).matches()
                                    && new File(file, SimulationInterface.FILENAME_PROP).exists();
                        }
                    });

                    matchFiles.addAll(Arrays.asList(dirNames));
                } else {
                    extraFlag = arg[i];
                }

            }

            folderNames = new String[matchFiles.size()];
            for (int i = 0; i < folderNames.length; i++) {
                folderNames[i] = matchFiles.get(i).getName();
            }
        }

        System.out.println(String.format("# folder matched  = %d", folderNames.length));

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
            sim.setSimExtraFlags(extraFlag);

            Integer type = (Integer) sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_SIM_TYPE];

            if (type != null && (type == 1 || type == 2 || type == 3)) {
                // Optimsiation 
                System.out.println("Performing optimisation from " + propFile.toFile().getParentFile().getAbsolutePath() + " of type = " + type);

                String[] rArg = new String[5];
                // ImportPath
                rArg[0] = propFile.toFile().getParentFile().getAbsolutePath();
                // Target
                rArg[1] = sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_CUSTOM_PARAMETER]
                        != null ? sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_CUSTOM_PARAMETER].toString() : "20,830,860";
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

                switch (type) {
                    case 1:
                        run = new Optimisation_MSM_TransProb_NumInfFit_NM(rArg);
                        break;
                    case 2:
                        run = new Optimisation_MSM_TranProb_NumInfFit_GA(rArg);
                        break;
                    default:
                        run = new Optimisation_MSM_TranProb_NumInfRangeFit_GA(rArg);

                }

                run.runOptimisation();

            } else {
                File simBaseDir = new File(baseDir, singleSimFolder);
                File zipOBJ = new File(simBaseDir, "OBJ_Files.zip");

                if (zipOBJ.exists()) {
                    System.out.println("Results Obj zip ("+ zipOBJ.getName()+ ") already existed for " + simBaseDir 
                            + ". Simulation skipped.");

                } else {

                    System.out.println("Generate results set for " + propFile.toFile().getParentFile().getAbsolutePath());

                    if (sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_CUSTOM_PARAMETER] != null) {
                        sim.setSimCustomParameterStr(sim.getPropVal()[Simulation_MSM_Population.PROP_MSM_CUSTOM_PARAMETER].toString());
                    }

                    sim.generateOneResultSet();
                    sim.decodeSnapCountFile();

                    Simulation_MSM_Population.decodeExportedPopulationFiles(new File(baseDir, singleSimFolder), null);

                    // Zipping and clean up
                    boolean zipOkay_DIR = false;
                    boolean zipOkay_CSV = false;
                    boolean zipOkay_OBJ = false;

                    // Zipping all directory
                    File[] exportDirs = simBaseDir.listFiles(new FileFilter() {
                        @Override
                        public boolean accept(File file) {
                            return file.isDirectory();
                        }
                    });
                    try {

                        for (File dir : exportDirs) {
                            util.FileZipper.zipFile(dir, new File(simBaseDir, dir.getName() + ".zip"));
                        }

                        zipOkay_DIR = true;
                    } catch (IOException ex) {
                        ex.printStackTrace(System.err);
                    }

                    // Zipping of CSV file
                    File[] csvFile = simBaseDir.listFiles(new FileFilter() {
                        @Override
                        public boolean accept(File pathname) {
                            return pathname.getName().endsWith(".csv");
                        }
                    });

                    File zipCSV = new File(simBaseDir, "CSV_Files.zip");
                    try (ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(zipCSV))) {
                        for (File csv : csvFile) {
                            ZipEntry zEnt = new ZipEntry(csv.getName());
                            zos.putNextEntry(zEnt);
                            byte[] data = Files.readAllBytes(csv.toPath());
                            zos.write(data, 0, data.length);
                            zos.closeEntry();
                        }

                        zipOkay_CSV = true;
                    } catch (IOException ex) {
                        ex.printStackTrace(System.err);
                    }

                    // Zipping of CSV file
                    File[] objFile = simBaseDir.listFiles(new FileFilter() {
                        @Override
                        public boolean accept(File pathname) {
                            return pathname.getName().endsWith(".obj");
                        }
                    });

                    try (ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(zipOBJ))) {
                        for (File obj : objFile) {
                            ZipEntry zEnt = new ZipEntry(obj.getName());
                            zos.putNextEntry(zEnt);
                            byte[] data = Files.readAllBytes(obj.toPath());
                            zos.write(data, 0, data.length);
                            zos.closeEntry();
                        }

                        zipOkay_OBJ = true;
                    } catch (IOException ex) {
                        ex.printStackTrace(System.err);
                    }

                    if (!sim.getSimExtraFlags().contains(FLAG_NO_REMOVAL)) {

                        File[] removeFile = simBaseDir.listFiles(new FileFilter() {
                            @Override
                            public boolean accept(File pathname) {
                                return pathname.getName().endsWith("_pre");
                            }
                        });

                        for (File prefile : removeFile) {
                            Files.delete(prefile.toPath());

                        }

                        if (zipOkay_CSV) {
                            for (File csv : csvFile) {
                                Files.delete(csv.toPath());
                            }
                        }

                        if (zipOkay_OBJ) {
                            for (File obj : objFile) {
                                Files.delete(obj.toPath());
                            }
                        }

                        if (zipOkay_DIR) {
                            for (File dir : exportDirs) {
                                File[] allFiles = dir.listFiles();
                                for (File f : allFiles) {
                                    Files.delete(f.toPath());
                                }
                                Files.delete(dir.toPath());
                            }
                        }

                    }

                }

            }

        }

    }

}
