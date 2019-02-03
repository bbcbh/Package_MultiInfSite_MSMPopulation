package sim;

import java.beans.PropertyChangeSupport;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Properties;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.math3.distribution.PoissonDistribution;
import person.AbstractIndividualInterface;
import population.AbstractRegCasRelMapPopulation;
import population.MSMPopulation;
import population.person.MultiSiteMultiStrainPersonInterface;
import population.person.RelationshipPerson_MSM;
import random.RandomGenerator;
import static sim.SinglePopRunnable.EXPORT_PREFIX;
import util.AppendableObjOutstream;
import util.FileZipper;
import util.PersonClassifier;
import util.StaticMethods;
import util.runnable.ExtractFieldRunnable;

/**
 *
 * @author Ben Hui
 *
 * @version 20190123
 *
 * <pre>
 * History
 * 20150528
 *     - Change the role and renaming of snap count
 * 20180530
 *     - Add support for skip thread based on snap count
 * 20180531
 *     - Adjustment snapshot count sampling frequency
 * 20181011
 *     - Combine optimisation and simulation into a single simulation interface
 * 20190114
 *     - Add thread support for decoding population zip file.
 * 20190116
 *     - Add support for exported individual snapshot
 * 20190122
 *  - Add support for survival analysis
 * 20190123
 *  - Add support for skip range and output print
 * 20190204
 *  - Add zipping function for folders generated during simulation 
 * </pre>
 */
public class Simulation_MSM_Population implements SimulationInterface {

    public static final String[] PROP_NAME_MSM_MIS = {
        "PROP_STRAINS_INTRO_AT", "PROP_STRAINS_COEXIST_MAT",
        "PROP_MSM_SIM_TYPE", "PROP_MSM_CUSTOM_PARAMETER",
        "PROP_MSM_SKIP_THREAD_RANGE",};

    public static final Class[] PROP_CLASS_MSM_MIS = {
        float[][].class, // float[]{globaltime, strainNum, site, likelihood to co-exist, number of infection to introduce, frequency (optional) }
        float[][].class, // float[exist_strain][new_strain]{likelihood to coexist}
        Integer.class, // 0 (default) = simulation, 1 = optimisation
        String.class, // For optimisation infection targer      
        int[].class,};

    public static final int PROP_STRAINS_INTRO_AT = PROP_NAME.length;
    public static final int PROP_STRAINS_COEXIST_MAT = PROP_STRAINS_INTRO_AT + 1;
    public static final int PROP_MSM_SIM_TYPE = PROP_STRAINS_COEXIST_MAT + 1;
    public static final int PROP_MSM_CUSTOM_PARAMETER = PROP_MSM_SIM_TYPE + 1;
    public static final int PROP_MSM_SKIP_THREAD_RANGE = PROP_MSM_CUSTOM_PARAMETER + 1;

    // Output filename
    public static final String[] FILE_NAMES_OBJ = {"endNumInf.obj", "extinctAt.obj", "snapCount.obj",
        "eventPt.obj", "incidenceCount.obj"};
    public static final int FILE_END_NUM_INF = 0;
    public static final int FILE_EXTINCT_AT = FILE_END_NUM_INF + 1;
    public static final int FILE_SNAPCOUNTS = FILE_EXTINCT_AT + 1;
    public static final int FILE_EVENT_POINTER = FILE_SNAPCOUNTS + 1;
    public static final int FILE_INCIDENT_COUNT = FILE_EVENT_POINTER + 1;

    public static final String[] FILE_NAMES_CSV = {"endNumInfPerson.csv", "infStatSnapshot.csv",
        "numPartnersInlast6Months.csv", "newStrainsHasRegPartners.csv", "strainCompositionActiveRange.csv"};
    public static final int FILE_END_NUM_INF_PERSON_CSV = 0;
    public static final int FILE_INFECTION_STAT_CSV = FILE_END_NUM_INF_PERSON_CSV + 1;
    public static final int FILE_NUM_PARTERS_IN_LAST_6_MONTHS = FILE_INFECTION_STAT_CSV + 1;
    public static final int FILE_NEW_STRAIN_HAS_REG_PARTNERS = FILE_NUM_PARTERS_IN_LAST_6_MONTHS + 1;
    public static final int FILE_STRAIN_COMPOSITION_ACTIVE_RANGE = FILE_NEW_STRAIN_HAS_REG_PARTNERS + 1;

    public static final String[] DIR_NAMES = {"output", "newStrainSpread"};
    public static final int DIR_NAMES_OUTPUT = 0;
    public static final int DIR_NEW_STRAIN_SPREAD = DIR_NAMES_OUTPUT + 1;

    public static final String POP_PROP_INIT_PREFIX = "POP_PROP_INIT_PREFIX_";

    protected Object[] propVal = new Object[PROP_NAME.length + PROP_NAME_MSM_MIS.length];
    protected String[] propModelInit = null;

    protected File baseDir = new File("");
    protected boolean stopNextTurn = false;
    // Snapshot
    private PersonClassifier[] snapshotCountClassifier = null;
    private boolean[] snapshotCountAccum = null;
    // PropertyChangeSupport
    protected PropertyChangeSupport progressSupport = null;
    // Number progress to use
    private int maxThreads = Math.max(Runtime.getRuntime().availableProcessors(), 1);
    // Prextract 
    private Object[][] preExtractField = null;

    private boolean useImportIOThread = true;

    public final static Pattern Pattern_importFile = Pattern.compile("pop_(\\d+).zip");
    public final static Pattern Pattern_popStat = Pattern.compile(SinglePopRunnable.EXPORT_INDIV_PREFIX + "(\\d+).csv");

    private String simCustomParameterStr = null;

    public String getSimCustomParameterStr() {
        return simCustomParameterStr;
    }

    public void setSimCustomParameterStr(String simCustomParameterStr) {
        this.simCustomParameterStr = simCustomParameterStr;
    }

    public void setUseImportIOThread(boolean useImportIOThread) {
        this.useImportIOThread = useImportIOThread;
    }

    @Override
    public void loadProperties(Properties prop) {
        for (int i = 0; i < PROP_NAME.length; i++) {
            String ent = prop.getProperty(PROP_NAME[i]);
            if (ent != null) {
                propVal[i] = StaticMethods.propStrToObject(ent, PROP_CLASS[i]);
            }
        }
        for (int i = PROP_NAME.length; i < propVal.length; i++) {
            String ent = prop.getProperty(PROP_NAME_MSM_MIS[i - PROP_NAME.length]);
            if (ent != null) {
                propVal[i] = StaticMethods.propStrToObject(ent, PROP_CLASS_MSM_MIS[i - PROP_NAME.length]);
            }
        }

        int maxFieldNum = 0;
        for (Iterator<Object> it = prop.keySet().iterator(); it.hasNext();) {
            String k = (String) it.next();
            if (k.startsWith(POP_PROP_INIT_PREFIX)) {
                if (prop.getProperty(k) != null) {
                    maxFieldNum = Math.max(maxFieldNum,
                            Integer.parseInt(k.substring(POP_PROP_INIT_PREFIX.length())));
                }
            }
        }
        if (maxFieldNum >= 0) {
            propModelInit = new String[maxFieldNum + 1];
            for (int i = 0; i < propModelInit.length; i++) {
                String res = prop.getProperty(POP_PROP_INIT_PREFIX + i);
                if (res != null) {
                    propModelInit[i] = res;
                }
            }
        }
    }

    @Override
    public Properties generateProperties() {
        Properties prop = new Properties();
        for (int i = 0; i < PROP_NAME.length; i++) {
            prop.setProperty(PROP_NAME[i], StaticMethods.objectToPropStr(propVal[i], PROP_CLASS[i]));
        }
        for (int i = PROP_CLASS.length; i < propVal.length; i++) {
            prop.setProperty(PROP_NAME_MSM_MIS[i - PROP_NAME.length],
                    StaticMethods.objectToPropStr(propVal[i], PROP_CLASS_MSM_MIS[i - PROP_CLASS.length]));
        }

        return prop;
    }

    public Object[][] getPreExtractField() {
        return preExtractField;
    }

    public void setPreExtractField(Object[][] preExtractField) {
        this.preExtractField = preExtractField;
    }

    @Override
    public void setBaseDir(File baseDir) {
        this.baseDir = baseDir;
    }

    @Override
    public void setStopNextTurn(boolean stopNextTurn) {
        this.stopNextTurn = stopNextTurn;
    }

    @Override
    public void setSnapshotSetting(PersonClassifier[] snapshotCountClassifier, boolean[] snapshotCountAccum) {
        this.snapshotCountClassifier = snapshotCountClassifier;
        this.snapshotCountAccum = snapshotCountAccum;
    }

    @Override
    public void generateOneResultSet() throws IOException, InterruptedException {
        int simSoFar = 0;

        int numProcess = maxThreads;

        if (Integer.parseInt(propVal[PROP_USE_PARALLEL].toString()) > 0) {
            numProcess = Math.min(Integer.parseInt(propVal[PROP_USE_PARALLEL].toString()),
                    numProcess);
        }

        int threadCounter = 0;
        int numSimTotal;
        RandomGenerator rng;

        int[] init_expose = null;

        int[][] eventPointers = null;

        if (MSMPopulation.class.getName().equals(propVal[PROP_POP_TYPE])) {
            snapshotCountClassifier = new PersonClassifier[]{
                new CLASSIFIER_PREVAL(-1),
                new CLASSIFIER_PREVAL(RelationshipPerson_MSM.SITE_G),
                new CLASSIFIER_PREVAL(RelationshipPerson_MSM.SITE_A),
                new CLASSIFIER_PREVAL(RelationshipPerson_MSM.SITE_R)
            };
            snapshotCountAccum = new boolean[snapshotCountClassifier.length];

            init_expose = new int[]{3, 12 * 30, 12 * 7};

        } else {
            throw new UnsupportedOperationException(getClass().getName()
                    + ".generateOneResultSet: Population class "
                    + propVal[PROP_POP_TYPE] + " not supported yet");
        }

        numSimTotal = ((Number) propVal[PROP_NUM_SIM_PER_SET]).intValue();

        int rngCallCounter = 0;
        rng = new random.MersenneTwisterRandomGenerator(((Number) propVal[PROP_BASESEED]).longValue());

        boolean usingImportPop = propVal[PROP_POP_IMPORT_PATH] != null && !((String) propVal[PROP_POP_IMPORT_PATH]).isEmpty();
        File impDir = null;
        File[] impPopFiles = null;

        // Prextract 
        if (usingImportPop) {
            if (preExtractField != null) {
                // Ensure the length of preExtractField matches with number of simulation runs
                preExtractField = Arrays.copyOf(preExtractField, numSimTotal);

            } else {
                preExtractField = new Object[numSimTotal][];
            }

            // Extract fields when required
            impDir = new File((String) propVal[PROP_POP_IMPORT_PATH]);
            impPopFiles = new File[numSimTotal];

            for (int tI = 0; tI < preExtractField.length; tI++) {
                if (preExtractField[tI] == null) {
                    // Extract fields from file                    
                    File zipFile = new File(impDir, SimulationInterface.POP_FILE_PREFIX + tI + ".zip");
                    if (zipFile.exists()) {
                        impPopFiles[tI] = zipFile;
                    }
                }
            }
        }

        File preSnapFile = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS]);

        int simToSkip = 0;
        if (preSnapFile.exists()) {

            ObjectInputStream inStr = new ObjectInputStream(new FileInputStream(preSnapFile));

            try {

                try {
                    boolean containNull = false;

                    while (!containNull) {
                        int[][][] ent = (int[][][]) inStr.readObject();
                        containNull = checkForNullArray(ent);
                        if (!containNull) {
                            simToSkip++;
                        }
                    }
                } catch (EOFException | NullPointerException ex) {
                    inStr.close();
                }
                System.out.println("Number of vaild snapshot in " + preSnapFile.getAbsolutePath() + " = " + simToSkip);
            } catch (ClassNotFoundException ex) {
                ex.printStackTrace(System.err);

            }

        }

        //System.out.println("7: " + propModelInit[7]);
        //System.out.println("8: " + propModelInit[8]);
        while (simSoFar < numSimTotal && !stopNextTurn) {
            int numThreads = Math.min(numProcess, numSimTotal - simSoFar);
            if (threadCounter < simToSkip) {
                long skipSeed = rng.nextLong();
                showStrStatus("Simulation #" + threadCounter + " with seed " + skipSeed
                        + " skipped as "
                        + simToSkip + " snapshot results already present at "
                        + preSnapFile.getAbsolutePath());

                simSoFar++;
                threadCounter++;
            } else if (propVal[PROP_MSM_SKIP_THREAD_RANGE] != null
                    && ((int[]) propVal[PROP_MSM_SKIP_THREAD_RANGE])[0] <= threadCounter
                    && ((int[]) propVal[PROP_MSM_SKIP_THREAD_RANGE])[1] >= threadCounter) {
                long skipSeed = rng.nextLong();
                showStrStatus("Simulation #" + threadCounter + " with seed " + skipSeed
                        + " skipped as threadId within PROP_MSM_SKIP_THREAD_RANGE of "
                        + Arrays.toString((int[]) propVal[PROP_MSM_SKIP_THREAD_RANGE]));

                simSoFar++;
                threadCounter++;

            } else {

                ExecutorService executor;
                SinglePopRunnable[] runnable = null;
                PrintWriter[] priWri = null;
                File outputDirFile = new File(baseDir, DIR_NAMES[DIR_NAMES_OUTPUT]);
                outputDirFile.mkdirs();

                if (usingImportPop) {
                    File[] importPopFileSim = Arrays.copyOfRange(impPopFiles, simSoFar, simSoFar + numThreads);

                    int importPopSoFar = 0;

                    while (importPopSoFar < importPopFileSim.length) {

                        // Extract field thread
                        ExecutorService ioExecutor;
                        ExtractFieldRunnable[] extThreads;
                        java.util.concurrent.Future<Object[]>[] futureFields;
                        File[] extractPopFilesSet = new File[importPopFileSim.length];

                        futureFields = new java.util.concurrent.Future[importPopFileSim.length];
                        ioExecutor = Executors.newFixedThreadPool(numProcess);
                        extThreads = new ExtractFieldRunnable[Math.min(numThreads, importPopFileSim.length - importPopSoFar)];

                        for (int tI = 0; tI < extThreads.length; tI++) {
                            if (importPopFileSim[importPopSoFar + tI] != null) {
                                extractPopFilesSet[tI] = util.FileZipper.unzipFile(importPopFileSim[importPopSoFar + tI], baseDir);
                                extThreads[tI] = new ExtractFieldRunnable(extractPopFilesSet[tI]);
                                if (numThreads > 1 && ((Integer) propVal[PROP_USE_PARALLEL]) != 0 && useImportIOThread) {
                                    System.out.println("Submitting extThread from " + importPopFileSim[importPopSoFar + tI].getAbsolutePath());
                                    futureFields[importPopSoFar + tI] = ioExecutor.submit(extThreads[tI]);
                                } else {
                                    System.out.println("Running extThread from " + importPopFileSim[importPopSoFar + tI].getAbsolutePath());
                                    extThreads[tI].run();
                                }
                            }
                        }
                        ioExecutor.shutdown();

                        try {
                            if (!ioExecutor.awaitTermination(1, TimeUnit.DAYS)) {
                                System.err.println("Thread time-out in extracting fields!");
                            }
                            for (int tI = 0; tI < extThreads.length; tI++) {
                                if (extThreads[tI] != null) {
                                    if (futureFields[importPopSoFar + tI] != null) {
                                        preExtractField[simSoFar + importPopSoFar + tI] = futureFields[importPopSoFar + tI].get();
                                    } else {
                                        preExtractField[simSoFar + importPopSoFar + tI] = extThreads[tI].getExtFields();
                                    }
                                }
                            }

                        } catch (InterruptedException | ExecutionException ex) {
                            StringWriter str = new StringWriter();
                            try (PrintWriter wri = new PrintWriter(str)) {
                                ex.printStackTrace(wri);
                            }
                            System.err.println(str.toString());
                        }
                        importPopSoFar += extThreads.length;

                        for (File extractPopFile : extractPopFilesSet) {
                            if (extractPopFile != null) {
                                extractPopFile.delete();
                            }
                        }
                    }

                }

                showStrStatus("Running S" + threadCounter + " to S" + (threadCounter + numThreads - 1) + "...");

                if (runnable == null || runnable.length != numThreads) {
                    runnable = new SinglePopRunnable[numThreads];
                }
                priWri = new PrintWriter[numThreads];

                executor = Executors.newFixedThreadPool(numThreads);

                for (int r = 0; r < runnable.length; r++) {

                    runnable[r] = new SinglePopRunnable(threadCounter,
                            ((Number) propVal[PROP_NUM_SNAP]).intValue(), ((Number) propVal[PROP_SNAP_FREQ]).intValue());

                    runnable[r].setBaseDir(baseDir);

                    // Set output 
                    priWri[r] = new PrintWriter(new FileWriter(new File(outputDirFile, DIR_NAMES[DIR_NAMES_OUTPUT] + "_" + threadCounter + ".txt")));

                    final PrintWriter pri_Output = priWri[r];
                    runnable[r].setProgressSupport(new PropertyChangeSupport(runnable) {

                        @Override
                        public void firePropertyChange(String key, Object notUsed, Object str) {
                            if (SimulationInterface.PROGRESS_MSG.equals(key)) {
                                pri_Output.println(str);
                                //System.out.println(str);
                            }
                        }

                    });

                    if (propVal[PROP_POP_EXPORT_AT] != null) {
                        runnable[r].setExportBurnInPop((int[]) propVal[PROP_POP_EXPORT_AT]);
                    }

                    if (propVal[PROP_STRAINS_INTRO_AT] != null) {
                        // Generate strain intro file file                 
                        float[][] allEnt = (float[][]) propVal[PROP_STRAINS_INTRO_AT];
                        for (float[] ent : allEnt) {
                            runnable[r].addStrainIntroEnt(ent);
                        }
                    }

                    if (propVal[PROP_STRAINS_COEXIST_MAT] != null) {
                        runnable[r].setCoexistMat((float[][]) propVal[PROP_STRAINS_COEXIST_MAT]);
                    }

                    runnable[r].setPrintPrevalence(((Integer) propVal[PROP_USE_PARALLEL]) <= 1);

                    threadCounter++;

                    if (progressSupport != null) {
                        runnable[r].setProgressSupport(progressSupport);
                    }

                    // New Pop                
                    if (MSMPopulation.class.getName().equals(propVal[PROP_POP_TYPE])) {
                        runnable[r].setPopulation(new MSMPopulation(rng.nextLong()));
                    } else {
                        throw new UnsupportedOperationException(getClass().getName()
                                + ".generateOneResultSet: Population class "
                                + propVal[PROP_POP_TYPE] + " not supported yet");
                    }

                    boolean useImport = false;
                    if (propVal[PROP_POP_IMPORT_PATH] != null) {
                        useImport = populationImport(runnable, r);
                    }

                    // Common snapshot count (if any)
                    runnable[r].setSnapShotOutput(snapshotCountClassifier, snapshotCountAccum);

                    if (!useImport) {
                        System.out.println("Thread #" + (threadCounter - 1) + " generated with seed of " + runnable[r].getPopulation().getSeed());
                        rngCallCounter++;
                        runnable[r].model_prop_initialise(((Number) propVal[PROP_BURNIN]).intValue(), propModelInit);

                    }
                    if (propVal[PROP_INFECTION_INTRO] != null) {

                        float[][] preval_intro = (float[][]) propVal[PROP_INFECTION_INTRO];

                        for (int infId = 0; infId < preval_intro.length; infId++) {
                            for (int t = 0; t < preval_intro[infId].length; t += 2) {
                                runnable[r].setInfectionIntroAt((int) preval_intro[infId][t],
                                        infId, preval_intro[infId][t + 1], init_expose[infId]);
                            }

                        }

                    }

                    if (eventPointers != null) {
                        eventPointers[runnable[r].getId()] = runnable[r].getEventsPointer();
                    }

                    if (getSimCustomParameterStr() != null) {
                        if (getSimCustomParameterStr().contains("Survival_Analysis")) {
                            runnable[r].set_patient_zero(true);
                        }
                    }

                    if (((Integer) propVal[PROP_USE_PARALLEL]) != 0) {
                        executor.submit(runnable[r]);
                    } else {
                        runnable[r].run();
                    }
                }

                executor.shutdown();
                if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
                    showStrStatus("Thread time-out!");
                }

                if (progressSupport != null) {
                    progressSupport.firePropertyChange(PROGRESS_SIM_STORED, simSoFar, simSoFar + runnable.length);
                }
                simSoFar += runnable.length;

                genertateOutputFiles(runnable);

                for (PrintWriter wri : priWri) {
                    wri.close();
                }

            }

        }

        finalise(simSoFar);

    }

    private static boolean checkForNullArray(Object arr) {
        if (arr == null) {
            return true;
        } else {

            if (arr.getClass().isArray()) {
                if (arr instanceof Object[]) {
                    Object[] arrObj = (Object[]) arr;
                    return checkForNullArray(arrObj[arrObj.length - 1]);
                } else {
                    return false;
                }
            } else {
                return arr != null;
            }
        }
    }

    protected int[] generateStrainIntroAt(float[][] ent, int i, RandomGenerator introRNG) {
        // ent[i] : {if +ive number or 1/per snap freq if -ive, stating at}
        //          {if +ive number or 1/per snap freq if -ive, stating at, count limit if +ive, snapshot time limit if -ive}
        int[] introAt = new int[((Number) propVal[PROP_NUM_SNAP]).intValue()];
        float numIntro = ent[i][0];
        int introTime = ((int) ent[i][1]);
        float limit = Float.POSITIVE_INFINITY;
        if (ent[i].length >= 3 && ent[i][2] > 0) {
            limit = ent[i][2];
        }
        if (numIntro > 0) {
            introAt[introTime] = (int) numIntro;
        } else {
            // Randomly based on Poisson dist
            PoissonDistribution poi = new PoissonDistribution(introRNG, -1f / numIntro,
                    PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
            int counter = 0;

            for (int st = introTime; st < introAt.length && counter < limit; st++) {
                introAt[st] = poi.sample();

                if (introAt[st] > 0 && limit != Float.POSITIVE_INFINITY) {
                    introAt[st] = Math.min(introAt[st], ((int) limit) - counter);
                    counter++;
                }

                if (ent[i].length >= 3 && ent[i][2] < 0 && (st - introTime) >= -ent[i][2]) {
                    limit = -1; // Exit at next turn
                }
            }
        }
        return introAt;
    }

    public int getMaxThreads() {
        return maxThreads;
    }

    public void setMaxThreads(int maxThreads) {
        this.maxThreads = maxThreads;
    }

    public PropertyChangeSupport getProgressSupport() {
        return progressSupport;
    }

    public void setProgressSupport(PropertyChangeSupport progressSupport) {
        this.progressSupport = progressSupport;
    }

    protected void showStrStatus(String str) {
        if (progressSupport == null) {
            System.out.println(str);
        } else {
            progressSupport.firePropertyChange(PROGRESS_MSG, null, str);
        }
    }

    protected void finalise(int simSoFar) throws IOException {
        try {
            StaticMethods.decodeResultObjFile(new File(baseDir, FILE_NAMES_OBJ[FILE_END_NUM_INF]), simSoFar);
            StaticMethods.decodeResultObjFile(new File(baseDir, FILE_NAMES_OBJ[FILE_EXTINCT_AT]), simSoFar);

            // Zipping all directory
            File[] exportDirs = baseDir.listFiles(new FileFilter() {
                @Override
                public boolean accept(File file) {
                    return file.isDirectory();
                }
            });

            for (File dir : exportDirs) {
                util.FileZipper.zipFile(dir, new File(baseDir, dir.getName() + ".zip"));
            }
            
            

        } catch (ClassNotFoundException ex) {
            System.err.println(getClass().getName() + ".generateOneResultSet: Error - corrupted data file");
        }
        if (progressSupport
                != null) {
            progressSupport.firePropertyChange(PROGRESS_ALL_DONE, null, simSoFar);
        }

    }

    protected void genertateOutputFiles(SinglePopRunnable[] runnableCollection) throws IOException {
        // Generating file
        ObjectOutputStream[] objS;

        objS = new ObjectOutputStream[FILE_NAMES_OBJ.length];

        for (int i = 0; i < objS.length; i++) {
            File f = new File(baseDir, FILE_NAMES_OBJ[i]);
            if (f.exists()) {
                Files.copy(f.toPath(), new File(baseDir, FILE_NAMES_OBJ[i] + "_pre").toPath(), StandardCopyOption.REPLACE_EXISTING);
            }
            objS[i] = AppendableObjOutstream.generateFromFile(f);
        }

        PrintWriter pri_numInfectPerson = new PrintWriter(new FileWriter(new File(baseDir, FILE_NAMES_CSV[FILE_END_NUM_INF_PERSON_CSV]), true));
        PrintWriter pri_strainCompositionStat = new PrintWriter(new FileWriter(new File(baseDir, FILE_NAMES_CSV[FILE_STRAIN_COMPOSITION_ACTIVE_RANGE]), true));

        for (SinglePopRunnable runnable : runnableCollection) {
            if (runnable != null) {
                objS[FILE_END_NUM_INF].writeObject(runnable.getPopulation().getNumInf());
                objS[FILE_END_NUM_INF].flush();
                objS[FILE_EXTINCT_AT].writeObject(runnable.getExtinctionAt());
                objS[FILE_EXTINCT_AT].flush();
                if (objS[FILE_SNAPCOUNTS] != null) {
                    int[][][] snapcount = runnable.getSnapCounts();
                    objS[FILE_SNAPCOUNTS].writeObject(snapcount);
                    objS[FILE_SNAPCOUNTS].flush();
                    System.out.println("S" + runnable.getId() + "'s end snap = " + Arrays.deepToString(snapcount[snapcount.length - 1]));

                }
                if (objS[FILE_EVENT_POINTER] != null) {
                    objS[FILE_EVENT_POINTER].writeObject(runnable.getEventsPointer());
                    objS[FILE_EVENT_POINTER].flush();
                }
                // PopSnapObj - Bug if there is more than one null in sim?
                int[][] incidentBySnapCount = runnable.getIncidentCounts();
                for (int i = 0; i < incidentBySnapCount.length; i++) {
                    if (incidentBySnapCount[i] == null) {
                        incidentBySnapCount[i] = new int[0];
                    }
                }

                objS[FILE_INCIDENT_COUNT].writeObject(incidentBySnapCount);
                objS[FILE_INCIDENT_COUNT].flush();

                // Number of infected person
                int[] sel = new int[]{-1, 0, 1, 2};
                int[] numInfected = runnable.getPopulation().getNumberOfInfected(sel);

                pri_numInfectPerson.print(runnable.getId());
                pri_numInfectPerson.print(',');
                pri_numInfectPerson.print(runnable.getPopulation().getPop().length);
                for (int i = 0; i < sel.length; i++) {
                    pri_numInfectPerson.print(',');
                    pri_numInfectPerson.print(numInfected[i]);
                }
                pri_numInfectPerson.println();

                int[][] strainComStat = runnable.getStrainCompositionActiveRange();
                pri_strainCompositionStat.print(runnable.getId());

                for (int sC = 1; sC < strainComStat.length; sC++) {
                    pri_strainCompositionStat.print(',');
                    pri_strainCompositionStat.print(sC);
                    for (int i = 0; i < strainComStat[sC].length; i++) {
                        pri_strainCompositionStat.print(',');
                        pri_strainCompositionStat.print(strainComStat[sC][i]);
                    }
                }
                pri_strainCompositionStat.println();

            }

        }
        for (ObjectOutputStream objS1 : objS) {
            if (objS1 != null) {
                objS1.close();
            }
        }
        pri_numInfectPerson.close();
        pri_strainCompositionStat.close();

    }

    public void decodeSnapCountFile() throws FileNotFoundException, IOException, ClassNotFoundException {
        decodeSnapCountFile(baseDir, propVal);
    }

    public static void decodeSnapCountFile(File baseDir, Object[] propVal) throws FileNotFoundException, IOException, ClassNotFoundException {
        File snapCountFile = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS]);

        if (snapCountFile.exists()) {

            ArrayList<int[][][]> snapcountBySim = new ArrayList<>();
            ObjectInputStream inStr = new ObjectInputStream(new FileInputStream(snapCountFile));

            try {
                while (true) {
                    int[][][] ent = (int[][][]) inStr.readObject();

                    if (!checkForNullArray(ent)) {
                        snapcountBySim.add(ent);
                    }
                }
            } catch (EOFException ex) {
                inStr.close();
            }
            System.out.println("Number of snapshot read = " + snapcountBySim.size());

            int[][][][] entCollection = snapcountBySim.toArray(new int[snapcountBySim.size()][][][]);

            int simCounter = 0;
            int[] newStrainExinctAtSnapshot = new int[entCollection.length];
            int[][][] prevalenceByStrainAt = new int[entCollection.length][4][3]; // 1 yr, 2 yr, 5 yr and 10 yr, no infection, base strain, import strain
            int[][][][] yearly_prevalenceBySiteAndStrainAt
                    = new int[entCollection.length][11][4][3]; // [year][site]{no infection, base strain, import strain}

            //Format: from snapshotCountClassifier
            final int PREVAL_ALL = 0;

            for (int[][][] entry : entCollection) {
                int importAtSnap = -1;

                int snapFreq = (int) propVal[PROP_SNAP_FREQ];
                final int snapPerYr = 360 / snapFreq;

                if (propVal[PROP_STRAINS_INTRO_AT] != null) {
                    float[] inFirstStrain = ((float[][]) propVal[PROP_STRAINS_INTRO_AT])[0];
                    int importTime = (int) inFirstStrain[0];
                    int burnIn = propVal[PROP_BURNIN] == null ? 0 : (int) propVal[PROP_BURNIN];
                    importAtSnap = ((importTime - burnIn) / snapFreq) - 1;
                }

                for (int snapCounter = 0; snapCounter < entry.length; snapCounter++) {
                    if (entry[snapCounter][PREVAL_ALL][2] + entry[snapCounter][PREVAL_ALL][3] > 0) {
                        newStrainExinctAtSnapshot[simCounter]++;
                    }
                    if (importAtSnap > 0) {
                        int timePt = -1;
                        int yearTimePt = -1;
                        if (snapCounter - importAtSnap == snapPerYr) {
                            timePt = 0;
                        } else if (snapCounter - importAtSnap == 2 * snapPerYr) {
                            timePt = 1;
                        } else if (snapCounter - importAtSnap == 5 * snapPerYr) {
                            timePt = 2;
                        } else if (snapCounter - importAtSnap == 10 * snapPerYr) {
                            timePt = 3;
                        }

                        if ((snapCounter >= importAtSnap) && (snapCounter - importAtSnap) % snapPerYr == 0) {
                            yearTimePt = (snapCounter - importAtSnap) / snapPerYr;
                        }

                        if (timePt >= 0) {
                            prevalenceByStrainAt[simCounter][timePt][0] = entry[snapCounter][PREVAL_ALL][0];
                            prevalenceByStrainAt[simCounter][timePt][1] = entry[snapCounter][PREVAL_ALL][1];
                            prevalenceByStrainAt[simCounter][timePt][2] = entry[snapCounter][PREVAL_ALL][2] + entry[snapCounter][PREVAL_ALL][3];
                        }

                        if (yearTimePt >= 0 && yearTimePt < yearly_prevalenceBySiteAndStrainAt[simCounter].length) {
                            for (int site = 0; site < entry[snapCounter].length; site++) {
                                yearly_prevalenceBySiteAndStrainAt[simCounter][yearTimePt][site][0] = entry[snapCounter][site][0];
                                yearly_prevalenceBySiteAndStrainAt[simCounter][yearTimePt][site][1] = entry[snapCounter][site][1];
                                yearly_prevalenceBySiteAndStrainAt[simCounter][yearTimePt][site][2] = entry[snapCounter][site][2] + entry[snapCounter][site][3];
                            }
                        }

                    }

                }
                simCounter++;
            }
            // Printing CSV
            PrintWriter pri;

            // Extinct at
            File newStrainExinctAt = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS] + "_newStrainExinctAtSnapshot.csv");
            pri = new PrintWriter(newStrainExinctAt);
            for (int s = 0; s < newStrainExinctAtSnapshot.length; s++) {
                pri.print(s);
                pri.print(',');
                pri.print(newStrainExinctAtSnapshot[s]);
                pri.println();
            }
            pri.close();

            //Prevalence at select time 
            File prevalenceAtSelTime = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS] + "_prevalenceAtSelectedTime.csv");
            pri = new PrintWriter(prevalenceAtSelTime);
            for (int s = 0; s < prevalenceByStrainAt.length; s++) {
                pri.print(s);
                for (int timePt = 0; timePt < prevalenceByStrainAt[s].length; timePt++) {
                    for (int strainType = 0; strainType < prevalenceByStrainAt[s][timePt].length; strainType++) {
                        pri.print(',');
                        pri.print(prevalenceByStrainAt[s][timePt][strainType]);
                    }
                }
                pri.println();
            }

            pri.close();

            // Yearly prevalene at time
            File prevalenceYearly_Any = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS] + "_prevalenceYearly_Any.csv");
            PrintWriter pri_Any = new PrintWriter(prevalenceYearly_Any);

            File prevalenceYearly_G = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS] + "_prevalenceYearly_G.csv");
            PrintWriter pri_G = new PrintWriter(prevalenceYearly_G);

            File prevalenceYearly_A = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS] + "_prevalenceYearly_A.csv");
            PrintWriter pri_A = new PrintWriter(prevalenceYearly_A);

            File prevalenceYearly_R = new File(baseDir, FILE_NAMES_OBJ[FILE_SNAPCOUNTS] + "_prevalenceYearly_R.csv");
            PrintWriter pri_R = new PrintWriter(prevalenceYearly_R);

            PrintWriter[] priEnt = new PrintWriter[]{pri_Any, pri_G, pri_A, pri_R};

            for (int site = 0; site < priEnt.length; site++) {

                if (yearly_prevalenceBySiteAndStrainAt.length > 0) {
                    StringBuilder header = new StringBuilder();
                    header.append("Sim");
                    for (int timePt = 0; timePt < yearly_prevalenceBySiteAndStrainAt[0].length; timePt++) {
                        header.append(',');
                        header.append("Yr_");
                        header.append(timePt);

                        for (int s = 0; s < yearly_prevalenceBySiteAndStrainAt[0][timePt].length - 2; s++) {
                            header.append(',');
                        }

                    }

                    priEnt[site].println(header.toString());
                    priEnt[site].flush();
                }
                for (int s = 0; s < yearly_prevalenceBySiteAndStrainAt.length; s++) {
                    priEnt[site].print(s);
                    for (int timePt = 0; timePt < yearly_prevalenceBySiteAndStrainAt[s].length; timePt++) {
                        for (int strainType = 0; strainType < yearly_prevalenceBySiteAndStrainAt[s][timePt][site].length; strainType++) {
                            priEnt[site].print(',');
                            priEnt[site].print(yearly_prevalenceBySiteAndStrainAt[s][timePt][site][strainType]);
                        }
                    }
                    priEnt[site].println();
                }

                priEnt[site].close();
            }

        } else {
            System.err.println("Simulation_MSM_Population.decodeSnapCountFile: " + snapCountFile.getAbsolutePath() + " doesn't exist");
        }

    }

    private static class Callable_decodeExportPop implements Callable<int[][]> {

        File popFile, exportDir;
        Pattern p;

        static final String inf_stat_header = "Sim,Total,Any,G,A,R,Any_10,G01,A01,R01,G10,A10,R10,G11,A11,R11";

        static final int OFFSET_ANY = 2;
        static final int OFFSET_SITE = 3;
        static final int OFFSET_ANY_10 = 6;
        static final int OFFSET_S01 = 7;
        static final int OFFSET_S10 = OFFSET_S01 + 3;
        static final int OFFSET_S11 = OFFSET_S10 + 3;

        Callable_decodeExportPop(File popFile, File exportDir, Pattern p) {
            this.popFile = popFile;
            this.exportDir = exportDir;
            this.p = p;

        }

        @Override
        public int[][] call() throws Exception {
            int[][] res = new int[5][]; // 0 = count, 1 = map_NumberCasual6Months, 2 = map_NumberCasual6MonthInfected, 3 = map_NumberCasual6MonthInfectedNewStrain
            // 5 = newStrainHasRegPartner

            Matcher m = p.matcher(popFile.getName());
            m.find();
            int popId = new Integer(m.group(1));

            int[] map_NumberCasual6Months = new int[20];
            int[] map_NumberCasual6MonthInfected = new int[20];
            int[] map_NumberCasual6MonthInfectedNewStrain = new int[20];
            int[] count = new int[inf_stat_header.split(",").length];

            ArrayList<Integer> newStrainHasReg = new ArrayList<>();

            count[0] = popId;

            // Check for decoded indivdual snapshot             
            File decodedSnapFile = new File(exportDir, SinglePopRunnable.EXPORT_INDIV_PREFIX + popId + ".csv");

            if (decodedSnapFile.exists()) {
                BufferedReader csv = new BufferedReader(new FileReader(decodedSnapFile));
                String line;

                while ((line = csv.readLine()) != null) {

                    // Id,Age,BehavType,# Reg,# Cas,# Cas in 6 month,Inf Stat_0,Strain Stat_0,Inf Stat_1,Strain Stat_1,Inf Stat_2,Strain Stat_2
                    if (!line.startsWith("Id")) {
                        String[] ent = line.split(",");
                        int INF_OFFSET = 6;
                        int numCasual6Months = Integer.parseInt(ent[5]);
                        int[] infStat = new int[(ent.length - INF_OFFSET) / 2];
                        int[] strainStat = new int[infStat.length];
                        int numReg = Integer.parseInt(ent[3]);

                        if (numCasual6Months >= map_NumberCasual6Months.length) {
                            map_NumberCasual6Months = Arrays.copyOf(map_NumberCasual6Months, numCasual6Months + 1);
                            map_NumberCasual6MonthInfected = Arrays.copyOf(map_NumberCasual6MonthInfected, numCasual6Months + 1);
                            map_NumberCasual6MonthInfectedNewStrain = Arrays.copyOf(map_NumberCasual6MonthInfectedNewStrain, numCasual6Months + 1);
                        }

                        for (int p = 0; p < infStat.length; p++) {
                            infStat[p] = Integer.parseInt(ent[INF_OFFSET + 2 * p]);
                            strainStat[p] = Integer.parseInt(ent[INF_OFFSET + 2 * p + 1]);
                        }

                        updateExportedPopCount(count, infStat, strainStat, numReg, numCasual6Months,
                                map_NumberCasual6Months, map_NumberCasual6MonthInfected, map_NumberCasual6MonthInfectedNewStrain, newStrainHasReg);

                    }

                }

            } else { // Do it directly from pop file

                File temp = FileZipper.unzipFile(popFile, exportDir);
                MSMPopulation pop;
                try (ObjectInputStream inStream = new ObjectInputStream(new FileInputStream(temp))) {
                    pop = MSMPopulation.importMSMPopulation(inStream);
                }
                temp.delete();

                if (pop != null) {
                    count[1] = pop.getPop().length;

                    for (AbstractIndividualInterface p : pop.getPop()) {

                        int[] casualRec = ((RelationshipPerson_MSM) p).getCasualRecord();
                        int numCasual = 0;

                        for (int i = 0; i < casualRec.length; i++) {
                            numCasual += (casualRec[i] != 0) ? 1 : 0;
                        }

                        if (numCasual >= map_NumberCasual6Months.length) {
                            map_NumberCasual6Months = Arrays.copyOf(map_NumberCasual6Months, numCasual + 1);
                            map_NumberCasual6MonthInfected = Arrays.copyOf(map_NumberCasual6MonthInfected, numCasual + 1);
                            map_NumberCasual6MonthInfectedNewStrain = Arrays.copyOf(map_NumberCasual6MonthInfectedNewStrain, numCasual + 1);
                        }

                        int[] infStat = p.getInfectionStatus();
                        int[] strainStat = null;
                        if (p instanceof MultiSiteMultiStrainPersonInterface) {
                            strainStat = ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite();
                        }

                        int numReg = pop.getRelMap()[MSMPopulation.MAPPING_REG].containsVertex(p.getId())
                                ? pop.getRelMap()[MSMPopulation.MAPPING_REG].edgesOf(p.getId()).size() : 0;

                        updateExportedPopCount(count, infStat, strainStat, numReg, numCasual,
                                map_NumberCasual6Months, map_NumberCasual6MonthInfected, map_NumberCasual6MonthInfectedNewStrain,
                                newStrainHasReg);
                    }

                }
            }

            res[0] = count;
            res[1] = map_NumberCasual6Months;
            res[2] = map_NumberCasual6MonthInfected;
            res[3] = map_NumberCasual6MonthInfectedNewStrain;
            res[4] = new int[newStrainHasReg.size()];
            for (int i = 0; i < res[4].length; i++) {
                res[4][i] = newStrainHasReg.get(i);
            }

            System.out.println("Analysing pop file " + popFile.getName() + " completed.");
            return res;
        }

        protected void updateExportedPopCount(int[] count, int[] infStat, int[] strainStat, int numReg, int numCasual,
                int[] map_NumberCasual6Months, int[] map_NumberCasual6MonthInfected, int[] map_NumberCasual6MonthInfectedNewStrain,
                ArrayList<Integer> newStrainHasReg) {
            // Size of population
            count[1]++;

            // Infection stat
            boolean hasInf = false;
            boolean hasNewStrain = false;
            for (int site = 0; site < infStat.length; site++) {
                if (infStat[site] != AbstractIndividualInterface.INFECT_S) {
                    count[OFFSET_SITE + site]++;
                    hasInf = true;
                    if (strainStat != null && strainStat[site] > 0) {
                        hasNewStrain |= (strainStat[site] == 0b10 || strainStat[site] == 0b11);
                    }
                }
                if (strainStat != null && strainStat[site] > 0) {
                    switch (strainStat[site]) {
                        case 0b01:
                            count[OFFSET_S01 + site]++;
                            break;
                        case 0b10:
                            count[OFFSET_S10 + site]++;
                            break;
                        case 0b11:
                            count[OFFSET_S11 + site]++;
                            break;
                        default:
                            System.err.println("Simulation_MSM_Population.decodeExportedPopulationFiles(): strain stat 0b"
                                    + Integer.toBinaryString(strainStat[site]) + " not supported");
                    }
                }
            }

            map_NumberCasual6Months[numCasual]++;
            if (hasInf) {
                count[OFFSET_ANY]++;
                map_NumberCasual6MonthInfected[numCasual]++;
            }

            if (hasNewStrain) {
                count[OFFSET_ANY_10]++;
                map_NumberCasual6MonthInfectedNewStrain[numCasual]++;
                newStrainHasReg.add(numReg);
            }

        }

    }

    public static void decodeExportedPopulationFiles(File baseDir) throws IOException, ClassNotFoundException {
        decodeExportedPopulationFiles(baseDir, new int[]{0, Integer.MAX_VALUE});
    }

    public static void decodeExportedPopulationFiles(File baseDir, int[] range) throws IOException, ClassNotFoundException {

        File[] exportFolders = baseDir.listFiles(new FileFilter() {
            @Override
            public boolean accept(File file) {
                return file.isDirectory() && file.getName().startsWith(EXPORT_PREFIX);
            }
        });

        for (File exportDir : exportFolders) {

            File[] popFiles;

            popFiles = exportDir.listFiles(new FileFilter() {
                @Override
                public boolean accept(File file) {
                    Matcher m = Pattern_popStat.matcher(file.getName());
                    return m.find();
                }
            });

            Arrays.sort(popFiles, new Comparator<File>() {
                @Override
                public int compare(File f1, File f2) {
                    Matcher m1 = Pattern_popStat.matcher(f1.getName());
                    Matcher m2 = Pattern_popStat.matcher(f2.getName());
                    m1.find();
                    m2.find();
                    Integer n1 = new Integer(m1.group(1));
                    Integer n2 = new Integer(m2.group(1));
                    return n1.compareTo(n2);
                }
            });

            if (popFiles.length == 0) {

                System.out.println("Pop snap shot files not found. Retry using popFile instead.");

                popFiles = exportDir.listFiles(new FileFilter() {
                    @Override
                    public boolean accept(File file) {
                        Matcher m = Pattern_importFile.matcher(file.getName());
                        return m.find();
                    }
                });

                Arrays.sort(popFiles, new Comparator<File>() {
                    @Override
                    public int compare(File f1, File f2) {
                        Matcher m1 = Pattern_importFile.matcher(f1.getName());
                        Matcher m2 = Pattern_importFile.matcher(f2.getName());
                        m1.find();
                        m2.find();
                        Integer n1 = new Integer(m1.group(1));
                        Integer n2 = new Integer(m2.group(1));
                        return n1.compareTo(n2);
                    }
                });

            }

            System.out.println("Analysing " + popFiles.length + " pop file"
                    + (popFiles.length == 1 ? "" : "s") + " at " + exportDir.getAbsolutePath());

            boolean[] prePrintExist = new boolean[FILE_NAMES_CSV.length];

            for (int i = 0; i < prePrintExist.length; i++) {
                prePrintExist[i] = new File(exportDir, FILE_NAMES_CSV[i]).exists();
            }

            PrintWriter pri_inf_stat = new PrintWriter(
                    new FileWriter(new File(exportDir, FILE_NAMES_CSV[FILE_INFECTION_STAT_CSV]), true));

            if (!prePrintExist[FILE_INFECTION_STAT_CSV]) {
                pri_inf_stat.println(Callable_decodeExportPop.inf_stat_header);
            }

            PrintWriter pri_numPartnerLast6Months = new PrintWriter(
                    new FileWriter(new File(exportDir, FILE_NAMES_CSV[FILE_NUM_PARTERS_IN_LAST_6_MONTHS]), true));

            if (!prePrintExist[FILE_NUM_PARTERS_IN_LAST_6_MONTHS]) {
                pri_numPartnerLast6Months.println("Sim, Num casual partners in last 6 months, Freq, Freq (infected), Freq (new strain)");
            }

            PrintWriter pri_newStrainHasRegPartner = new PrintWriter(
                    new FileWriter(new File(exportDir, FILE_NAMES_CSV[FILE_NEW_STRAIN_HAS_REG_PARTNERS]), true));
            if (!prePrintExist[FILE_NEW_STRAIN_HAS_REG_PARTNERS]) {
                pri_newStrainHasRegPartner.println("Sim, Num of Reg for those with new strain");
            }

            PrintWriter[] writers = new PrintWriter[]{pri_inf_stat, pri_numPartnerLast6Months, pri_newStrainHasRegPartner};

            ExecutorService threadpool = null;
            int inThreadPool = 0;
            java.util.concurrent.Future<int[][]>[] result_Map = new java.util.concurrent.Future[popFiles.length];
            int exportedSoFar = 0;
            int THREAD_POOLSIZE = Runtime.getRuntime().availableProcessors();

            for (File popFile : popFiles) {
                Pattern p;
                if (Pattern_popStat.matcher(popFile.getName()).matches()) {
                    p = Pattern_popStat;
                } else {
                    p = Pattern_importFile;
                }
                Matcher m = p.matcher(popFile.getName());
                m.find();
                int popId = new Integer(m.group(1));

                if (!prePrintExist[FILE_INFECTION_STAT_CSV]
                        || popId >= range[0] && popId < range[1]) {

                    System.out.println("Submitting thread to analysis pop file " + popFile.getName());

                    if (threadpool == null) {
                        threadpool = Executors.newFixedThreadPool(THREAD_POOLSIZE);
                        inThreadPool = 0;
                    }

                    Callable_decodeExportPop thread = new Callable_decodeExportPop(popFile, exportDir, p);
                    result_Map[popId] = threadpool.submit(thread);
                    inThreadPool++;

                    if (inThreadPool == THREAD_POOLSIZE) {
                        exportedSoFar = exeutePopDecodeThread(threadpool, exportedSoFar,
                                result_Map, writers);

                        threadpool = null;
                        inThreadPool = 0;
                    }

                } else {
                    System.out.println("Analysing pop file " + popFile.getName() + " skipped.");
                }
            }
            if (inThreadPool > 0) {
                exportedSoFar = exeutePopDecodeThread(threadpool, exportedSoFar,
                        result_Map, writers);
                threadpool = null;
                inThreadPool = 0;
            }

            for (PrintWriter writer : writers) {
                writer.close();
            }

        }

    }

    private static int exeutePopDecodeThread(ExecutorService threadpool,
            int exportedSoFar, Future<int[][]>[] result_Map,
            PrintWriter[] writers) {

        PrintWriter pri_inf_stat, pri_numPartnerLast6Months, pri_newStrainHasRegPartner;

        pri_inf_stat = writers[0];
        pri_numPartnerLast6Months = writers[1];
        pri_newStrainHasRegPartner = writers[2];

        threadpool.shutdown();

        try {
            if (!threadpool.awaitTermination(1, TimeUnit.DAYS)) {
                System.err.println("Thread time-out in decoding exported pop!");
            }
        } catch (InterruptedException ex) {
            StringWriter str = new StringWriter();
            try (PrintWriter wri = new PrintWriter(str)) {
                ex.printStackTrace(wri);
            }
            System.err.println(str.toString());
        }

        while (exportedSoFar < result_Map.length && result_Map[exportedSoFar] != null) {
            java.util.concurrent.Future<int[][]> resFuture = result_Map[exportedSoFar];
            try {
                int[][] res = resFuture.get();
                int[] count = res[0];
                int[] map_NumberCasual6Months = res[1];
                int[] map_NumberCasual6MonthInfected = res[2];
                int[] map_NumberCasual6MonthInfectedNewStrain = res[3];
                int[] newStrainHasReg = res[4];

                if (newStrainHasReg.length == 0) {
                    newStrainHasReg = new int[]{-1};
                }

                pri_inf_stat.print(exportedSoFar);
                for (int c = 1; c < count.length; c++) {
                    pri_inf_stat.print(',');
                    pri_inf_stat.print(count[c]);

                }
                pri_inf_stat.println();
                pri_inf_stat.flush();

                for (int c = 1; c < map_NumberCasual6Months.length; c++) {
                    pri_numPartnerLast6Months.print(exportedSoFar);
                    pri_numPartnerLast6Months.print(',');
                    pri_numPartnerLast6Months.print(c);
                    pri_numPartnerLast6Months.print(',');
                    pri_numPartnerLast6Months.print(map_NumberCasual6Months[c]);
                    pri_numPartnerLast6Months.print(',');
                    pri_numPartnerLast6Months.print(map_NumberCasual6MonthInfected[c]);
                    pri_numPartnerLast6Months.print(',');
                    pri_numPartnerLast6Months.println(map_NumberCasual6MonthInfectedNewStrain[c]);
                }

                pri_numPartnerLast6Months.flush();

                pri_newStrainHasRegPartner.print(exportedSoFar);
                for (int k = 0; k < newStrainHasReg.length; k++) {
                    pri_newStrainHasRegPartner.print(',');
                    pri_newStrainHasRegPartner.print(newStrainHasReg[k]);
                }
                pri_newStrainHasRegPartner.println();
                pri_newStrainHasRegPartner.flush();

            } catch (InterruptedException | ExecutionException ex) {
                StringWriter str = new StringWriter();
                try (PrintWriter wri = new PrintWriter(str)) {
                    ex.printStackTrace(wri);
                }
                System.err.println(str.toString());

            }

            exportedSoFar++;
        }

        return exportedSoFar;
    }

    protected boolean populationImport(SinglePopRunnable[] runnable, int r) {
        File importResFile;
        boolean useImport = false;
        if (preExtractField != null) {
            if (runnable[r].getId() < preExtractField.length
                    && preExtractField[runnable[r].getId()] != null) {
                Object[] fields = preExtractField[runnable[r].getId()];
                for (int f = 0; f < fields.length; f++) {
                    Object ent = fields[f];
                    if (propModelInit[f] != null && !propModelInit[f].isEmpty()) {
                        ent = StaticMethods.propStrToObject(propModelInit[f],
                                ((AbstractRegCasRelMapPopulation) runnable[r].getPopulation()).getFieldClass(f));

                    }
                    runnable[r].getPopulation().setParameter(" ", f, ent);
                }

                runnable[r].setModelBurnIn(0);
                runnable[r].getPopulation().initialiseInfection(0);
                useImport = true;
            } else {
                System.out.println("Pre extracted import pop for Thread #" + runnable[r].getId()
                        + " not found. A new pop with default burn in will be generated");
            }

        } else {
            File importResDir = new File((String) propVal[PROP_POP_IMPORT_PATH]);
            importResFile = new File(importResDir, SimulationInterface.POP_FILE_PREFIX + runnable[r].getId() + ".zip");

            if (importResFile.isFile()) {
                useImport = true;
                runnable[r].setImportFile(importResFile, propModelInit);
                runnable[r].getPopulation().initialiseInfection(1); // Dummy
            } else {
                System.out.println("Import pop at '" + importResFile.getAbsolutePath()
                        + "' not found. A new pop with default burn in will be generated");
            }
        }
        return useImport;
    }

    public Object[] getPropVal() {
        return propVal;
    }

    public String[] getPropModelInit() {
        return propModelInit;

    }

    private final class CLASSIFIER_PREVAL implements PersonClassifier {

        int siteId;

        public CLASSIFIER_PREVAL(int siteId) {
            this.siteId = siteId;
        }

        @Override
        public int classifyPerson(AbstractIndividualInterface p) {

            if (p instanceof MultiSiteMultiStrainPersonInterface) {
                if (siteId >= 0) {
                    return ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite()[siteId];
                } else {
                    int combineSite = 0;
                    int[] strainStat = ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite();
                    for (int i = 0; i < strainStat.length; i++) {
                        combineSite = combineSite | ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite()[i];
                    }
                    return combineSite;
                }
            } else {
                return -1;
            }
        }

        @Override
        public int numClass() {
            return 4; // 0 - not infected, 0b01, 0b10, 0b11 
        }
    }

}
