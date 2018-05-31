package sim;

import java.beans.PropertyChangeSupport;
import java.io.EOFException;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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
 * @version 20180531
 *
 * <pre>
 * History
 * 20150528
 *     - Change the role and renaming of snap count
 * 20180530
 *     - Add support for skip thread based on snap count
 * 20180531
 *     - Adjustment snapshot count sampling frequency
 * </pre>
 */
public class Simulation_MSM_Population implements SimulationInterface {

    public static final String[] PROP_NAME_MSM_MIS = {
        "PROP_STRAINS_INTRO_AT", "PROP_STRAINS_COEXIST_MAT"
    };

    public static final Class[] PROP_CLASS_MSM_MIS = {
        float[][].class, // float[]{globaltime, strainNum, site, likelihood to co-exist, number of infection to introduce, frequency (optional) }
        float[][].class, // float[exist_strain][new_strain]{likelihood to coexist}
    };

    public static final int PROP_STRAINS_INTRO_AT = PROP_NAME.length;
    public static final int PROP_STRAINS_COEXIST_MAT = PROP_STRAINS_INTRO_AT + 1;

    // Output filename
    public static final String[] FILE_NAMES_OBJ = {"endNumInf.obj", "extinctAt.obj", "snapCount.obj",
        "eventPt.obj", "incidenceCount.obj"};
    public static final int FILE_END_NUM_INF = 0;
    public static final int FILE_EXTINCT_AT = FILE_END_NUM_INF + 1;
    public static final int FILE_SNAPCOUNTS = FILE_EXTINCT_AT + 1;
    public static final int FILE_EVENT_POINTER = FILE_SNAPCOUNTS + 1;
    public static final int FILE_INCIDENT_COUNT = FILE_EVENT_POINTER + 1;

    public static final String[] FILE_NAMES_CSV = {"endNumInfPerson.csv", "infStatSnapshot.csv"};
    public static final int FILE_END_NUM_INF_PERSON_CSV = 0;
    public static final int FILE_INFECTION_STAT_CSV = FILE_END_NUM_INF_PERSON_CSV + 1;

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
    // Skipping seed
    private boolean[] skipNThSeed = new boolean[0];

    private boolean useImportIOThread = true;

    public boolean[] getSkipNThSeed() {
        return skipNThSeed;
    }

    public void setSkipNThSeed(boolean[] skipNThSeed) {
        this.skipNThSeed = skipNThSeed;
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

            } else {

                ExecutorService executor;
                SinglePopRunnable[] runnable = null;

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
                
                executor = Executors.newFixedThreadPool(numThreads);

                for (int r = 0; r < runnable.length; r++) {

                    runnable[r] = new SinglePopRunnable(threadCounter,
                            ((Number) propVal[PROP_NUM_SNAP]).intValue(), ((Number) propVal[PROP_SNAP_FREQ]).intValue());

                    runnable[r].setBaseDir(baseDir);

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

                    runnable[r].setPrintPrevalence(((Integer) propVal[PROP_USE_PARALLEL]) == 0);

                    threadCounter++;

                    if (progressSupport != null) {
                        runnable[r].setProgressSupport(progressSupport);
                    }

                    while (rngCallCounter < skipNThSeed.length && skipNThSeed[rngCallCounter]) {
                        System.out.println("Thread #" + (threadCounter - 1) + " skip seed of " + rng.nextLong());
                        rngCallCounter++;
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
            }

        }

        finalise(simSoFar);

    }

    private boolean checkForNullArray(Object arr) {
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

        PrintWriter pri = new PrintWriter(new FileWriter(new File(baseDir, FILE_NAMES_CSV[FILE_END_NUM_INF_PERSON_CSV])));

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

                int[] sel = new int[]{-1, 0, 1, 2};
                int[] numInfected = runnable.getPopulation().getNumberOfInfected(sel);

                pri.print(runnable.getId());
                pri.print(',');
                pri.print(runnable.getPopulation().getPop().length);
                for (int i = 0; i < sel.length; i++) {
                    pri.print(',');
                    pri.print(numInfected[i]);
                }
                pri.println();
            }

        }
        for (ObjectOutputStream objS1 : objS) {
            if (objS1 != null) {
                objS1.close();
            }
        }
        pri.close();
    }

    public void decodeSnapCountFile() throws FileNotFoundException, IOException, ClassNotFoundException {
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

            //Format: from snapshotCountClassifier
            final int PREVAL_ALL = 0;
            final int PREVAL_G = 0;
            final int PREVAL_A = 0;
            final int PREVAL_R = 0;
                                   
            
            
            

            for (int[][][] entry : entCollection) {
                int importAtSnap = -1;
                
                int snapFreq = (int) propVal[PROP_SNAP_FREQ];
                final int snapPerYr = 360/snapFreq;
                
                
                if(propVal[PROP_STRAINS_INTRO_AT] != null){
                    float[] inFirstStrain  = ((float[][]) propVal[PROP_STRAINS_INTRO_AT])[0];                    
                    int importTime = (int) inFirstStrain[0];       
                    int burnIn = propVal[PROP_BURNIN] == null?  0: (int) propVal[PROP_BURNIN];                    
                    importAtSnap = ((importTime - burnIn) / snapFreq) -1;
                }
                

                for (int snapCounter = 0; snapCounter < entry.length; snapCounter++) {
                    if (entry[snapCounter][PREVAL_ALL][2] + entry[snapCounter][PREVAL_ALL][3] > 0) {
                        newStrainExinctAtSnapshot[simCounter]++;                        
                    }
                    if (importAtSnap > 0) {
                        int timePt = -1;
                        if(snapCounter - importAtSnap == snapPerYr){
                            timePt = 0;
                        }else if (snapCounter - importAtSnap == 2 * snapPerYr){
                            timePt = 1;
                        }else if (snapCounter - importAtSnap == 5 * snapPerYr){
                             timePt = 2;
                        }else if (snapCounter - importAtSnap == 10 * snapPerYr){
                             timePt = 3;
                        }                                                
                       
                        if (timePt >= 0) {
                            prevalenceByStrainAt[simCounter][timePt][0] = entry[snapCounter][PREVAL_ALL][0];
                            prevalenceByStrainAt[simCounter][timePt][1] = entry[snapCounter][PREVAL_ALL][1];
                            prevalenceByStrainAt[simCounter][timePt][2] = entry[snapCounter][PREVAL_ALL][2] + entry[snapCounter][PREVAL_ALL][3];
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

        } else {
            System.err.println(getClass().getName() + "decodeSnapCountFile: " + snapCountFile.getAbsolutePath() + " doesn't exist");
        }

    }

    public void decodeExportedPopulationFiles() throws IOException, ClassNotFoundException {
        final Pattern Pattern_importFile = Pattern.compile("pop_(\\d+).zip");
        File[] exportFolders = baseDir.listFiles(new FileFilter() {
            @Override
            public boolean accept(File file) {
                return file.isDirectory() && file.getName().startsWith(EXPORT_PREFIX);
            }
        });

        for (File exportDir : exportFolders) {
            File[] popZip = exportDir.listFiles(new FileFilter() {
                @Override
                public boolean accept(File file) {
                    Matcher m = Pattern_importFile.matcher(file.getName());
                    return m.find();
                }
            });

            Arrays.sort(popZip, new Comparator<File>() {
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

            System.out.println("Analysing " + popZip.length + " pop file"
                    + (popZip.length == 1 ? "" : "s") + " at " + exportDir.getAbsolutePath());

            final int OFFSET_ANY = 2;
            final int OFFSET_SITE = 3;
            final int OFFSET_ANY_10 = 6;
            final int OFFSET_S01 = 7;
            final int OFFSET_S10 = OFFSET_S01 + 3;
            final int OFFSET_S11 = OFFSET_S10 + 3;

            try (PrintWriter preOutputFile = new PrintWriter(new File(exportDir, FILE_NAMES_CSV[FILE_INFECTION_STAT_CSV]))) {
                String header = "Sim,Total,Any,G,A,R,Any_10,G01,A01,R01,G10,A10,R10,G11,A11,R11";
                preOutputFile.println(header);

                for (File popFile : popZip) {
                    Matcher m = Pattern_importFile.matcher(popFile.getName());
                    m.find();
                    int popId = new Integer(m.group(1));
                    int[] count = new int[header.split(",").length];
                    count[0] = popId;

                    File temp = FileZipper.unzipFile(popFile, exportDir);
                    MSMPopulation pop;
                    try (ObjectInputStream inStream = new ObjectInputStream(new FileInputStream(temp))) {
                        pop = MSMPopulation.importMSMPopulation(inStream);
                    }
                    temp.delete();

                    if (pop != null) {
                        count[1] = pop.getPop().length;
                        for (AbstractIndividualInterface p : pop.getPop()) {
                            int[] infStat = p.getInfectionStatus();
                            int[] strainStat = null;
                            if (p instanceof MultiSiteMultiStrainPersonInterface) {
                                strainStat = ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite();
                            }

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
                                            System.err.println(getClass().getName() + ".decodeExportedPopulationFiles(): strain stat 0b"
                                                    + Integer.toBinaryString(strainStat[site]) + " not supported");
                                    }
                                }
                            }
                            if (hasInf) {
                                count[OFFSET_ANY]++;
                            }

                            if (hasNewStrain) {
                                count[OFFSET_ANY_10]++;
                            }
                        }
                    }

                    preOutputFile.print(popId);
                    for (int c = 1; c < count.length; c++) {
                        preOutputFile.print(',');
                        preOutputFile.print(count[c]);

                    }
                    preOutputFile.println();
                    preOutputFile.flush();

                    System.out.println("Analysing pop file " + popFile.getName() + " completed.");
                }
            }

        }

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
