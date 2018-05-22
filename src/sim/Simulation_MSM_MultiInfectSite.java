package sim;

import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.commons.math3.distribution.PoissonDistribution;
import person.AbstractIndividualInterface;
import population.AbstractRegCasRelMapPopulation;
import population.MSMPopulation_MultiType_Strains;
import population.person.RelationshipPerson;
import population.person.RelationshipPerson_MSM;
import random.MersenneTwisterFastRandomGenerator;
import random.RandomGenerator;
import util.AppendableObjOutstream;
import util.FileZipper;
import util.PersonClassifier;
import util.StaticMethods;

/**
 *
 * @author Ben Hui
 */
public class Simulation_MSM_MultiInfectSite implements SimulationInterface {

    public static final String[] PROP_NAME_MSM_MIS = {
        "PROP_INIT_RESIST_PERCENT", "PROP_RESIST_STRAINS_INTRO_AT",};

    public static final Class[] PROP_CLASS_MSM_MIS = {
        float[][][].class, float[][].class, // float[infNum]{m, k} = m of them for strain i at k-th snapshot or float[infNum]{m, k, freq}
    };

    public static final int PROP_INIT_RESIST_PERCENT = PROP_NAME.length;
    public static final int PROP_RESIST_STRAINS_INTRO_AT = PROP_INIT_RESIST_PERCENT + 1;

    // Output filename
    public static final String[] FILE_NAMES = {"endNumInf.obj", "extinctAt.obj", "outputCount.obj",
        "eventPt.obj", "popSnapCount.obj", "strainIntroAtSnap.csv", "strainIntroCountAtSnap.csv"};
    public static final int FILE_END_NUM_INF = 0;
    public static final int FILE_EXTINCT_AT = FILE_END_NUM_INF + 1;
    public static final int FILE_SNAPCOUNTS = FILE_EXTINCT_AT + 1;
    public static final int FILE_EVENT_POINTER = FILE_SNAPCOUNTS + 1;
    public static final int FILE_POP_SNAPCOUNT = FILE_EVENT_POINTER + 1;
    public static final int FILE_STRAIN_INTRO_TIME_CSV = FILE_POP_SNAPCOUNT + 1;
    public static final int FILE_STRAIN_INTRO_COUNT_CSV = FILE_STRAIN_INTRO_TIME_CSV + 1;

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
                propVal[i] = StaticMethods.propStrToObject(ent, PROP_CLASS_MSM_MIS[i]);
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
        if (!MSMPopulation_MultiType_Strains.class.getName().equals(propVal[PROP_POP_TYPE])) {
            throw new UnsupportedOperationException("Class " + propVal[PROP_POP_TYPE]
                    + " not supported in " + getClass().getName() + ".");
        }

        // Snapshot classifier
        PersonClassifier[] snapshotCountClassifer = new PersonClassifier[]{
            new STRAIN_SITE_COUNTER(RelationshipPerson_MSM.SITE_G, 2),
            new STRAIN_SITE_COUNTER(RelationshipPerson_MSM.SITE_A, 2),
            new STRAIN_SITE_COUNTER(RelationshipPerson_MSM.SITE_R, 2)
        };
        boolean[] snapshotAccum = new boolean[snapshotCountClassifer.length]; // All non-accumulative
        Arrays.fill(snapshotAccum, false);

        setSnapshotSetting(snapshotCountClassifer, snapshotAccum);
        // Initial prevalence - seeding
        float[][] init_preval = new float[][]{{100f}, {460f}, {400f}}; // init_preval[site][class=0] = number of infected
        int[] init_expose = new int[]{3, 12 * 30, 12 * 7};             // init_expose[site]
        float[][][] init_strain_decomp = new float[][][]{ // init_strain_decomp[site][class=0][strain decompositon]
            {{0, 100f}}, {{0, 460f}}, {{0, 400f}}};                    // All with strain type of 0b01 
        PersonClassifier[] init_preval_classifier = new PersonClassifier[]{
            new ANY_PERSON(RelationshipPerson_MSM.SITE_G),
            new ANY_PERSON(RelationshipPerson_MSM.SITE_A),
            new ANY_PERSON(RelationshipPerson_MSM.SITE_R),};

        if (propVal[PROP_INIT_RESIST_PERCENT] != null) {
            init_strain_decomp = (float[][][]) propVal[PROP_INIT_RESIST_PERCENT];
        }

        // Simulation setting        
        int numSimTotal = ((Number) propVal[PROP_NUM_SIM_PER_SET]).intValue();
        boolean useParallel = ((Integer) propVal[PROP_USE_PARALLEL]) != 0;
        int simSoFar = 0;
        int threadCounter = 0;
        MersenneTwisterFastRandomGenerator rng = new random.MersenneTwisterFastRandomGenerator(((Number) propVal[PROP_BASESEED]).longValue());
        int numProcess = Math.min(getMaxThreads(), Runtime.getRuntime().availableProcessors());

        if (Integer.parseInt(propVal[PROP_USE_PARALLEL].toString()) > 0) {
            numProcess = Math.min(Integer.parseInt(propVal[PROP_USE_PARALLEL].toString()),
                    numProcess);
        }

        int[][] eventPointers = null; // Not used

        while (simSoFar < numSimTotal && !stopNextTurn) {

            int numThreads = Math.min(numProcess, numSimTotal - simSoFar);
            // Generate strain intro file file
            PrintWriter strainIntroTimeWri = null;
            PrintWriter strainIntroCountWri = null;

            if (propVal[PROP_RESIST_STRAINS_INTRO_AT] != null) {
                // Generate strain intro file file                 
                strainIntroTimeWri = new PrintWriter(new FileWriter(new File(baseDir, FILE_NAMES[FILE_STRAIN_INTRO_TIME_CSV]), true));
                strainIntroCountWri = new PrintWriter(new FileWriter(new File(baseDir, FILE_NAMES[FILE_STRAIN_INTRO_COUNT_CSV]), true));
            }

            if (useParallel) {
                showStrStatus("Running S" + threadCounter + " to S" + (threadCounter + numThreads - 1) + " in parallel ...");
            }

            // Creating threads
            SinglePopRunnable[] runnable;
            runnable = new SinglePopRunnable[numThreads];

            for (int threadIndex = 0; threadIndex < runnable.length; threadIndex++) {
                runnable[threadIndex] = new SinglePopRunnable(threadCounter,
                        ((Number) propVal[PROP_NUM_SNAP]).intValue(), ((Number) propVal[PROP_SNAP_FREQ]).intValue());

                runnable[threadIndex].setBaseDir(baseDir);

                if (propVal[PROP_INF_HIST_PREFIX] != null) {
                    runnable[threadIndex].setInfectioHistoryPrefix((String) propVal[PROP_INF_HIST_PREFIX]);
                }

                threadCounter++;

                if (progressSupport != null) {
                    runnable[threadIndex].setProgressSupport(progressSupport);
                }

                runnable[threadIndex].setPopulation(new MSMPopulation_MultiType_Strains(rng.nextLong()));

                boolean useImport = false;
                if (propVal[PROP_POP_IMPORT_PATH] != null) {
                    useImport = populationImport(runnable, threadIndex);
                }

                if (!useImport) {
                    System.out.println("Thread #" + (threadCounter - 1) + " generated with seed of " + runnable[threadIndex].getPopulation().getSeed());

                    runnable[threadIndex].model_prop_initialise(((Number) propVal[PROP_BURNIN]).intValue(), propModelInit);

                    for (int siteId = 0; siteId < init_preval_classifier.length; siteId++) {
                        runnable[threadIndex].setInfectionIntroAt(0, siteId,
                                init_preval_classifier[siteId], init_preval[siteId], init_expose[siteId],
                                init_strain_decomp == null ? null : init_strain_decomp[siteId]);
                    }
                }

                runnable[threadIndex].setSnapShotOutput(snapshotCountClassifer, snapshotCountAccum);

                // Intro strain                 
                if (propVal[PROP_RESIST_STRAINS_INTRO_AT] != null) {
                    float[][] ent = (float[][]) propVal[PROP_RESIST_STRAINS_INTRO_AT];
                    random.RandomGenerator introRNG = new random.MersenneTwisterFastRandomGenerator(runnable[threadIndex].getPopulation().getSeed());

                    int[][] preIntroAt = null; // For periodic random site allocations

                    // Random site allocations 
                    if (ent.length == 1) {
                        float[] entOrg = ent[0];
                        float[] entDist = Arrays.copyOfRange(entOrg, entOrg.length - 3, entOrg.length); // 3 sites
                        int maxDist = (int) entDist[entDist.length - 1];
                        ent = new float[entDist.length][];

                        if (entOrg[0] < 0) { // for 1/per snap freq option
                            float[][] genEnt = new float[][]{Arrays.copyOfRange(entOrg, 0, entOrg.length - 3)};
                            int[] introAtTotal = generateStrainIntroAt(genEnt, 0, introRNG);

                            preIntroAt = new int[entDist.length][introAtTotal.length];

                            // Allocate site by distribution
                            for (int i = 0; i < introAtTotal.length; i++) {
                                int count = introAtTotal[i];
                                while (count > 0) { // Allocate by site
                                    int pSite = introRNG.nextInt(maxDist);
                                    int siteId = 0;
                                    while (entDist[siteId] <= pSite && siteId < entDist.length) {
                                        siteId++;
                                    }
                                    preIntroAt[siteId][i]++;
                                    count--;
                                }
                            }

                            for (int s = 0; s < entDist.length; s++) {

                                if (preIntroAt.length > 0) {
                                    // Dummy value so it will be trigger in the 
                                    // if (ent[siteId] != null && ent[siteId].length >= 2) loop
                                    ent[s] = entOrg;
                                } else {
                                    ent[s] = null;
                                }
                            }

                        } else {
                            for (int count = 0; count < entOrg[0]; count++) {
                                int pSite = introRNG.nextInt(maxDist);
                                int siteId = 0;
                                while (entDist[siteId] <= pSite && siteId < entDist.length) {
                                    siteId++;
                                }
                                if (ent[siteId] == null) {
                                    ent[siteId] = Arrays.copyOfRange(entOrg, 0, entOrg.length - 3);
                                    ent[siteId][0] = 0;
                                }
                                ent[siteId][0]++;
                            }
                        }
                    }

                    for (int siteId = 0; siteId < ent.length; siteId++) {
                        if (ent[siteId] != null && ent[siteId].length >= 2) {
                            int[] introAt;
                            if (preIntroAt != null) {
                                introAt = preIntroAt[siteId];
                            } else {
                                introAt = generateStrainIntroAt(ent, siteId, introRNG);
                            }
                            StringBuilder strainIntroTimeStr = new StringBuilder(Integer.toString(runnable[threadIndex].getId()));
                            StringBuilder strainIntroCountStr = new StringBuilder(Integer.toString(runnable[threadIndex].getId()));

                            for (int snapNum = 0; snapNum < introAt.length; snapNum++) {
                                if (introAt[snapNum] != 0) {
                                    // Strain index 1 or change accordingly
                                    runnable[threadIndex].setStrainIntroAt(snapNum, siteId,
                                            new ANY_PERSON(false, siteId), new float[]{introAt[snapNum]}, 1);
                                    strainIntroTimeStr.append(',');
                                    strainIntroTimeStr.append(snapNum);
                                    strainIntroTimeStr.append(',');
                                    strainIntroTimeStr.append(siteId);
                                    strainIntroCountStr.append(',');
                                    strainIntroCountStr.append(introAt[snapNum]);
                                    strainIntroCountStr.append(',');
                                    strainIntroCountStr.append(siteId);
                                }
                            }

                            if (strainIntroTimeWri != null && strainIntroCountWri != null) {
                                strainIntroTimeWri.println(strainIntroTimeStr.toString());
                                strainIntroCountWri.println(strainIntroCountStr.toString());
                            }
                        }
                    }
                }

                if (eventPointers != null) {
                    eventPointers[runnable[threadIndex].getId()] = runnable[threadIndex].getEventsPointer();
                }
            }

            if (strainIntroTimeWri != null && strainIntroCountWri != null) {
                strainIntroTimeWri.close();
                strainIntroCountWri.close();
            }

            // Start threads
            ExecutorService executor = null;
            if (useParallel) {
                executor = Executors.newFixedThreadPool(numThreads);
            }

            for (SinglePopRunnable singlePopRunnable : runnable) {
                if (useParallel && executor != null) {
                    executor.submit(singlePopRunnable);
                } else {
                    showStrStatus("Running S" + singlePopRunnable.getId() + "...");
                    singlePopRunnable.run();
                }
            }

            if (useParallel && executor != null) {
                executor.shutdown();
                if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
                    showStrStatus("Thread time-out!");
                }
            }

            // ExportPop if needed 
            if (propVal[PROP_POP_EXPORT] != null && ((Boolean) propVal[PROP_POP_EXPORT])) {
                populationExport(runnable);
            }

            generateOutputFiles(eventPointers, runnable);

            if (progressSupport != null) {
                progressSupport.firePropertyChange(PROGRESS_SIM_STORED, simSoFar, simSoFar + runnable.length);
            }
            simSoFar += runnable.length;
        }

        finalise(simSoFar);

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
            StaticMethods.decodeResultObjFile(new File(baseDir, FILE_NAMES[FILE_END_NUM_INF]), simSoFar);
            StaticMethods.decodeResultObjFile(new File(baseDir, FILE_NAMES[FILE_EXTINCT_AT]), simSoFar);
        } catch (ClassNotFoundException ex) {
            System.err.println(getClass().getName() + ".generateOneResultSet: Error - corrupted data file");
        }
        if (progressSupport
                != null) {
            progressSupport.firePropertyChange(PROGRESS_ALL_DONE, null, simSoFar);
        }
    }
    
    protected void generateOutputFiles(int[][] eventPointers, SinglePopRunnable[] runnable) throws IOException {
        // Generating file
        ObjectOutputStream[] objS;

        objS = new ObjectOutputStream[FILE_NAMES.length];

        objS[FILE_END_NUM_INF] = AppendableObjOutstream.generateFromFile(new File(baseDir, FILE_NAMES[FILE_END_NUM_INF]));
        objS[FILE_EXTINCT_AT] = AppendableObjOutstream.generateFromFile(new File(baseDir, FILE_NAMES[FILE_EXTINCT_AT]));
        if (snapshotCountClassifier.length > 0) {
            objS[FILE_SNAPCOUNTS] = AppendableObjOutstream.generateFromFile(new File(baseDir, FILE_NAMES[FILE_SNAPCOUNTS]));
        }
        if (eventPointers != null) {
            objS[FILE_EVENT_POINTER] = AppendableObjOutstream.generateFromFile(new File(baseDir, FILE_NAMES[FILE_EVENT_POINTER]));
        }
        objS[FILE_POP_SNAPCOUNT] = AppendableObjOutstream.generateFromFile(new File(baseDir, FILE_NAMES[FILE_POP_SNAPCOUNT]));

        for (SinglePopRunnable runnable1 : runnable) {

            objS[FILE_END_NUM_INF].writeObject(runnable1.getPopulation().getNumInf());
            objS[FILE_END_NUM_INF].flush();
            objS[FILE_EXTINCT_AT].writeObject(runnable1.getExtinctionAt());
            objS[FILE_EXTINCT_AT].flush();
            if (objS[FILE_SNAPCOUNTS] != null) {
                objS[FILE_SNAPCOUNTS].writeObject(runnable1.getSnapCounts());
                objS[FILE_SNAPCOUNTS].flush();
            }
            if (objS[FILE_EVENT_POINTER] != null) {
                objS[FILE_EVENT_POINTER].writeObject(runnable1.getEventsPointer());
                objS[FILE_EVENT_POINTER].flush();
            }
            // PopSnapObj - Bug if there is more than one null in sim?
            int[][] popSnap = runnable1.getPopSnapCounts();
            for (int i = 0; i < popSnap.length; i++) {
                if (popSnap[i] == null) {
                    popSnap[i] = new int[0];
                }
            }

            objS[FILE_POP_SNAPCOUNT].writeObject(popSnap);
            objS[FILE_POP_SNAPCOUNT].flush();
        }
        for (ObjectOutputStream objS1 : objS) {
            if (objS1 != null) {
                objS1.close();
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
            importResDir = new File(importResDir, SimulationInterface.EXPORT_DIR);
            importResFile = new File(importResDir, SimulationInterface.EXPORT_FILE_PREFIX + runnable[r].getId() + ".zip");

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

    protected void populationExport(SinglePopRunnable[] runnable) throws IOException {
        File popDir = new File(baseDir, EXPORT_DIR);
        popDir.mkdirs();
        ObjectOutputStream objStr;
        for (SinglePopRunnable r : runnable) {
            File targetFile = new File(popDir, EXPORT_FILE_PREFIX + r.getId() + ".obj");
            File zippedPopFile = new File(popDir, EXPORT_FILE_PREFIX + r.getId() + ".zip");
            objStr = new ObjectOutputStream(new FileOutputStream(targetFile));
            r.getPopulation().exportPop(objStr);
            objStr.close();
            FileZipper.zipFile(targetFile, zippedPopFile);
            targetFile.delete();
            System.out.println("Exporting population to " + zippedPopFile.getAbsolutePath() + " completed.");
        }
    }

    private class STRAIN_SITE_COUNTER implements PersonClassifier {

        final int siteId;
        final int maxNumStrain;

        public STRAIN_SITE_COUNTER(int siteId, int maxStrain) {
            this.siteId = siteId;
            this.maxNumStrain = maxStrain;
        }

        @Override
        public int classifyPerson(AbstractIndividualInterface p) {
            if (((RelationshipPerson) p).getCurrentStrainsAtSite() != null
                    && ((RelationshipPerson) p).getCurrentStrainsAtSite()[siteId] != Double.NaN) {
                return ((RelationshipPerson) p).getCurrentStrainsAtSite()[siteId];
            } else {
                return -1;
            }
        }

        @Override
        public int numClass() {
            return (1 << maxNumStrain); // Binary
        }
    }

    protected class ANY_PERSON implements PersonClassifier {

        boolean infectedOnly = false;
        int siteId = -1;

        public ANY_PERSON(int siteId) {
            this.siteId = siteId;
        }

        public ANY_PERSON(boolean infectedOnly, int siteId) {
            this.infectedOnly = infectedOnly;
            this.siteId = siteId;
        }

        @Override
        public int classifyPerson(AbstractIndividualInterface p) {
            boolean incl;
            if (infectedOnly) {
                if (siteId < 0) {
                    incl = false;
                    for (int i = 0; i < p.getInfectionStatus().length && !incl; i++) {
                        incl |= p.getInfectionStatus()[i] != AbstractIndividualInterface.INFECT_S;
                    }
                } else {
                    incl = p.getInfectionStatus()[siteId] != AbstractIndividualInterface.INFECT_S;
                }
            } else {
                incl = true;
            }

            if (incl && p instanceof RelationshipPerson_MSM) {
                // Only those who are not immune 
                String str;
                switch (siteId) {
                    case RelationshipPerson_MSM.SITE_G:
                        str = "PARAM_IMMUNE_ACT_SITE_G";
                        break;
                    case RelationshipPerson_MSM.SITE_A:
                        str = "PARAM_IMMUNE_ACT_SITE_A";
                        break;
                    case RelationshipPerson_MSM.SITE_R:
                        str = "PARAM_IMMUNE_ACT_SITE_R";
                        break;
                    default:
                        str = null;
                }
                if (str != null) {
                    incl = (int) p.getParameter(str) == 0;
                }
            }

            return incl ? 0 : -1;

        }

        @Override
        public int numClass() {
            return 1;
        }
    }

}
