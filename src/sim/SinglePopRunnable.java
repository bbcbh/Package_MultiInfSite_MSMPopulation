package sim;

import infection.AbstractInfection;
import infection.GonorrhoeaSiteInfection;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import population.AbstractRegCasRelMapPopulation;
import population.MSMPopulation;
import person.AbstractIndividualInterface;
import relationship.RelationshipMap;

import util.FileZipper;
import util.PersonClassifier;
import util.StaticMethods;
import infection.MultiStrainInfectionInterface;
import infection.vaccination.SiteSpecificVaccination;
import java.nio.file.Files;
import population.person.MultiSiteMultiStrainPersonInterface;
import population.person.RelationshipPerson;
import relationship.SingleRelationship;

/**
 *
 * @author Ben Hui
 * @version 20190620
 *
 * History:  <pre>
 * 20150313
 *  - Add infection history support
 * 20150523
 *  - Add additional setInfectionIntroAt function for import infection from prop file
 * 20180528
 *  - Simplfy strain intro add method.
 * 20181012
 *  - Remove export population if basedir is null
 * 20190116
 *  - Add support for exporting indivdual behavior
 * 20190122
 *  - Add support for survival analysis
 * 20190123
 *  - Add relMap index for new strain spread summary CSV
 *  - Remove automatic exportAt at end of simulation run
 * 20190125
 *  - Switch export all option to person stat only by default
 *  - Add partner history stat for patient zero
 *  - Add behavior option for intro at
 * 20190127
 *  - Create a new RNG for selecting candidate within introStrain
 * 20190523
 *  - Add support for vaccine setting
 * 20190620
 *  - Add support for initial pop size
 * </pre>
 *
 */
public class SinglePopRunnable implements Runnable {

    AbstractRegCasRelMapPopulation pop;
    private final int id;
    protected final int numSnaps;
    protected final int snapFreq;
    private int lastReportedStep = 0;
    private int modelBurnIn = 0;
    private transient PropertyChangeSupport progressSupport = null;
    protected final transient InstantInfectionEvent[] infectionEvents;  // infectionEvents[snapNum]
    private final transient PersonClassifier[][] snapClassifiers;     // snapClassifiers[snapNum][classNum]
    private final transient int[][][] snapCounts;                     // snapCount[snapNum][classNum][numByclassTypeId]
    private final transient boolean[][] cummulativeSnap;              // cummulativeSnap[snapNum][classNum]
    private transient int[] extinctionAt;                             // -1 = not extinct
    private final transient int[][] incidenceCounts;                    // incidenceCounts[snapNum][countType]

    private transient int[][] strainCompositionActiveRange; // [strainsCompostion]{start_time, end_time}

    private File importFile = null;
    private String[] importSetting = null;
    private File baseDir = null;

    private String infectioHistoryPrefix = null;  // If not null, then export infection history

    private final HashMap<Integer, float[]> infectionIntroMap = new HashMap();  // Global time, prevalence
    private final HashMap<Integer, Integer> infectionPreExposeMap = new HashMap(); // infection_id, days

    public static final String EXPORT_PREFIX = "export_";
    public static final String EXPORT_FORCED_PREFIX = "forced_";
    public static final String EXPORT_INDIV_PREFIX = "IndivdualSnap_";
    private int[] exportAt = null;
    private int exportAtPt = 0;
    private int printPrevalenceAtFreq = -1;

    ArrayList<float[]> strainIntroEnt = new ArrayList<>();
    int strainIntroPt = 0;
    public static final int STRAIN_INTRO_ENTRY_INTRO_AT = 0; //Globaltime
    public static final int STRAIN_INTRO_ENTRY_STRAIN_NUM = STRAIN_INTRO_ENTRY_INTRO_AT + 1; //strainNum
    public static final int STRAIN_INTRO_ENTRY_SITE = STRAIN_INTRO_ENTRY_STRAIN_NUM + 1;// Site
    public static final int STRAIN_INTRO_ENTRY_NUM_INF = STRAIN_INTRO_ENTRY_SITE + 1; // Number of infection to introduce
    public static final int STRAIN_INTRO_ENTRY_FREQ = STRAIN_INTRO_ENTRY_NUM_INF + 1; // Frequency (optional, or set to <=0 if not needed)
    public static final int STRAIN_INTRO_ENTRY_BEHAVOR = STRAIN_INTRO_ENTRY_FREQ + 1; // Candidate behavior (optional, or set to <0 if not used)

    float[][] coexistMat;

    // For survial analysis
    boolean use_patient_zero = false;
    RelationshipPerson patient_zero = null;
    int patient_zero_start_time = -1;
    int patient_zero_end_time = -1;

    HashMap<Integer, int[]> patient_zero_partners_past_strain = null; // partner id: 
    HashMap<Integer, int[]> patient_zero_partners_record = null; // partner id: start time,end time (if ended), partnership type, partner behavor

    public static final int PATIENT_ZERO_PARNTER_REC_START_TIME = 0; // < 0  = as record started
    public static final int PATIENT_ZERO_PARNTER_REC_RELATIONSHIP_LASTED_FOR = PATIENT_ZERO_PARNTER_REC_START_TIME + 1;
    public static final int PATIENT_ZERO_PARNTER_REC_TYPE = PATIENT_ZERO_PARNTER_REC_RELATIONSHIP_LASTED_FOR + 1;
    public static final int PATIENT_ZERO_PARNTER_REC_BEHAVIOR = PATIENT_ZERO_PARNTER_REC_TYPE + 1;
    public static final int LENGTH_PATIENT_ZERO_PARNTER_REC = PATIENT_ZERO_PARNTER_REC_BEHAVIOR + 1;

    ArrayList<int[]> patient_zero_newStrainSpreadSummary = new ArrayList();

    // Vaccination
    double[][] vaccineSettingAll;
    public static final int VACCINE_SETTING_LENGTH = 0;
    public static final int VACCINE_START_TIME = VACCINE_SETTING_LENGTH + 1;
    public static final int VACCINE_END_TIME = VACCINE_START_TIME + 1;

    public void setVaccineSetting(double[][] vaccineSetting) {
        this.vaccineSettingAll = vaccineSetting;
    }

    public void set_patient_zero(boolean tsa_patient_zero) {
        this.use_patient_zero = tsa_patient_zero;
    }

    public void setCoexistMat(float[][] coexistMat) {
        this.coexistMat = coexistMat;
    }

    public void setBaseDir(File baseDir) {
        this.baseDir = baseDir;
    }

    public void setPrintPrevalenceAtFreq(int printPrevalenceAtFreq) {
        this.printPrevalenceAtFreq = printPrevalenceAtFreq;
    }

    public void setInfectioHistoryPrefix(String infectioHistoryPrefix) {
        this.infectioHistoryPrefix = infectioHistoryPrefix;
    }

    public void setImportFile(File importPath, String[] importSetting) {
        this.importFile = importPath;
        this.importSetting = importSetting;
    }

    public int[] getEventsPointer() {
        return pop.getEventsPointer();
    }

    public void setExportBurnInPop(int[] exportAt) {
        this.exportAt = exportAt;
    }

    public void setEventsPointer(int[] eventsPointer) {
        pop.setEventsPointer(eventsPointer);
    }

    public void addStrainIntroEnt(float[] ent) {
        strainIntroEnt.add(ent);
    }

    // <editor-fold defaultstate="collapsed" desc="InstantInfectionEvent and subclasses">    
    protected abstract class InstantInfectionEvent {

        PersonClassifier[] infClassiferArr;
        float[][] prevalArr;

        public InstantInfectionEvent(int numInf) {
            infClassiferArr = new PersonClassifier[numInf];
            prevalArr = new float[numInf][];

        }

        protected void setIntantInfectionEvent(int infId,
                PersonClassifier infClassifer, float[] preval) {
            infClassiferArr[infId] = infClassifer;
            prevalArr[infId] = preval;
        }

        public PersonClassifier getInfClassifer(int infId) {
            return infClassiferArr[infId];
        }

        public float[] getPreval(int infId) {
            return prevalArr[infId];
        }

        public boolean hasEvent(int infId) {
            return infClassiferArr[infId] != null && prevalArr[infId] != null;
        }

        public abstract void setEvent(int infId, PersonClassifier infClassifer,
                float[] preval, Object[] extra);
    }

    protected class InstantInfectionIntroEvent extends InstantInfectionEvent {

        float[][][] strainDecomposition; // new float[infId][classIndex][strain decompositon]
        int[] preExposeMax;

        public InstantInfectionIntroEvent(int numInf) {
            super(numInf);
            strainDecomposition = new float[numInf][][];
            preExposeMax = new int[numInf];
        }

        @Override
        public void setEvent(int infId, PersonClassifier infClassifer,
                float[] preval, Object[] extra) {
            setEvent(infId, infClassifer, preval, ((Number) extra[0]).intValue(),
                    (float[][]) extra[1]);
        }

        public void setEvent(int infId, PersonClassifier infClassifer,
                float[] preval, int preExposeMax, float[][] strainDecomp) {
            super.setIntantInfectionEvent(infId, infClassifer, preval);
            this.preExposeMax[infId] = preExposeMax;
            strainDecomposition[infId] = strainDecomp;
        }

        public float[][] getStrainDecomp(int infId) {
            return strainDecomposition[infId];
        }

        public int getPreExpose(int infId) {
            return preExposeMax[infId];
        }
    }

    protected class StrainIntroEvent extends InstantInfectionEvent {

        int[] strainNum;

        public StrainIntroEvent(int numInf) {
            super(numInf);
            strainNum = new int[numInf];
        }

        public void setEvent(int infId, PersonClassifier infClassifer,
                float[] preval, int strainNum) {
            super.setIntantInfectionEvent(infId, infClassifer, preval);
            this.strainNum[infId] = strainNum;
        }

        public int getStrainNum(int infId) { // Either max-preExpose or strain number
            return strainNum[infId];
        }

        @Override
        public void setEvent(int infId, PersonClassifier infClassifer, float[] preval, Object[] extra) {
            setEvent(infId, infClassifer, preval, ((Number) extra[0]).intValue());
        }
    }

    // </editor-fold>
    public SinglePopRunnable(int id, int numSnaps, int snapFreq) {
        this.id = id;
        this.numSnaps = numSnaps;
        this.snapFreq = snapFreq;
        this.snapCounts = new int[numSnaps][][];
        this.snapClassifiers = new PersonClassifier[numSnaps][];
        this.cummulativeSnap = new boolean[numSnaps][];
        this.infectionEvents = new InstantInfectionEvent[numSnaps];
        this.incidenceCounts = new int[numSnaps][];
    }

    public void setProgressSupport(PropertyChangeSupport progressSupport) {
        this.progressSupport = progressSupport;
    }

    public int getId() {
        return id;
    }

    public int getLastReportedStep() {
        return lastReportedStep;
    }

    public AbstractRegCasRelMapPopulation getPopulation() {
        return pop;
    }

    public void setPopulation(AbstractRegCasRelMapPopulation pop) {
        this.pop = pop;
    }

    public int[][][] getSnapCounts() {
        return snapCounts;
    }

    public int[][] getIncidentCounts() {
        return incidenceCounts;
    }

    public int getModelBurnIn() {
        return modelBurnIn;
    }

    public void setModelBurnIn(int modelBurnIn) {
        this.modelBurnIn = modelBurnIn;
    }

    private void modelBurnIn() {
        if (modelBurnIn > 0) {
            //long tic = System.currentTimeMillis();

            //showStrStatus("S" + getId() + ": Model burn-in for " + modelBurnIn + " steps");
            for (int t = 0; t < modelBurnIn; t++) {
                getPopulation().advanceTimeStep(1);

                //showStrStatus("S" + getId() + ": t = " + t + ", time = " + (System.currentTimeMillis() - tic) / 1000f + "s");
            }
            //showStrStatus("S" + getId() + ": Model burn-in complete.");

        }
    }

    // Export popualtion if needed
    public void exportPopAt() {
        exportPopAt(false);
    }

    public void exportPopAt(boolean forced) {

        if (baseDir != null) {

            if (forced || (exportAt != null && exportAtPt < exportAt.length)) {

                if (!forced) {
                    while (Math.abs(exportAt[exportAtPt]) < getPopulation().getGlobalTime()) {
                        exportAtPt++;
                    }
                }

                if (forced || (exportAtPt < exportAt.length
                        && Math.abs(exportAt[exportAtPt]) == getPopulation().getGlobalTime())) {
                    
                    String prefix = forced?  EXPORT_FORCED_PREFIX:"";                                                            

                    File exportPopFileDir = new File(baseDir, prefix + EXPORT_PREFIX + Integer.toString(getPopulation().getGlobalTime()));
                    exportPopFileDir.mkdirs();

                    File exportPopFileZip = new File(exportPopFileDir, SimulationInterface.POP_FILE_PREFIX + getId() + ".zip");

                    // Export pop
                    if (exportAt[exportAtPt] < 0) {

                        try {

                            File exportPopFileRaw = new File(exportPopFileDir, SimulationInterface.POP_FILE_PREFIX + getId());

                            try (ObjectOutputStream outStr = new ObjectOutputStream(new FileOutputStream(exportPopFileRaw))) {
                                getPopulation().exportPop(outStr);
                            }
                            util.FileZipper.zipFile(exportPopFileRaw, exportPopFileZip);
                            exportPopFileRaw.delete();

                        } catch (IOException ex) {
                            ex.printStackTrace(System.err);
                            showStrStatus("Error in exporting pop file " + exportPopFileZip.getAbsolutePath());

                        }
                    }

                    // Export indivdual snapshot 
                    File behaviourFile = new File(exportPopFileDir,
                            EXPORT_INDIV_PREFIX + getId() + ".csv");

                    if (getPopulation() instanceof MSMPopulation) {

                        try {

                            PrintWriter wri = new PrintWriter(new java.io.FileWriter(behaviourFile));

                            StringBuilder header = new StringBuilder("Id,Age,BehavType,# Reg,# Cas, # Cas in 6 months");

                            HashMap<Integer, int[]> currentlyVaccinated = null;

                            if (getPopulation() instanceof MSMPopulation) {
                                MSMPopulation msmPop = (MSMPopulation) getPopulation();
                                currentlyVaccinated = (HashMap<Integer, int[]>) msmPop.getValue("", MSMPopulation.MSM_SITE_CURRENTLY_VACCINATED);

                            }

                            if (currentlyVaccinated != null) {
                                header.append(',');
                                header.append("Vaccinated Until");
                            }

                            for (int p = 0; p < getPopulation().getPop().length; p++) {
                                RelationshipPerson person = (RelationshipPerson) getPopulation().getPop()[p];
                                int inReg
                                        = getPopulation().getRelMap()[MSMPopulation.MAPPING_REG].containsVertex(person.getId())
                                        ? getPopulation().getRelMap()[MSMPopulation.MAPPING_REG].degreeOf(person.getId()) : 0;

                                int inCas
                                        = getPopulation().getRelMap()[MSMPopulation.MAPPING_CAS].containsVertex(person.getId())
                                        ? getPopulation().getRelMap()[MSMPopulation.MAPPING_CAS].degreeOf(person.getId()) : 0;

                                int numCasualIn6month = 0;
                                if (person instanceof population.person.RelationshipPerson_MSM) {
                                    int[] casualRec = ((population.person.RelationshipPerson_MSM) person).getCasualRecord();
                                    for (int i = 0; i < casualRec.length; i++) {
                                        numCasualIn6month += (casualRec[i] != 0) ? 1 : 0;
                                    }
                                }

                                int[] infStat = person.getInfectionStatus();
                                int[] strainStat = person.getCurrentStrainsAtSite();

                                if (p == 0) {
                                    for (int i = 0; i < infStat.length; i++) {
                                        header.append(',');
                                        header.append("Inf Stat_" + i);
                                        header.append(',');
                                        header.append("Strain Stat_" + i);
                                    }
                                    wri.println(header.toString());
                                }

                                StringBuilder numPartnStr = new StringBuilder();
                                numPartnStr.append(person.getId());
                                numPartnStr.append(',');
                                numPartnStr.append((int) person.getAge());
                                numPartnStr.append(',');
                                numPartnStr.append(((Number) person.getParameter(person.indexToParamName(population.person.RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue());
                                numPartnStr.append(',');
                                numPartnStr.append(inReg);
                                numPartnStr.append(',');
                                numPartnStr.append(inCas);
                                numPartnStr.append(',');
                                numPartnStr.append(numCasualIn6month);
                                if (currentlyVaccinated != null) {
                                    numPartnStr.append(',');
                                    int[] res = currentlyVaccinated.get(person.getId());
                                    int vr = 0;
                                    if (res != null) {
                                        vr = res[MSMPopulation.VACC_SETTING_AGE_EXPIRY];
                                    }
                                    numPartnStr.append(vr);
                                }

                                for (int i = 0; i < infStat.length; i++) {
                                    numPartnStr.append(',');
                                    numPartnStr.append(infStat[i]);
                                    numPartnStr.append(',');
                                    numPartnStr.append(strainStat[i]);
                                }
                                wri.println(numPartnStr.toString());
                            }

                            wri.close();

                            // Zipping and remove behaviour file                            
                            File zipFile = new File(behaviourFile.getParentFile(), behaviourFile.getName() + ".zip");
                            util.FileZipper.zipFile(behaviourFile, zipFile);
                            Files.delete(behaviourFile.toPath());

                        } catch (IOException ex) {
                            ex.printStackTrace(System.err);
                        }
                    }

                    if (!forced) {
                        exportAtPt++;
                    }
                }
            }
        }
    }

    public void model_prop_initialise(int modelBurnIn, String[] model_init_val) {
        if (model_init_val != null) {
            Object[] propInitVal = new Object[model_init_val.length];
            for (int i = 0; i < model_init_val.length; i++) {
                if (model_init_val[i] != null) {
                    propInitVal[i] = util.PropValUtils.propStrToObject(model_init_val[i],
                            ((AbstractRegCasRelMapPopulation) getPopulation()).getFieldClass(i));
                }
            }
            ((AbstractRegCasRelMapPopulation) getPopulation()).loadPropertiesToPop(propInitVal);
        }
        setModelBurnIn(modelBurnIn);

        getPopulation().initialise();
        showStrStatus("Model initialised.");
    }

    public int[] getExtinctionAt() {
        return extinctionAt;
    }

    public int[][] getStrainCompositionActiveRange() {
        return strainCompositionActiveRange;
    }

    public void setSnapShotOutput(PersonClassifier[] cf, boolean[] cummulStep) {
        for (int i = 0; i < snapClassifiers.length; i++) {
            setSnapShotOutput(cf, cummulStep, i);
        }
    }

    public void setSnapShotOutput(PersonClassifier[] cf, boolean[] cummulStep, int snapShotNum) {
        snapClassifiers[snapShotNum] = cf;
        snapCounts[snapShotNum] = null;
        cummulativeSnap[snapShotNum] = null;
        if (cf != null) {
            snapCounts[snapShotNum] = new int[cf.length][];
            cummulativeSnap[snapShotNum] = new boolean[cf.length];
            System.arraycopy(cummulStep, 0, cummulativeSnap[snapShotNum], 0, cf.length);
        }

    }

    public void setInfectionIntroAt(int snapNum, int infId,
            PersonClassifier infClassifer, float[] preval, int preExposeMax,
            float[][] strainDecomposition) { // new float[classIndex][strain decompositon]
        if (infectionEvents[snapNum] == null) {
            infectionEvents[snapNum] = new InstantInfectionIntroEvent(getPopulation().getInfList().length);
        }
        infectionEvents[snapNum].setEvent(infId, infClassifer, preval,
                new Object[]{preExposeMax, strainDecomposition});
    }

    public void setInfectionIntroAt(int globalTime, int infId, float preval, int preExposeMax) {
        float[] preval_ent = infectionIntroMap.get(globalTime);
        if (preval_ent == null) {
            preval_ent = new float[infId + 1];
            infectionIntroMap.put(globalTime, preval_ent);
        } else if (preval_ent.length <= infId) {
            preval_ent = Arrays.copyOf(preval_ent, infId + 1);
            infectionIntroMap.put(globalTime, preval_ent);
        }
        preval_ent[infId] = preval;
        infectionPreExposeMap.put(infId, preExposeMax);
    }

    public void setStrainIntroAt(int snapNum, int infId,
            PersonClassifier infClassifer, float[] preval, int strainNum) {
        if (infectionEvents[snapNum] == null) {
            infectionEvents[snapNum] = new StrainIntroEvent(getPopulation().getInfList().length);
        }
        infectionEvents[snapNum].setEvent(infId, infClassifer, preval, new Object[]{new Integer(strainNum)});
    }

    private void showStrStatus(String str) {
        if (progressSupport != null) {
            progressSupport.firePropertyChange(SimulationInterface.PROGRESS_MSG, null, str);
        } else {
            System.out.println(str);
        }
    }

    private void reportStepStatus(int reportStep) {
        if (progressSupport != null) {
            progressSupport.firePropertyChange(SimulationInterface.PROGESS_CURRENT_STEP, lastReportedStep, reportStep);
            lastReportedStep = reportStep;
        }
    }

    @Override
    public void run() {
        int t;
        try {
            if (importFile != null) {
                if (importFile.isFile()) {
                    showStrStatus("S" + getId() + ": Importing population from " + importFile.getAbsolutePath());

                    File temp;

                    if (importFile.getName().endsWith(".zip")) {  // Otherwise should already be preset
                        temp = FileZipper.unzipFile(importFile, importFile.getParentFile());

                        try (java.io.ObjectInputStream objIn = new java.io.ObjectInputStream(new java.io.FileInputStream(temp))) {

                            Object[] fields = (Object[]) objIn.readObject();
                            for (int f = 0; f < fields.length; f++) {
                                if (importSetting[f] == null) { // Skip those which are overwritten
                                    getPopulation().setParameter(" ", f, fields[f]);
                                } else {
                                    Object ent = StaticMethods.propStrToObject(importSetting[f],
                                            ((AbstractRegCasRelMapPopulation) getPopulation()).getFieldClass(f));

                                    getPopulation().setParameter(" ", f, ent);
                                }
                            }

                        } catch (Exception ex) {
                            java.io.StringWriter sWri = new StringWriter();
                            java.io.PrintWriter pWri = new PrintWriter(sWri);
                            ex.printStackTrace(pWri);
                            System.err.println(sWri.toString());
                        }

                        temp.delete();
                    }

                    setModelBurnIn(0); // No burn in                       
                    getPopulation().initialiseInfection(0); // 0 = using orginal RNG

                }
            }
            // Set coexist matrix            

            AbstractInfection[] infList = getPopulation().getInfList();
            for (AbstractInfection inf : infList) {
                if (coexistMat != null) {
                    ((MultiStrainInfectionInterface) inf).setStrainCoexistMatrix(coexistMat);
                } else {
                    ((MultiStrainInfectionInterface) inf).setStrainCoexistMatrix(new float[1][1]); // Default single strain
                }
            }

            extinctionAt = new int[getPopulation().getInfList().length];
            Arrays.fill(extinctionAt, -1);

            strainCompositionActiveRange = new int[3][2]; // Start time, end time
            for (int[] ent : strainCompositionActiveRange) {
                ent[0] = Integer.MAX_VALUE;
                ent[1] = 0;
            }

            long tic = System.currentTimeMillis();
            modelBurnIn();
            showStrStatus("S" + getId() + ": Model burn-in finished after "
                    + (System.currentTimeMillis() - tic) / 1000f + "s.");

            reportStepStatus(modelBurnIn);
            t = modelBurnIn;
            showStrStatus("S" + getId() + ": Simulation in progress");

            HashMap<String, int[]> currentInfectedAt = new HashMap();
            // Str[person id, site Id], [infection status, startTime]                      

            runSim:
            for (int s = 0; s < numSnaps; s++) {
                if (infectionEvents[s] != null) {

                    for (int infId = 0; infId < getPopulation().getInfList().length; infId++) {
                        if (infectionEvents[s].hasEvent(infId)) {
                            if (infectionEvents[s] instanceof InstantInfectionIntroEvent) {
                                getPopulation().setInstantInfection(infId,
                                        infectionEvents[s].getInfClassifer(infId),
                                        infectionEvents[s].getPreval(infId),
                                        ((InstantInfectionIntroEvent) infectionEvents[s]).getPreExpose(infId),
                                        ((InstantInfectionIntroEvent) infectionEvents[s]).getStrainDecomp(infId));

                            } else if (infectionEvents[s] instanceof StrainIntroEvent) {

                                getPopulation().introduceStrains(infId,
                                        infectionEvents[s].getInfClassifer(infId),
                                        infectionEvents[s].getPreval(infId),
                                        ((StrainIntroEvent) infectionEvents[s]).getStrainNum(infId));

                            } else {
                                System.err.println(getClass().getName() + ".run: infectionEvents of class "
                                        + infectionEvents[s].getClass().getName() + "not support yet");
                            }
                        }
                    }

                }

                // Set snapshot classifiers                
                getPopulation().setSnapshotClassifier(snapClassifiers[s]);

                // Step up cummulative if needed               
                if (snapClassifiers[s] != null && snapClassifiers[s].length > 0) {
                    for (int c = 0; c < snapClassifiers[s].length; c++) {
                        if (cummulativeSnap[s][c]) {
                            snapCounts[s][c] = new int[snapClassifiers[s][c].numClass()];
                        } else {
                            snapCounts[s][c] = null;
                        }
                    }
                }

                // Snap count pre snapshot (or null if not cummulative)
                getPopulation().setSnapshotCount(snapCounts[s]);

                for (int f = 0; f < snapFreq; f++) {

                    // Set vaccine only at start and end point
                    if (vaccineSettingAll != null) {
                        for (double[] vaccSetting : vaccineSettingAll) {
                            if (getPopulation().getGlobalTime() == vaccSetting[VACCINE_START_TIME]) {
                                SiteSpecificVaccination vacc;
                                int settingStart = vaccSetting.length - (int) vaccSetting[VACCINE_SETTING_LENGTH];
                                vacc = new SiteSpecificVaccination(Arrays.copyOfRange(vaccSetting, settingStart, vaccSetting.length),
                                        new int[]{(int) vaccSetting[VACCINE_START_TIME], (int) vaccSetting[VACCINE_END_TIME]});
                                getPopulation().setParameter("Vaccine_site_specific", MSMPopulation.MSM_SITE_SPECIFIC_VACCINATION, vacc);
                            }

                        }
                    }

                    // Check for intro infection
                    if (infectionIntroMap.containsKey(getPopulation().getGlobalTime())) {
                        float[] prevalCount = infectionIntroMap.get(getPopulation().getGlobalTime());

                        // Count how many in relationship
                        int pI = 0;
                        int numInRel = 0;
                        AbstractIndividualInterface p;

                        while (pI < getPopulation().getPop().length) {
                            if (inRelationship(getPopulation().getPop()[pI])) {
                                numInRel++;
                            }
                            pI++;
                        }

                        for (int infId = 0; infId < prevalCount.length; infId++) {
                            if (prevalCount[infId] > 0) {
                                if (prevalCount[infId] < 1) {
                                    prevalCount[infId] = Math.round(prevalCount[infId] * numInRel);
                                }
                            }
                        }

                        pI = 0;

                        while (pI < getPopulation().getPop().length) {
                            p = getPopulation().getPop()[pI];

                            if (inRelationship(p)) {
                                for (int infId = 0; infId < prevalCount.length; infId++) {
                                    int preExpose = infectionPreExposeMap.get(infId) == null ? 0 : infectionPreExposeMap.get(infId).intValue();

                                    if (prevalCount[infId] > 0) {
                                        if (getPopulation().getInfList()[infId].isInfected(p)) {
                                            prevalCount[infId]--;
                                        } else {
                                            int numAvail = numInRel;
                                            if (getPopulation().getInfList()[infId].getRNG().nextInt(numAvail) < prevalCount[infId]) {
                                                getPopulation().getInfList()[infId].infecting(p);

                                                if (preExpose > 0) {
                                                    preExpose = getPopulation().getInfList()[infId].getRNG().nextInt(preExpose);
                                                }

                                                double infectAt = p.getAge() - preExpose;
                                                p.setLastInfectedAtAge(getPopulation().getInfList()[infId].getInfectionIndex(), infectAt);
                                                // Determine status immediately
                                                int stateStart = -preExpose;
                                                double cumulStageTime = p.getTimeUntilNextStage(getPopulation().getInfList()[infId].getInfectionIndex());
                                                p.setTimeUntilNextStage(getPopulation().getInfList()[infId].getInfectionIndex(), cumulStageTime + stateStart);

                                                while ((p.getTimeUntilNextStage(getPopulation().getInfList()[infId].getInfectionIndex())) < 0) {
                                                    cumulStageTime += Math.round(getPopulation().getInfList()[infId].advancesState(p));
                                                    p.setTimeUntilNextStage(getPopulation().getInfList()[infId].getInfectionIndex(), cumulStageTime + stateStart);
                                                }

                                                if (p.getInfectionStatus()[infId] != AbstractIndividualInterface.INFECT_S) {
                                                    if (getPopulation().getInfList()[infId] instanceof MultiStrainInfectionInterface
                                                            && p instanceof MultiSiteMultiStrainPersonInterface) {
                                                        ((MultiSiteMultiStrainPersonInterface) p).setCurrentStrainAtSite(infId, 1); // default - set all strain as 1

                                                    }
                                                    prevalCount[infId]--;
                                                }

                                            }
                                        }

                                    }
                                }
                                numInRel--;
                            }
                            pI++;
                        }
                    }

                    if (snapClassifiers[s] != null && f == snapFreq - 1) {
                        // Set up instant count if needed
                        for (int c = 0; c < snapClassifiers[s].length; c++) {
                            if (!cummulativeSnap[s][c]) {
                                snapCounts[s][c] = new int[snapClassifiers[s][c].numClass()];
                            }
                        }
                        getPopulation().setSnapshotCount(snapCounts[s]);
                    }

                    introStrain();
                    getPopulation().advanceTimeStep(1);
                    // Check if need to export pop
                    exportPopAt();

                    if (printPrevalenceAtFreq > 0 && (getPopulation().getGlobalTime() % printPrevalenceAtFreq == 0)) {
                        StringBuilder output = new StringBuilder();
                        output.append(this.getId());
                        output.append(',');
                        output.append(getPopulation().getGlobalTime());
                        output.append(',');
                        output.append(getPopulation().getNumInf()[0]);
                        output.append(',');
                        output.append(getPopulation().getNumInf()[1]);
                        output.append(',');
                        output.append(getPopulation().getNumInf()[2]);

                        for (int siteId = 0; siteId < getPopulation().getInfList().length; siteId++) {
                            if (getPopulation().getInfList()[siteId] instanceof MultiStrainInfectionInterface) {
                                int numStrains
                                        = ((MultiStrainInfectionInterface) (getPopulation().getInfList()[siteId])).getStrainCoexistMatrix().length;

                                int[] strainCounter = new int[numStrains];

                                for (AbstractIndividualInterface p : getPopulation().getPop()) {
                                    if (p instanceof MultiSiteMultiStrainPersonInterface) {
                                        if (p.getInfectionStatus()[siteId] != AbstractIndividualInterface.INFECT_S) {
                                            int currentS = ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite()[siteId];
                                            for (int strain = 0; strain < strainCounter.length; strain++) {
                                                if ((currentS & (1 << strain)) != 0) {
                                                    strainCounter[strain]++;
                                                }
                                            }
                                        }
                                    }
                                }

                                for (int strain = 0; strain < strainCounter.length; strain++) {
                                    output.append(',');
                                    output.append(strainCounter[strain]);
                                }

                            }

                        }

                        showStrStatus(output.toString());

                        if (printPrevalenceAtFreq == 1) {
                            System.out.println(output.toString());
                        }
                    }
                    t++;

                    if (infectioHistoryPrefix != null) {

                        int[] currentInfectedAtKey = new int[2];

                        // 20150313 - Check infection history
                        for (AbstractIndividualInterface person : getPopulation().getPop()) {

                            //HashMap<String,int[]> currentInfectedAt = new HashMap(); 
                            // String[person id, site Id], [startTime,infection status]                
                            int[] infectStat = person.getInfectionStatus();

                            currentInfectedAtKey[0] = person.getId();

                            for (int siteId = 0; siteId < infectStat.length; siteId++) {

                                currentInfectedAtKey[1] = siteId;

                                int infStat = infectStat[siteId];

                                String keyStr = Arrays.toString(currentInfectedAtKey).replaceAll("\\[", "").replaceAll("\\]", "");

                                // Has infection
                                int[] ent;
                                if (!currentInfectedAt.containsKey(keyStr)
                                        && infStat != AbstractIndividualInterface.INFECT_S
                                        && infStat != GonorrhoeaSiteInfection.STATUS_IMM) {
                                    // Not in record before                                                                
                                    ent = Arrays.copyOf(new int[]{getPopulation().getGlobalTime(), infStat}, 2);
                                    currentInfectedAt.put(keyStr, ent);

                                } else if (currentInfectedAt.containsKey(keyStr)) {

                                    if (infStat == GonorrhoeaSiteInfection.STATUS_SYM
                                            || infStat == GonorrhoeaSiteInfection.STATUS_ASY) {
                                        // Update  syptoms
                                        ent = currentInfectedAt.get(keyStr);
                                        ent[1] = infStat;

                                    } else if (infStat == AbstractIndividualInterface.INFECT_S) {
                                        // Person just recovered
                                        int[] entry = currentInfectedAt.get(keyStr);
                                        File infectionRecordFile = new File(baseDir,
                                                infectioHistoryPrefix + Integer.toString(this.getId())
                                                + ".obj");

                                        boolean fileExisted = infectionRecordFile.exists();

                                        try {

                                            java.io.FileOutputStream fileOutStream = new java.io.FileOutputStream(infectionRecordFile,
                                                    fileExisted);
                                            java.io.BufferedOutputStream bufferedOutStream = new java.io.BufferedOutputStream(fileOutStream);

                                            java.io.DataOutputStream dataOutStream = new java.io.DataOutputStream(bufferedOutStream);

                                            dataOutStream.writeInt(getPopulation().getGlobalTime());
                                            dataOutStream.writeInt(person.getId());
                                            dataOutStream.writeInt(siteId);
                                            dataOutStream.writeInt(entry[1]); // Inf Stat 
                                            dataOutStream.writeInt(entry[0]); // Start time, with duration = Global time - Start time                                         

                                            dataOutStream.close();
                                            bufferedOutStream.close();
                                            fileOutStream.close();

                                        } catch (java.io.IOException ex) {
                                            ex.printStackTrace();
                                        }

                                        currentInfectedAt.remove(currentInfectedAtKey);
                                    }
                                }

                            }

                        }
                    }

                    // Check for extinction    
                    // Only valid if there already some infection before                                                                                              
                    if (getPopulation().getNumInf() != null) {
                        boolean allExtinct = true;
                        for (int infId = 0; infId < extinctionAt.length; infId++) {
                            if (extinctionAt[infId] == -1 && getPopulation().getNumInf()[infId] == 0) {
                                extinctionAt[infId] = getPopulation().getGlobalTime();
                            } else {
                                allExtinct &= false;
                            }
                        }
                        if (allExtinct) {
                            showStrStatus("S" + getId() + ": extinction at " + t);
                            reportStepStatus(modelBurnIn + numSnaps * snapFreq);
                            // Fill the remaining step count for early extinction
                            while (s < numSnaps) {
                                if (snapClassifiers[s] != null) {
                                    for (int c = 0; c < snapClassifiers[s].length; c++) {
                                        if (!cummulativeSnap[s][c]) {
                                            snapCounts[s][c] = new int[snapClassifiers[s][c].numClass()];
                                            Arrays.fill(snapCounts[s][c], 0);
                                        } else {
                                            if (s > 1) {
                                                snapCounts[s][c] = Arrays.copyOf(snapCounts[s - 1][c], snapCounts[s - 1][c].length);
                                            } else {
                                                AbstractIndividualInterface[] allPop = getPopulation().getPop();
                                                if (snapClassifiers[s][c] != null) {
                                                    snapCounts[s][c] = new int[snapClassifiers[s][c].numClass()];
                                                    for (int p = 0; p < allPop.length; p++) {
                                                        int cType = snapClassifiers[s][c].classifyPerson(allPop[p]);
                                                        if (cType > 0 && cType < snapClassifiers[s][c].numClass()) {
                                                            snapCounts[s][c][cType]++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                s++;
                            }
                            // Force export population
                            if (exportAt != null) {
                                exportPopAt(true);
                            }
                            //break runSim;

                        }

                    }
                    // Check new strain and their elmination 

                    AbstractIndividualInterface[] allPop = getPopulation().getPop();
                    boolean[] hasStrainComposition = new boolean[strainCompositionActiveRange.length];
                    Arrays.fill(hasStrainComposition, false);

                    for (AbstractIndividualInterface p : allPop) {
                        if (p instanceof MultiSiteMultiStrainPersonInterface) {
                            int[] strainCompositionBySite = ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite();

                            for (int strainComposition : strainCompositionBySite) {
                                if (strainComposition != 0) {
                                    if (strainComposition >= strainCompositionActiveRange.length) {
                                        int org_len = strainCompositionActiveRange.length;
                                        strainCompositionActiveRange = Arrays.copyOf(strainCompositionActiveRange, strainComposition + 1);
                                        hasStrainComposition = Arrays.copyOf(hasStrainComposition, strainComposition + 1);
                                        for (int i = org_len; i < strainCompositionActiveRange.length; i++) {
                                            strainCompositionActiveRange[i] = new int[]{Integer.MAX_VALUE, 0};
                                            hasStrainComposition[i] = false;
                                        }

                                    }
                                    hasStrainComposition[strainComposition] = true;
                                    // First in pop
                                    if (strainCompositionActiveRange[strainComposition][0] > pop.getGlobalTime()) {
                                        strainCompositionActiveRange[strainComposition][0] = pop.getGlobalTime();
                                    }
                                }
                            }
                        }
                    }

                    for (int i = 0; i < hasStrainComposition.length; i++) {
                        // Extinction time
                        if (!hasStrainComposition[i]
                                && strainCompositionActiveRange[i][0] < pop.getGlobalTime()
                                && strainCompositionActiveRange[i][1] == 0) {
                            strainCompositionActiveRange[i][1] = pop.getGlobalTime();
                        }
                    }

                    // Check for surival analysis 
                    if (use_patient_zero && patient_zero != null) {

                        HashMap<Integer, int[]> patient_zero_partners_current_strain = new HashMap<>();

                        for (int r = 0; r < pop.getRelMap().length; r++) {
                            RelationshipMap relMap = pop.getRelMap()[r];
                            if (relMap.containsVertex(patient_zero.getId())) {
                                SingleRelationship[] relArr = relMap.edgesOf(patient_zero.getId()).toArray(new SingleRelationship[0]);

                                for (SingleRelationship rel1 : relArr) {
                                    int[] partners = rel1.getLinksValues();

                                    RelationshipPerson partner;

                                    if (patient_zero.getId() == partners[0]) {
                                        partner = (RelationshipPerson) pop.getPersonById(partners[1]);
                                    } else {
                                        partner = (RelationshipPerson) pop.getPersonById(partners[0]);
                                    }

                                    // Check if already in partnership record
                                    if (!patient_zero_partners_record.containsKey(partner.getId())) {
                                        int[] totalCollectionEnt = new int[LENGTH_PATIENT_ZERO_PARNTER_REC];
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_START_TIME] = pop.getGlobalTime();
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_RELATIONSHIP_LASTED_FOR] = 0;
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_TYPE] = r;

                                        if (partner instanceof population.person.RelationshipPerson_MSM) {
                                            totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_BEHAVIOR]
                                                    = ((Number) partner.getParameter(
                                                            partner.indexToParamName(population.person.RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue();
                                        }

                                        patient_zero_partners_record.put(partner.getId(), totalCollectionEnt);
                                    } else {
                                        int[] totalCollectionEnt = patient_zero_partners_record.get(partner.getId());
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_RELATIONSHIP_LASTED_FOR]++;
                                        patient_zero_partners_record.put(partner.getId(), totalCollectionEnt);
                                    }

                                    int[] mapEnt = Arrays.copyOf(partner.getCurrentStrainsAtSite(),
                                            partner.getCurrentStrainsAtSite().length);
                                    patient_zero_partners_current_strain.put(partner.getId(), mapEnt);

                                    if (hasNewStrain(partner)) {

                                        boolean newInfectFromPatientZero = false;

                                        if (patient_zero_partners_past_strain.containsKey(partner.getId())) {

                                            int[] strainStat = partner.getCurrentStrainsAtSite();
                                            int[] pastStat = patient_zero_partners_past_strain.get(partner.getId());

                                            for (int k = 0; k < strainStat.length; k++) {
                                                newInfectFromPatientZero |= (strainStat[k] & 0b10) > 0 && (pastStat[k] & 0b10) == 0;
                                            }

                                        } else {

                                            newInfectFromPatientZero = true;

                                        }

                                        if (newInfectFromPatientZero) {

                                            int[] newEnt = newEnt = new int[1 // time
                                                    + 1 // relationship type
                                                    + 2 // patient_zero_id, patient_zero_age
                                                    + patient_zero.getCurrentStrainsAtSite().length
                                                    + 2 // partner_id, partner_age
                                                    + partner.getCurrentStrainsAtSite().length];

                                            int sp = 0;
                                            newEnt[sp] = pop.getGlobalTime();
                                            sp++;
                                            newEnt[sp] = r;
                                            sp++;
                                            newEnt[sp] = patient_zero.getId();
                                            sp++;
                                            newEnt[sp] = (int) patient_zero.getAge();
                                            sp++;
                                            for (int st = 0; st < patient_zero.getCurrentStrainsAtSite().length; st++) {
                                                newEnt[sp] = patient_zero.getCurrentStrainsAtSite()[st];
                                                sp++;
                                            }
                                            newEnt[sp] = partner.getId();
                                            sp++;
                                            newEnt[sp] = (int) partner.getAge();
                                            sp++;
                                            for (int st = 0; st < partner.getCurrentStrainsAtSite().length; st++) {
                                                newEnt[sp] = partner.getCurrentStrainsAtSite()[st];
                                                sp++;
                                            }

                                            patient_zero_newStrainSpreadSummary.add(newEnt);
                                        }

                                    }

                                }

                            }
                        }

                        patient_zero_partners_past_strain = patient_zero_partners_current_strain;

                        if (!hasNewStrain(patient_zero)) {
                            patient_zero = null;
                            patient_zero_end_time = pop.getGlobalTime();

                            if (exportAt == null || exportAt.length == 0
                                    || exportAt[exportAt.length - 1] < pop.getGlobalTime()) {
                                showStrStatus("S" + getId() + ": Patient zero has no new strain. Simulation terminated.");
                                break runSim;
                            }
                        }

                    }

                }

                reportStepStatus(t);

                if (getPopulation() instanceof MSMPopulation) {
                    incidenceCounts[s] = ((MSMPopulation) getPopulation()).cumulativeIncidencesBySitesCount();
                }

                //showStrStatus("S" + getId() + ", " + getPopulation().getGlobalTime() 
                //        + ":" + Arrays.toString(getPopulation().getNumInf()));
            }

            showStrStatus("S" + getId() + ": Simulation complete."
                    + " Num infected at end = " + Arrays.toString(getPopulation().getNumInf()) + " ");
            showStrStatus("S" + getId() + ": Time required (total) " + (System.currentTimeMillis() - tic) / 1000f + "s.");

            if (use_patient_zero && patient_zero_newStrainSpreadSummary != null) {

                File strainSpreadSummaryCSV = new File(baseDir, Simulation_MSM_Population.DIR_NAMES[Simulation_MSM_Population.DIR_NEW_STRAIN_SPREAD]);

                strainSpreadSummaryCSV.mkdirs();

                strainSpreadSummaryCSV = new File(strainSpreadSummaryCSV,
                        Simulation_MSM_Population.DIR_NAMES[Simulation_MSM_Population.DIR_NEW_STRAIN_SPREAD]
                        + "_" + this.getId() + "_" + patient_zero_start_time + ".csv");

                try (PrintWriter wri = new PrintWriter(strainSpreadSummaryCSV)) {
                    wri.println("Start time," + patient_zero_start_time);
                    wri.println("End time," + patient_zero_end_time);
                    wri.println("Days with new strain," + (patient_zero_end_time - patient_zero_start_time));
                    wri.println("PATIENT_ZERO_PARNTER_ID,"
                            + "PATIENT_ZERO_PARNTER_REC_START_TIME,"
                            + "PATIENT_ZERO_PARNTER_REC_RELATIONSHIP_LASTED_FOR,"
                            + "PATIENT_ZERO_PARNTER_REC_TYPE,"
                            + "PATIENT_ZERO_PARNTER_REC_PARNTER_BEHAVIOR");

                    if (patient_zero_partners_record != null) {

                        for (Integer partner_id : patient_zero_partners_record.keySet()) {
                            int[] ent = patient_zero_partners_record.get(partner_id);
                            wri.print(partner_id);
                            for (int i = 0; i < ent.length; i++) {
                                wri.print(',');
                                wri.print(ent[i]);
                            }
                            wri.println();
                        }
                        wri.println();
                    }

                    wri.println("Time,Relationship_Type,Patient_zero_id,Patient_zero_age,Patient_zero_strain_stat,,,Target_id,Target_age,Targer_strain_stat,,,");
                    for (int[] ent : patient_zero_newStrainSpreadSummary) {
                        for (int i = 0; i < ent.length; i++) {
                            if (i != 0) {
                                wri.print(',');
                            }
                            wri.print(ent[i]);
                        }
                        wri.println();
                    }
                }
            }

            if (pop instanceof MSMPopulation) {
                MSMPopulation mPop = (MSMPopulation) pop;
                showStrStatus("S" + getId() + ": Aver rel dur = "
                        + Arrays.toString(new float[]{
                    1f * mPop.relLen[0] / mPop.relTotal[0],
                    1f * mPop.relLen[1] / mPop.relTotal[1],})
                );

                int[] casualPart_Stat = new int[10];
                int[] sum = new int[4]; // 0, 1-9, 10+, All
                String[] key = new String[]{"  0", "1-9", "10+", ""};
                for (AbstractIndividualInterface p : pop.getPop()) {
                    RelationshipPerson person = (RelationshipPerson) p;
                    int behaviour = ((Number) person.getParameter(person.indexToParamName(population.person.RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue();

                    if (behaviour != population.person.RelationshipPerson_MSM.BEHAV_REG_ONLY) {
                        int[] casualRec = ((population.person.RelationshipPerson_MSM) person).getCasualRecord();
                        int numCasual = 0;
                        for (int i = 0; i < casualRec.length; i++) {
                            numCasual += (casualRec[i] != 0) ? 1 : 0;
                        }

                        if (numCasual >= casualPart_Stat.length) {
                            casualPart_Stat = Arrays.copyOf(casualPart_Stat, numCasual + 1);
                        }

                        casualPart_Stat[numCasual]++;

                        sum[3]++;
                        if (numCasual == 0) {
                            sum[0]++;
                        } else if (numCasual >= 1 && numCasual <= 9) {
                            sum[1]++;
                        } else if (numCasual >= 10) {
                            sum[2]++;
                        }
                    }
                }

                showStrStatus("S" + getId() + ": Casual partnership distribution at the end = " + Arrays.toString(casualPart_Stat) + " Max = " + casualPart_Stat.length);
                for (int s = 0; s < sum.length - 1; s++) {
                    showStrStatus("S" + getId() + ":      " + key[s] + " = " + (100f * sum[s]) / sum[sum.length - 1] + "%");
                    //System.out.println("S" + getId() + ":      " + key[s] + " = " + (100f * sum[s]) / sum[sum.length - 1] + "%");
                }

            }

        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            try {
                java.io.File errFile = new java.io.File(baseDir, "error.log");
                java.io.FileOutputStream errStr;
                java.io.PrintWriter errWri;

                errStr = new java.io.FileOutputStream(errFile, true);
                errWri = new java.io.PrintWriter(errStr, true);
                ex.printStackTrace(errWri);
                errWri.close();
                errStr.close();
            } catch (IOException ex1) {
            }
        }

    }

    protected static boolean hasNewStrain(MultiSiteMultiStrainPersonInterface p) {
        int[] strainStat = p.getCurrentStrainsAtSite();
        boolean hasNewStrain = false;
        for (int i = 0; i < strainStat.length; i++) {
            hasNewStrain |= (strainStat[i] & 0b10) > 0;
        }
        return hasNewStrain;
    }

    public void introStrain() {
        if (strainIntroPt < strainIntroEnt.size()) {
            float[] ent = strainIntroEnt.get(strainIntroPt);

            if (ent[STRAIN_INTRO_ENTRY_INTRO_AT] == getPopulation().getGlobalTime() + 1) { // Intro in next step
                int strainId = (int) ent[STRAIN_INTRO_ENTRY_STRAIN_NUM];
                int site = (int) ent[STRAIN_INTRO_ENTRY_SITE];
                int numInfect = (int) ent[STRAIN_INTRO_ENTRY_NUM_INF];

                int behavior = -1;

                if (STRAIN_INTRO_ENTRY_BEHAVOR < ent.length) {
                    behavior = (int) ent[STRAIN_INTRO_ENTRY_BEHAVOR];
                }

                ArrayList<AbstractIndividualInterface> candidate = new ArrayList<>();

                for (AbstractIndividualInterface p : getPopulation().getPop()) {
                    boolean validCandidate = p.getInfectionStatus()[site] != AbstractIndividualInterface.INFECT_S
                            && p.getTimeUntilNextStage(site) > 1;

                    if (p instanceof population.person.RelationshipPerson_MSM) {

                        validCandidate &= (behavior < 0)
                                || behavior == ((Number) p.getParameter(
                                        ((RelationshipPerson) p).indexToParamName(
                                                population.person.RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue();
                    }

                    if (validCandidate) {
                        if (p instanceof MultiSiteMultiStrainPersonInterface) {
                            int strain = ((MultiSiteMultiStrainPersonInterface) p).getCurrentStrainsAtSite()[site];
                            if ((strain & (1 << strainId)) == 0) {
                                candidate.add(p);
                            }
                        }
                    }
                }

                int numNewInfected = Math.min(numInfect, candidate.size());

                if (numNewInfected == 0) {
                    showStrStatus("S" + getId() + " Warning: No valid candidate for imporation");
                }

                random.RandomGenerator candidSelRNG = new random.MersenneTwisterRandomGenerator(getPopulation().getInfList()[site].getRNG().nextInt());
                AbstractIndividualInterface[] canArr = candidate.toArray(new AbstractIndividualInterface[candidate.size()]);

                for (int i = 0; i < canArr.length && numNewInfected > 0; i++) {
                    if (candidSelRNG.nextInt(canArr.length - i) < numNewInfected) {

                        ((MultiSiteMultiStrainPersonInterface) canArr[i]).setCurrentStrainAtSite(site, 1 << strainId);
                        numNewInfected--;

                        showStrStatus("S" + getId() + ": switching #" + canArr[i].getId() + "'s infection at site " + site
                                + " to strain 0b" + Integer.toBinaryString(1 << strainId) + " at " + getPopulation().getGlobalTime());

                        if (use_patient_zero && patient_zero == null) {

                            // Set survival analysis                            
                            patient_zero = (RelationshipPerson) canArr[i];
                            patient_zero_start_time = pop.getGlobalTime();
                            patient_zero_partners_past_strain = new HashMap<>();
                            patient_zero_partners_record = new HashMap<>();

                            for (int r = 0; r < pop.getRelMap().length; r++) {
                                RelationshipMap relMap = pop.getRelMap()[r];
                                if (relMap.containsVertex(patient_zero.getId())) {
                                    SingleRelationship[] relArr = relMap.edgesOf(patient_zero.getId()).toArray(new SingleRelationship[0]);

                                    for (SingleRelationship rel1 : relArr) {
                                        int[] partners = rel1.getLinksValues();

                                        RelationshipPerson partner;

                                        if (patient_zero.getId() == partners[0]) {
                                            partner = (RelationshipPerson) pop.getPersonById(partners[1]);
                                        } else {
                                            partner = (RelationshipPerson) pop.getPersonById(partners[0]);
                                        }

                                        int[] mapEnt = Arrays.copyOf(partner.getCurrentStrainsAtSite(),
                                                partner.getCurrentStrainsAtSite().length);

                                        patient_zero_partners_past_strain.put(partner.getId(), mapEnt);

                                        int[] totalCollectionEnt = new int[LENGTH_PATIENT_ZERO_PARNTER_REC];
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_START_TIME] = -pop.getGlobalTime(); // - negative -> start before this
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_RELATIONSHIP_LASTED_FOR] = 0;
                                        totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_TYPE] = r;

                                        if (partner instanceof population.person.RelationshipPerson_MSM) {
                                            totalCollectionEnt[PATIENT_ZERO_PARNTER_REC_BEHAVIOR]
                                                    = ((Number) partner.getParameter(partner.indexToParamName(
                                                            population.person.RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue();
                                        }

                                        patient_zero_partners_record.put(partner.getId(), totalCollectionEnt);

                                    }

                                }
                            }

                        }
                    }
                }

                if (STRAIN_INTRO_ENTRY_FREQ < ent.length) {
                    if (ent[STRAIN_INTRO_ENTRY_FREQ] > 0) {
                        ent = Arrays.copyOf(ent, ent.length);
                        ent[STRAIN_INTRO_ENTRY_INTRO_AT] += ent[STRAIN_INTRO_ENTRY_FREQ];
                        strainIntroEnt.add(ent);
                    }
                }
                strainIntroPt++;
                introStrain(); // Recursive call if necessary
            }

        }

    }

    public boolean inRelationship(AbstractIndividualInterface p) {
        boolean inRel = false;
        RelationshipMap[] relMap = getPopulation().getRelMap();
        for (int m = 0; m < relMap.length; m++) {
            if (relMap[m].containsVertex(p.getId())) {
                inRel |= relMap[m].degreeOf(p.getId()) > 0;
            }
        }
        return inRel;
    }

    long tstamp = 0;
}
