package sim;

import infection.GonorrhoeaSiteInfection;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.HashMap;
import population.AbstractRegCasRelMapPopulation;
import population.MSMPopulation;
import person.AbstractIndividualInterface;

import util.FileZipper;
import util.PersonClassifier;
import util.StaticMethods;

/**
 *
 * @author Ben Hui
 * @version 20150407
 *
 * History:
 *
 * 20150313 - Add infection history support
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
    private final transient int[][] popSnapCounts;                    // popSnapCounts[snapNum][countType]

    private File importFile = null;
    private String[] importSetting = null;
    private File baseDir = null;

    private String infectioHistoryPrefix = null;  // If not null, then export infection history

    public void setBaseDir(File baseDir) {
        this.baseDir = baseDir;
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

    public void setEventsPointer(int[] eventsPointer) {
        pop.setEventsPointer(eventsPointer);
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
        this.popSnapCounts = new int[numSnaps][];
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

    public int[][] getPopSnapCounts() {
        return popSnapCounts;
    }

    public int getModelBurnIn() {
        return modelBurnIn;
    }

    public void setModelBurnIn(int modelBurnIn) {
        this.modelBurnIn = modelBurnIn;
    }

    private void modelBurnIn() {
        if (modelBurnIn > 0) {
            showStrStatus("S" + getId() + ": Model burn-in for " + modelBurnIn + " steps");
            for (int t = 0; t < modelBurnIn; t++) {
                getPopulation().advanceTimeStep(1);
            }
            showStrStatus("S" + getId() + ": Model burn-in complete.");
        }
    }

    public void model_prop_initialise(int modelBurnIn, String[] model_init_val) {
        if (model_init_val != null) {
            Object[] propInitVal = new Object[model_init_val.length];
            for (int i = 0; i < model_init_val.length; i++) {
                if (model_init_val[i] != null) {
                    propInitVal[i] = StaticMethods.propStrToObject(model_init_val[i],
                            ((AbstractRegCasRelMapPopulation) getPopulation()).getFieldClass(i));
                }
            }
            ((AbstractRegCasRelMapPopulation) getPopulation()).loadPropertiesToPop(propInitVal);
        }
        setModelBurnIn(modelBurnIn);
        getPopulation().initialise();
    }

    public int[] getExtinctionAt() {
        return extinctionAt;
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

            extinctionAt = new int[getPopulation().getInfList().length];
            Arrays.fill(extinctionAt, -1);

            modelBurnIn();
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
                    if (snapClassifiers[s] != null && f == snapFreq - 1) {
                        // Set up instant count if needed
                        for (int c = 0; c < snapClassifiers[s].length; c++) {
                            if (!cummulativeSnap[s][c]) {
                                snapCounts[s][c] = new int[snapClassifiers[s][c].numClass()];
                            }
                        }
                        getPopulation().setSnapshotCount(snapCounts[s]);
                    }

                    getPopulation().advanceTimeStep(1);
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
                            break runSim;
                        }
                    }

                }

                reportStepStatus(t);
                popSnapCounts[s] = getPopulation().generatedPopSnapCount();
            }

            showStrStatus("S" + getId() + ": Simulation complete."
                    + " Num infected at end = " + Arrays.toString(getPopulation().getNumInf()));

            if (pop instanceof MSMPopulation) {
                MSMPopulation mPop = (MSMPopulation) pop;
                showStrStatus("S" + getId() + ": Aver rel dur = "
                        + Arrays.toString(new float[]{
                            1f * mPop.relLen[0] / mPop.relTotal[0],
                            1f * mPop.relLen[1] / mPop.relTotal[1],})
                );
            }

        } catch (Exception ex) {
            ex.printStackTrace(System.err);
            try {
                java.io.File errFile = new java.io.File("error.log");
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

    long tstamp = 0;
}
