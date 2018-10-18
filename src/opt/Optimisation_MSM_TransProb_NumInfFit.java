package opt;

import java.beans.PropertyChangeSupport;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import optimisation.AbstractParameterOptimiser;
import optimisation.AbstractResidualFunc;
import optimisation.NelderMeadOptimiser;
import transform.ParameterConstraintTransform;
import population.MSMPopulation;
import population.person.RelationshipPerson_MSM;
import static sim.SimulationInterface.PROP_BURNIN;
import static sim.SimulationInterface.PROP_INFECTION_INTRO;
import static sim.SimulationInterface.PROP_NUM_SNAP;
import static sim.SimulationInterface.PROP_SNAP_FREQ;
import sim.Simulation_MSM_Population;
import sim.SinglePopRunnable;

/**
 * Perform parameter optimisation of MSM model base on number of infection
 *
 * @author Ben Hui
 * @version 20181018
 *
 * <pre>
 * History:
 * 20181011 - Renaming of class and add method for 7 tranmission parameter
 * 20181012 - Add dummy PropertyChangeSupport to reduce textual output
 * 20181018 - Add target weight option 
 * </pre>
 *
 */
public class Optimisation_MSM_TransProb_NumInfFit {

    public static final String FILENAME_PARAM_CONSTRIANTS = "ParamConstriants.csv";
    public static final String FILENAME_P0 = "Pre_P0.csv";
    public static final String FILENAME_OPT_RESULTS_CSV = "ParamOpt.csv";
    public static final String FILENAME_OPT_RESULTS_OBJ = "ParamOpt.obj";
    public static final String FILENAME_OPT_SIMPLEX = "ParamSimplex.obj";

    public File baseDir = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Optimsation");
    public double[] NUM_INF_TARGET = new double[]{100, 460, 400};
    public double[] FITTING_WEIGHT = new double[]{1,1,1};
    public int NUM_SIM = 1;
    public long BASESEED = 6109418859537162492l;
    public int NUM_THREAD = Runtime.getRuntime().availableProcessors();

    public Optimisation_MSM_TransProb_NumInfFit(String[] arg) {
        if (arg.length > 0) {
            // ImportPath
            if (!arg[0].isEmpty()) {
                baseDir = new File(arg[0]);
            }
            // Target
            if (!arg[1].isEmpty()) {
                String[] val = arg[1].split(",");
                for (int i = 0; i < NUM_INF_TARGET.length; i++) {
                    NUM_INF_TARGET[i] = Double.parseDouble(val[i]);
                }
                if(val.length > NUM_INF_TARGET.length){
                    for (int i = 0; i < FITTING_WEIGHT.length; i++) {
                        FITTING_WEIGHT[i] = Double.parseDouble(val[NUM_INF_TARGET.length + i]);
                    }
                    
                }
                
            }
            // Number of sim 
            if (!arg[2].isEmpty()) {
                NUM_SIM = Integer.parseInt(arg[2]);
            }
            // Base seed
            if (!arg[3].isEmpty()) {
                BASESEED = Long.parseLong(arg[3]);
            }
            // Number of tread to use 
            if (!arg[4].isEmpty()) {
                NUM_THREAD = Integer.parseInt(arg[4]);
            }
        }
    }

    public void runOptimisation() throws IOException, ClassNotFoundException {

        ParameterConstraintTransform[] constraints;
        AbstractResidualFunc optimisationFunc;

        //<editor-fold defaultstate="collapsed" desc="Intialise constraints">   
        File costrainFile = new File(baseDir, FILENAME_PARAM_CONSTRIANTS);
        try (BufferedReader constraintReader = new BufferedReader(new FileReader(costrainFile))) {
            int lnNum = 0;
            String line;
            while (constraintReader.readLine() != null) {
                lnNum++;
            }
            constraints = new ParameterConstraintTransform[lnNum];
            
            lnNum = 0;
            BufferedReader constraintReader2 = new BufferedReader(new FileReader(costrainFile));

            while ((line = constraintReader2.readLine()) != null) {
                String[] ent = line.split(",");
                
                constraints[lnNum] = new transform.ParameterConstraintTransformSineCurve(new double[]{
                    Double.parseDouble(ent[0]), Double.parseDouble(ent[1])});
                
                
                /*
                constraints[lnNum] = new transform.ParameterConstraintTransformLinear(new double[]{
                    Double.parseDouble(ent[0]), Double.parseDouble(ent[1])});
                */
                
                lnNum++;
            }
        }
        //</editor-fold>
  

        optimisationFunc = new Opt_ResidualFunc(baseDir);

        AbstractParameterOptimiser opt = new NelderMeadOptimiser(optimisationFunc);

        // Initial value              
        double[] p0 = null;

        // Default p0, to be replace by imported if necessary
        File preP0 = new File(baseDir, FILENAME_P0);

        if (preP0.exists()) {
            ArrayList<String> p0_Arr = new ArrayList<>();
            try (BufferedReader p0Reader = new BufferedReader(new FileReader(preP0))) {
                String line;
                while ((line = p0Reader.readLine()) != null) {
                    p0_Arr.add(line);
                }
            }
            p0 = new double[p0_Arr.size()];
            int index = 0;
            for (String ent : p0_Arr) {
                p0[index] = Double.parseDouble(ent);
                index++;
            }

            System.out.println("P0 from " + preP0.getAbsolutePath() + " imported");
        }
        
        System.out.println("Fitting target = " + Arrays.toString(NUM_INF_TARGET));
        System.out.println("Weighting = " + Arrays.toString(FITTING_WEIGHT));

        //<editor-fold defaultstate="collapsed" desc="Optimisation process">    
        double[] r0 = null;
        opt.setResOptions(false, AbstractParameterOptimiser.RES_OPTIONS_PRINT);
        opt.setFilename(baseDir.getAbsolutePath() + File.separator + FILENAME_OPT_RESULTS_CSV, AbstractParameterOptimiser.FILEPATH_CSV);
        opt.setFilename(baseDir.getAbsolutePath() + File.separator + FILENAME_OPT_RESULTS_OBJ, AbstractParameterOptimiser.FILEPATH_OBJ);
        opt.setFilename(baseDir.getAbsolutePath() + File.separator + FILENAME_OPT_SIMPLEX, NelderMeadOptimiser.FILEPATH_SIMPLEX);

        File preSimplexFile = new File(baseDir, FILENAME_OPT_SIMPLEX);

        double[][] sX = null;
        double[][] sR = null;
        if (preSimplexFile.exists()) {
            System.out.print("Reading previous simplex....");
            try (ObjectInputStream objStr = new ObjectInputStream(new FileInputStream(preSimplexFile))) {
                sX = (double[][]) objStr.readObject();
                sR = (double[][]) objStr.readObject();
            }
            preSimplexFile.renameTo(new File(baseDir, preSimplexFile.getName() + "_" + System.currentTimeMillis()));
            System.out.println(" done");
        }

        if (sX != null) {
            r0 = sR[0];
            if (p0 == null) {
                p0 = new double[sX[0].length];
            }
            for (int i = 0; i < p0.length; i++) {
                p0[i] = sX[0][i];
                if (constraints[i] != null) {
                    p0[i] = constraints[i].toContrainted(p0[i]);
                }
            }
        }

        System.out.println(
                "P0 = " + Arrays.toString(p0));
        if (r0 != null) {
            System.out.println("R0 = " + Arrays.toString(r0));
        }

        opt.setP0(p0, constraints);

        opt.setR0(r0);

        if (sX != null) {
            for (int i = 0; i < sX.length; i++) {
                if (sX[i] != null) {
                    System.out.println("Loading simplex vertex #" + i);
                    System.out.println(i + ": X = " + Arrays.toString(sX[i]));
                    if (sR[i] == null) {
                        System.out.println(i + ": R to be generated");
                    }
                    sR[i] = ((NelderMeadOptimiser) opt).setPreDefineSimplex(sX[i], sR[i], i);
                    System.out.println(i + ": R = " + Arrays.toString(sR[i]));
                }
            }
        }

        //</editor-fold>
        opt.initialise();
        opt.optimise();

    }

    private class Opt_ResidualFunc extends AbstractResidualFunc {

        File baseDir;
        Object[] defaultProbVal;
        String[] defaultModelInit;

        private Opt_ResidualFunc(File baseDir) throws IOException {
            this.baseDir = baseDir;
            Simulation_MSM_Population sim = new Simulation_MSM_Population();

            Properties prop;
            prop = new Properties();

            Path propFile = new File(baseDir, Simulation_MSM_Population.FILENAME_PROP).toPath();
            try (InputStream inStr = java.nio.file.Files.newInputStream(propFile)) {
                prop.loadFromXML(inStr);
            }
            sim.loadProperties(prop);

            defaultProbVal = sim.getPropVal();
            defaultModelInit = sim.getPropModelInit();

        }

        @Override
        public double[] generateResidual(double[] param) {

            Object[] propVal = Arrays.copyOf(defaultProbVal, defaultProbVal.length);
            long[] seed;

            if (NUM_SIM <= 1) {
                seed = new long[]{BASESEED};
            } else {
                seed = new long[NUM_SIM];
                random.MersenneTwisterRandomGenerator RNG = new random.MersenneTwisterRandomGenerator(BASESEED);
                for (int i = 0; i < seed.length; i++) {
                    seed[i] = RNG.nextLong();
                }
            }

            SinglePopRunnable[] runnable = new SinglePopRunnable[seed.length];
            int[][] runInfected = new int[runnable.length][];
            int[] init_expose = new int[]{3, 12 * 30, 12 * 7};

            ExecutorService exec = null;

            int numInExe = 0;

            for (int r = 0; r < runnable.length; r++) {
                if (NUM_THREAD > 1 && exec == null) {
                    exec = Executors.newFixedThreadPool(NUM_THREAD);
                }

                runnable[r] = new SinglePopRunnable(r, ((Number) propVal[PROP_NUM_SNAP]).intValue(),
                        ((Number) propVal[PROP_SNAP_FREQ]).intValue());
                runnable[r].setPopulation(new MSMPopulation(seed[r]));

                String[] model_init_val = Arrays.copyOf(defaultModelInit, defaultModelInit.length);

                double[][] field_transmit = (double[][]) util.PropValUtils.propStrToObject(model_init_val[MSMPopulation.FIELDS_TRANSMIT], double[][].class);
                double[][] field_suscept = (double[][]) util.PropValUtils.propStrToObject(model_init_val[MSMPopulation.FIELDS_SUSCEPT], double[][].class);

                // Modified tranmission and susceptibilty behavoir
                // 0: G to A
                // 1: A to G
                // 2: G to R
                // 3: R to G
                // 4: A to R
                // 5: R to A
                // 6: R to R
                // Full 7                
                field_transmit[RelationshipPerson_MSM.SITE_G][0] = 1;  // Baseline
                field_suscept[RelationshipPerson_MSM.SITE_G][0] = 1;

                // 0: G to A
                field_suscept[RelationshipPerson_MSM.SITE_A][0] = param[0];
                // 1: A to G
                field_transmit[RelationshipPerson_MSM.SITE_A][0] = param[1];
                // 2: G to R
                field_suscept[RelationshipPerson_MSM.SITE_R][0] = param[2];
                // 3: R to G
                field_transmit[RelationshipPerson_MSM.SITE_R][0] = param[3];
                // 4: A to R
                field_transmit[MSMPopulation.TRAN_SUSC_INDEX_RIMMING_ANAL][0] = param[4];
                field_suscept[MSMPopulation.TRAN_SUSC_INDEX_RIMMING_ORAL][0] = 1;
                // 5: R to A
                field_transmit[MSMPopulation.TRAN_SUSC_INDEX_RIMMING_ORAL][0] = param[5];
                field_suscept[MSMPopulation.TRAN_SUSC_INDEX_RIMMING_ANAL][0] = 1;
                // 6: R to R
                field_transmit[MSMPopulation.TRAN_SUSC_INDEX_KISSING][0] = param[6];
                field_suscept[MSMPopulation.TRAN_SUSC_INDEX_KISSING][0] = 1;

                model_init_val[MSMPopulation.FIELDS_TRANSMIT] = util.PropValUtils.objectToPropStr(field_transmit, field_transmit.getClass());
                model_init_val[MSMPopulation.FIELDS_SUSCEPT] = util.PropValUtils.objectToPropStr(field_suscept, field_suscept.getClass());

                runnable[r].model_prop_initialise(((Number) propVal[PROP_BURNIN]).intValue(), model_init_val);

                if (propVal[PROP_INFECTION_INTRO] != null) {

                    float[][] preval_intro = (float[][]) propVal[PROP_INFECTION_INTRO];

                    for (int infId = 0; infId < preval_intro.length; infId++) {
                        for (int t = 0; t < preval_intro[infId].length; t += 2) {
                            runnable[r].setInfectionIntroAt((int) preval_intro[infId][t],
                                    infId, preval_intro[infId][t + 1], init_expose[infId]);
                        }
                    }
                }

                if (NUM_THREAD <= 1 || exec == null) {
                    runnable[r].run();

                } else {
                    runnable[r].setProgressSupport(
                            new PropertyChangeSupport(runnable[r]) {
                        @Override
                        public void firePropertyChange(String string, Object o, Object o2) {
                            // Dummy - do nothing
                        }

                    }
                    );

                    exec.submit(runnable[r]);
                    numInExe++;
                }

                if (numInExe == NUM_THREAD) {
                    try {
                        exec.shutdown();
                        if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
                            System.err.println("Inf Thread time-out!");
                        }
                    } catch (InterruptedException ex) {
                        ex.printStackTrace(System.err);
                    }
                    exec = null;
                }
            }

            // Final run (if needed)
            if (exec != null) {
                try {
                    exec.shutdown();
                    if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
                        System.err.println("Inf Thread time-out!");
                    }
                } catch (InterruptedException ex) {
                    ex.printStackTrace(System.err);
                }
                exec = null;
            }

            for (int r = 0; r < runnable.length; r++) {
                runInfected[r] = runnable[r].getPopulation().getNumInf();
            }

            double[] res = new double[3];

            for (int r = 0; r < runnable.length; r++) {
                for (int i = 0; i < res.length; i++) {
                    res[i] += ((runInfected[r][i] - NUM_INF_TARGET[i]) * FITTING_WEIGHT[i]) / runnable.length;
                }
            }

            return res;

        }

    }

}
