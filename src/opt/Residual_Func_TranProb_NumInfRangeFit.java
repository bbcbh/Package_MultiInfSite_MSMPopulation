package opt;

import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import optimisation.AbstractResidualFunc;
import population.MSMPopulation;
import population.person.RelationshipPerson_MSM;
import static sim.SimulationInterface.PROP_BURNIN;
import static sim.SimulationInterface.PROP_INFECTION_INTRO;
import static sim.SimulationInterface.PROP_NUM_SNAP;
import static sim.SimulationInterface.PROP_SNAP_FREQ;
import sim.Simulation_MSM_Population;
import static sim.Simulation_MSM_Population.PROP_MSM_INIT_POP_SIZE;
import sim.SinglePopRunnable;

/**
 *
 * @author Ben Hui
 */
public class Residual_Func_TranProb_NumInfRangeFit extends AbstractResidualFunc {

    File baseDir;
    Object[] defaultProbVal;
    String[] defaultModelInit;
    public int NUM_SIM_PER_RESIDUAL = 1;
    public long BASESEED = 6109418859537162492l;
    public int THREAD_POOL_SIZE_PER_RESIDUAL = Runtime.getRuntime().availableProcessors();
    Abstract_Optimisation_MSM src;

    public Residual_Func_TranProb_NumInfRangeFit(Abstract_Optimisation_MSM src) throws IOException {

        this.src = src;
        this.baseDir = src.getBaseDir();
        Simulation_MSM_Population sim = new Simulation_MSM_Population();

        Properties prop;
        prop = new Properties();

        Path propFile = new File(baseDir, Simulation_MSM_Population.FILENAME_PROP).toPath();
        try (InputStream inStr = java.nio.file.Files.newInputStream(propFile)) {
            prop.loadFromXML(inStr);
        }
        sim.loadProperties(prop);

        NUM_SIM_PER_RESIDUAL = sim.getPropVal()[Simulation_MSM_Population.PROP_NUM_SIM_PER_SET]
                != null ? (Integer) (sim.getPropVal()[Simulation_MSM_Population.PROP_NUM_SIM_PER_SET]) : NUM_SIM_PER_RESIDUAL;
        BASESEED = sim.getPropVal()[Simulation_MSM_Population.PROP_BASESEED]
                != null ? ((Long) sim.getPropVal()[Simulation_MSM_Population.PROP_BASESEED]) : BASESEED;

        defaultProbVal = sim.getPropVal();
        defaultModelInit = sim.getPropModelInit();

    }

    public Residual_Func_TranProb_NumInfRangeFit(Abstract_Optimisation_MSM src,
            int numSim_forced, int num_sim_per_residual) throws IOException {
        this(src);
        this.NUM_SIM_PER_RESIDUAL = numSim_forced;
        this.THREAD_POOL_SIZE_PER_RESIDUAL = num_sim_per_residual;
    }

    @Override
    public double[] generateResidual(double[] param) {
        Object[] propVal = Arrays.copyOf(defaultProbVal, defaultProbVal.length);

        long[] seed;

        if (preset_Seed == null) {

            if (NUM_SIM_PER_RESIDUAL <= 1) {
                seed = new long[]{BASESEED};
            } else {
                seed = new long[NUM_SIM_PER_RESIDUAL];
                random.MersenneTwisterRandomGenerator RNG = new random.MersenneTwisterRandomGenerator(BASESEED);
                for (int i = 0; i < seed.length; i++) {
                    seed[i] = RNG.nextLong();
                }
            }
        } else {
            seed = Arrays.copyOf(preset_Seed, NUM_SIM_PER_RESIDUAL);
        }

        SinglePopRunnable[] runnable = new SinglePopRunnable[seed.length];
        int[][] runInfected = new int[runnable.length][];
        int[] init_expose = new int[]{3, 12 * 30, 12 * 7};

        ExecutorService exec = null;

        int numInExe = 0;

        for (int r = 0; r < runnable.length; r++) {
            if (THREAD_POOL_SIZE_PER_RESIDUAL > 1 && exec == null) {
                exec = Executors.newFixedThreadPool(THREAD_POOL_SIZE_PER_RESIDUAL);
            }

            runnable[r] = new SinglePopRunnable(r, ((Number) propVal[PROP_NUM_SNAP]).intValue(),
                    ((Number) propVal[PROP_SNAP_FREQ]).intValue());
            runnable[r].setPopulation(new MSMPopulation(seed[r]));

            if (propVal[PROP_MSM_INIT_POP_SIZE] != null) {
                ((MSMPopulation) runnable[r].getPopulation()).setInitNumInPop((Integer) propVal[PROP_MSM_INIT_POP_SIZE]);
            }

            runnable[r].setPrintPrevalenceAtFreq(-1);
            runnable[r].setProgressSupport(new PropertyChangeSupport(runnable[r])); // Null response

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

            if (THREAD_POOL_SIZE_PER_RESIDUAL <= 1 || exec == null) {
                runnable[r].run();

            } else {
                runnable[r].setProgressSupport(
                        new PropertyChangeSupport(runnable[r]) {
                    /**
							 * 
							 */
							private static final long serialVersionUID = 1316198029768833713L;

					@Override
                    public void firePropertyChange(String string, Object o, Object o2) {
                        // Dummy - do nothing
                    }

                }
                );

                exec.submit(runnable[r]);
                numInExe++;
            }

            if (numInExe == THREAD_POOL_SIZE_PER_RESIDUAL) {
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

        double[] res = new double[3 * runnable.length];

        int p = 0;
        for (int r = 0; r < runnable.length; r++) {
            for (int i = 0; i < runInfected[r].length; i++) {
                res[p] = (runInfected[r][i] - src.getNUM_INF_TARGET()[i]) * src.getFITTING_WEIGHT()[i];
                p++;
            }
        }

        // System.out.println("P = " + Arrays.toString(param) + " R = " + Arrays.toString(res));
        return res;

    }

}
