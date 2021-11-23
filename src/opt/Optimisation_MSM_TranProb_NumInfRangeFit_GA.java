package opt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import optimisation.AbstractParameterOptimiser;
import optimisation.AbstractResidualFunc;
import optimisation.GeneticAlgorithmOptimiser;
import transform.ParameterConstraintTransform;

/**
 * Perform parameter optimisation using GA Optimisiser
 *
 * @author Ben Hui
 * @version 20181026
 *
 */
public class Optimisation_MSM_TranProb_NumInfRangeFit_GA extends Abstract_Optimisation_MSM {

    private int POOL_SIZE = 1000;
    public static final String FILENAME_OPT_GA_STORE = "GA_POP.obj";

    public Optimisation_MSM_TranProb_NumInfRangeFit_GA(String[] arg) {
        super(arg);
        if (OPT_SPECIFIC_PARAM != null) {
            POOL_SIZE = (int) OPT_SPECIFIC_PARAM[0];
        }

    }

    @Override
    public void runOptimisation() throws IOException, ClassNotFoundException {
        ParameterConstraintTransform[] constraints;
        AbstractResidualFunc optimisationFunc;

        constraints = initaliseContriants();
        optimisationFunc = new Residual_Func_TranProb_NumInfRangeFit(this, NUM_SIM, NUM_THREAD); // Using 1 threads per parameter sample                

        AbstractParameterOptimiser opt = new GeneticAlgorithmOptimiser(optimisationFunc);

        File parentDir = this.getBaseDir();
        parentDir.mkdirs();

        opt.setFilename(parentDir.getAbsolutePath() + File.separator + FILENAME_OPT_RESULTS_CSV, AbstractParameterOptimiser.FILEPATH_CSV);
        opt.setFilename(parentDir.getAbsolutePath() + File.separator + FILENAME_OPT_RESULTS_OBJ, AbstractParameterOptimiser.FILEPATH_OBJ);
        opt.setParameter(GeneticAlgorithmOptimiser.PARAM_GA_OPT_POP_FILE, new File(parentDir, FILENAME_OPT_GA_STORE));
        opt.setParameter(GeneticAlgorithmOptimiser.PARAM_GA_OPT_USE_PARALLEL, NUM_THREAD);
        opt.setParameter(GeneticAlgorithmOptimiser.PARAM_GA_OPT_POP_SIZE, POOL_SIZE);
        opt.setParameter(GeneticAlgorithmOptimiser.PARAM_GA_NUM_SEED_PER_GA_POP_ENT, NUM_SIM);

        opt.setResOptions(false, AbstractParameterOptimiser.RES_OPTIONS_PRINT);

        System.out.println("# Parameters = " + constraints.length);

        // Initial value             
        double[] p0 = null, r0 = null;

        // Default p0, to be replace by imported if necessary
        p0 = new double[constraints.length];
        r0 = new double[this.getNUM_INF_TARGET().length * NUM_SIM];

        Arrays.fill(r0, Double.NaN);

        // Set R0 tol - if all fit within tol then optimisation stop
        double[][] r0_bounds = new double[2][r0.length];

        for (int i = 0; i < r0.length; i += getNUM_INF_TARGET().length) {
            for (int t = 0; t < getNUM_INF_TARGET().length; t++) {
                r0_bounds[0][i + t] = (getNUM_INF_LB()[t] - getNUM_INF_TARGET()[t]) * getFITTING_WEIGHT()[t];
                r0_bounds[1][i + t] = (getNUM_INF_UB()[t] - getNUM_INF_TARGET()[t]) * getFITTING_WEIGHT()[t];
            }
        }
        opt.setParameter(GeneticAlgorithmOptimiser.PARAM_GA_RO_TOL, r0_bounds);

        // Default p0, to be replace by imported if necessary
        File preGAPop = new File(parentDir, FILENAME_OPT_GA_STORE);
        if (!preGAPop.exists()) {
            File preP0 = new File(parentDir, FILENAME_P0);

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
                //System.out.println("R0 to be generated within GA.");

            }

        }

        opt.setP0(p0, constraints);
        opt.setR0(r0);

        //< /editor-fold>
        opt.initialise();

        opt.optimise();

    }

}
