package opt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import optimisation.AbstractParameterOptimiser;
import optimisation.AbstractResidualFunc;
import optimisation.NelderMeadOptimiser;
import transform.ParameterConstraintTransform;

/**
 * Perform parameter optimisation of MSM model base on number of infection
 *
 * @author Ben Hui
 * @version 20181026
 *
 * <pre>
 * History:
 * 20181011 - Renaming of class and add method for 7 tranmission parameter
 * 20181012 - Add dummy PropertyChangeSupport to reduce textual output
 * 20181018 - Add target weight option
 * 20181026 - Renaming and add abstract class for optimisation
 * </pre>
 *
 */
public class Optimisation_MSM_TransProb_NumInfFit_NM extends Abstract_Optimisation_MSM {

    public static final String FILENAME_OPT_SIMPLEX = "ParamSimplex.obj";

    public Optimisation_MSM_TransProb_NumInfFit_NM(String[] arg) {
        super(arg);
    }

    @Override
    public void runOptimisation() throws IOException, ClassNotFoundException {

        ParameterConstraintTransform[] constraints;
        AbstractResidualFunc optimisationFunc;
        
        constraints = initaliseContriants();
        optimisationFunc = new Residual_Func_TranProb_NumInfFit(this);

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


}
