/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package opt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import transform.ParameterConstraintTransform;
import transform.ParameterConstraintTransformSineCurve;

/**
 *
 * @author Bhui
 */
public abstract class Abstract_Optimisation_MSM {

    public static final String FILENAME_PARAM_CONSTRIANTS = "ParamConstriants.csv";
    public static final String FILENAME_OPT_RESULTS_CSV = "ParamOpt.csv";
    public static final String FILENAME_OPT_RESULTS_OBJ = "ParamOpt.obj";
    public static final String FILENAME_P0 = "Pre_P0.csv";
    public File baseDir = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Optimsation");
    public double[] NUM_INF_TARGET = new double[]{20, 830, 860};
    public double[] FITTING_WEIGHT = new double[]{1, 1, 1};
    public double[] NUM_INF_LB = null;
    public double[] NUM_INF_UB = null;
    public double[] OPT_SPECIFIC_PARAM = null; // Currently jsut pool size

    public int NUM_THREAD = Runtime.getRuntime().availableProcessors();
    public int NUM_SIM = 1;
    public long BASESEED = 6109418859537162492l;

    public Abstract_Optimisation_MSM(String[] arg) {
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
                int offset = NUM_INF_TARGET.length;

                if (val.length > offset) {
                    for (int i = 0; i < FITTING_WEIGHT.length; i++) {
                        FITTING_WEIGHT[i] = Double.parseDouble(val[offset + i]);
                    }
                }

                offset += FITTING_WEIGHT.length;
                if (val.length > offset) {
                    NUM_INF_LB = Arrays.copyOf(NUM_INF_TARGET, NUM_INF_TARGET.length);
                    for (int i = 0; i < NUM_INF_LB.length; i++) {
                        NUM_INF_LB[i] = Double.parseDouble(val[offset + i]);
                    }
                }

                offset += NUM_INF_LB.length;
                if (val.length > offset) {
                    NUM_INF_UB = Arrays.copyOf(NUM_INF_TARGET, NUM_INF_TARGET.length);
                    for (int i = 0; i < NUM_INF_LB.length; i++) {
                        NUM_INF_UB[i] = Double.parseDouble(val[offset + i]);
                    }
                }
                
                offset += NUM_INF_UB.length;
                if (val.length > offset) {
                    OPT_SPECIFIC_PARAM = new double[val.length - offset];
                    for (int i = 0; i < OPT_SPECIFIC_PARAM.length; i++) {
                        OPT_SPECIFIC_PARAM[i] = Double.parseDouble(val[offset + i]);
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

    public abstract void runOptimisation() throws IOException, ClassNotFoundException;

    public File getBaseDir() {
        return baseDir;
    }

    public double[] getNUM_INF_TARGET() {
        return NUM_INF_TARGET;
    }

    public double[] getFITTING_WEIGHT() {
        return FITTING_WEIGHT;
    }

    public double[] getNUM_INF_LB() {
        return NUM_INF_LB;
    }

    public double[] getNUM_INF_UB() {
        return NUM_INF_UB;
    }

    protected ParameterConstraintTransform[] initaliseContriants() throws IOException, NumberFormatException {
        ParameterConstraintTransform[] constraints;
        //<editor-fold defaultstate="collapsed" desc="Intialise constraints">
        File costrainFile = new File(baseDir, FILENAME_PARAM_CONSTRIANTS);
        try (final BufferedReader constraintReader = new BufferedReader(new FileReader(costrainFile))) {
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
                constraints[lnNum] = new ParameterConstraintTransformSineCurve(new double[]{Double.parseDouble(ent[0]), Double.parseDouble(ent[1])});
                /*
                constraints[lnNum] = new transform.ParameterConstraintTransformLinear(new double[]{
                Double.parseDouble(ent[0]), Double.parseDouble(ent[1])});
                 */
                lnNum++;
            }
        }
        //</editor-fold>
        return constraints;
    }

}
