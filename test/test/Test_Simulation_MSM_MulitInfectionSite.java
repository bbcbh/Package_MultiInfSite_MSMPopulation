package test;

import java.io.IOException;
import sim.Simulation_MSM_Population_BatchRun;

public class Test_Simulation_MSM_MulitInfectionSite {

    public static void main(String[] arg) throws IOException, InterruptedException, ClassNotFoundException {

        String[] grpDirNames = new String[]{
            "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test",
            //"Baseline",
            //"Vacc_Ideal","Vacc_50TR","Vacc_50SU",            
            //"Vacc_Dur_1800", "Vacc_Dur_720",
            //"Vacc_50R","Vacc_00R","Vacc_00R00A",            
            //"Import_1_A_B",
            //"Import_1_A_R_10000",
            //"Import_1_A_R50", 
            //"Import_1_A_R_B2",
            //"Import_1_A_R_B1",
            //"Import_1_A_R_B0", 
            //"Optim_CLG",                                   
            //"Optim_Preval_Range_Fit",       
            "Vacc_test",
            //"Vacc_C30_Ideal_D720", "Vacc_C30_R100_D720"
            
        };

        Simulation_MSM_Population_BatchRun.main(grpDirNames);

    }

}
