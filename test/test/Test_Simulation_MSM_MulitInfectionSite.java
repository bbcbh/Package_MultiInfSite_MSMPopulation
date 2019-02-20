package test;

import java.io.IOException;
import sim.Simulation_MSM_Population_BatchRun;

public class Test_Simulation_MSM_MulitInfectionSite {

    public static void main(String[] arg) throws IOException, InterruptedException, ClassNotFoundException {

        String[] grpDirNames = new String[]{
            "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test",
            //"Import_1_A_B",
            //"Import_1_A_R",
            //"Import_1_A_R50", 
            "Import_1_A_R_B2",
            "Import_1_A_R_B1",
            "Import_1_A_R_B0", 
            //"Optim_CLG",                                   
        };

        Simulation_MSM_Population_BatchRun.main(grpDirNames);

    }

}
