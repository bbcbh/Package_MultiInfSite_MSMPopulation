package test;

import java.io.IOException;
import sim.Simulation_MSM_Population_BatchRun;

public class Test_Simulation_MSM_MulitInfectionSite {

    public static void main(String[] arg) throws IOException, InterruptedException, ClassNotFoundException {

   
        String[] grpDirNames = new String[]{
            "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test",             
            //"Import_1_B",
            "Import_1_R",
            "Import_1_50R"                       
            //"Import_1_No_Treatment_R_B2",   
            //"Import_1_No_Treatment_R_B1",
            //"Import_1_No_Treatment_R_B0",
            //"Optim_CLG",               
            //"C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\1000 Runs",                
        };        
       
        Simulation_MSM_Population_BatchRun.main(grpDirNames);
        
        
    }

}
