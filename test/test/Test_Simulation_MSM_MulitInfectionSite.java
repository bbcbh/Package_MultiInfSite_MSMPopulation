package test;

import java.io.IOException;
import sim.Simulation_MSM_Population_BatchRun;

public class Test_Simulation_MSM_MulitInfectionSite {

    public static void main(String[] arg) throws IOException, InterruptedException, ClassNotFoundException {

   
        String[] grpDirNames = new String[]{
            "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test", 
            //"Import_1_BaseModel",
            //"Optim",               
            //"C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\1000 Runs",                
        };        
       
        Simulation_MSM_Population_BatchRun.main(grpDirNames);
        
        
    }

}
