package test;

import java.io.IOException;
import sim.Simulation_MSM_Population_BatchRun;

public class Test_Simulation_MSM_MulitInfectionSite {

    public static void main(String[] arg) throws IOException, InterruptedException, ClassNotFoundException {

   
        String[] grpDirNames = new String[]{
            "C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test", 
            //"C:\\Users\\bbcbh\\Documents\\NetbeanProjects\\Test", 
            //"Import_1_BaseModel_Test",
            //"Import_1_No_Treatment_R",
            //"Import_1_50Eff_Treatment_R"
            "Import_1_No_Treatment_R_TSA_Opt7L",            
            "Optim_7",               
            //"C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\1000 Runs",                
        };        
       
        Simulation_MSM_Population_BatchRun.main(grpDirNames);
        
        
    }

}
