package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 *
 * @author Ben Hui
 */
public class Util_CopySimSpecificProp {

    public static final File FILE_BASEPROP = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\MSM_MulitSite\\Test\\BasePropFile\\HPC\\simSpecificSim.prop");

    public static final File FILE_TARGET_DIR = new File("C:\\Users\\Bhui\\Desktop\\FTP\\MSM");

    public static void main(String[] arg) throws FileNotFoundException, IOException {

        ArrayList<String> lines = new ArrayList();
        BufferedReader reader = new BufferedReader(new FileReader(FILE_BASEPROP));
        String line;

        int counter = 0;

        while ((line = reader.readLine()) != null) {
            lines.add(line);
            counter++;
        }

        System.out.println("Number of lines read = " + counter);

        // Vaccine offer protection 
        String replace_line_src  = "<entry key=\"PROP_MSM_VACCINE_SETTING\">[[8,7920,7920,-0.3,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,720]]</entry>";
        
        
        for (int sus = 0; sus < 100; sus += 10) {
            String dirName;
            for (int tran : new int[]{0, 25, 50, 75, 100}) {
                dirName = String.format("Vacc_C30_S%03d_T%03d", sus, tran);                
                File targetDirectory = new File(FILE_TARGET_DIR, dirName);
                
                int[] replace_lines_number = new int[]{109};
                String[] replace_lines = new String[]{String.format(
                        replace_line_src, 
                        tran/100.0, tran/100.0, tran/100.0,
                        sus/100.0, sus/100.0, sus/100.0)};                     
                genNewPropFile(targetDirectory, lines, replace_lines_number, replace_lines);
            }
        }                        
    }

    protected static void genNewPropFile(File targetDirectory, ArrayList<String> lines,
            int[] replace_lines_number, String[] replace_lines) throws FileNotFoundException {

        targetDirectory.mkdirs();
        PrintWriter pri = new PrintWriter(new File(targetDirectory, FILE_BASEPROP.getName()));
        int lineCounter = 0;
        int rPt = 0;

        for (String org_line : lines.toArray(new String[lines.size()])) {
            if (rPt < replace_lines_number.length 
                    && replace_lines_number[rPt] == lineCounter) {
                org_line = replace_lines[rPt];
                rPt++;
            }
            pri.println(org_line);
            lineCounter++;
        }

        pri.close();

        System.out.println("New propFile generated at " + targetDirectory.getAbsolutePath());
    }

}
