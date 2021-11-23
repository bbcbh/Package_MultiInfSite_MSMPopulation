package test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import static sim.SinglePopRunnable.VACCINE_END_TIME;
import static sim.SinglePopRunnable.VACCINE_SETTING_LENGTH;
import static sim.SinglePopRunnable.VACCINE_START_TIME;
import util.PropValUtils;

/**
 *
 * @author Ben Hui
 */
public class Util_Population_Gen_Prop_Vaccine {

    public static final File FILE_BASEPROP = new File("C:\\Users\\bhui\\OneDrive - UNSW\\MSM_MulitSite\\NG Vacc\\Vacc_Blank");
    public static final File FILE_TARGET_DIR = new File("C:\\Users\\bhui\\OneDrive - UNSW\\MSM_MulitSite\\NG Vacc\\Vacc_Gen");
    private static final String PROP_FILE_NAME = "simSpecificSim.prop";
    private static final String KEY_MSM_VACCINE_SETTING = "PROP_MSM_VACCINE_SETTING";

    public static void main(String[] arg) throws
            FileNotFoundException, IOException, ParserConfigurationException,
            SAXException, TransformerException {

        File propFile = new File(FILE_BASEPROP, PROP_FILE_NAME);
        Document xml_src = PropValUtils.parseXMLFile(propFile);

        NodeList nList_entry = xml_src.getElementsByTagName("entry");
        Element src_vaccine = null;

        for (int entId = 0; entId < nList_entry.getLength(); entId++) {
            Element entryElement = (Element) nList_entry.item(entId);
            if (KEY_MSM_VACCINE_SETTING.equals(entryElement.getAttribute("key"))) {
                src_vaccine = entryElement;
            }
        }
        if (src_vaccine == null) {
            src_vaccine = xml_src.createElement("entry");
            src_vaccine.setAttribute("key", KEY_MSM_VACCINE_SETTING);
        }

        /*
        PROP_MSM_VACCINE_SETTING:
        
        0:VACCINE_SETTING_LENGTH
        1:VACCINE_START_TIME
        2:VACCINE_END_TIME -  re-vaccine booster time if end time is less than 0,
        
        3:EFFECT_INDEX_PROPORTION_VACC_COVERAGE_SETTING
        4:EFFECT_INDEX_TRANMISSION_EFFICACY_G
        5:EFFECT_INDEX_TRANMISSION_EFFICACY_A 
        6:EFFECT_INDEX_TRANMISSION_EFFICACY_R

        7:EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_G
        8:EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_A
        9:EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_R    
        
        10:OPTIONAL_EFFECT_VACCINE_DURATION_DEFAULT
        11:OPTIONAL_EFFECT_ADJ_INF_DUR_DEFAULT    
        
        12:OPTIONAL_EFFECT_REMOVE_SYM_RATE_DEFAULT
        13:OPTIONAL_EFFECT_REMOVE_SYM_STATE_DEFAULT
        14:OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_MEDIAN
        15:OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_SD
         */
        double[] COVERAGE_RANGE = new double[]{0.30};
        double[] SUSCEPTIBLE_ADJ = new double[]{0, 0.25, 0.50, 0.75, 1};
        double[] TRANMISSION_ADJ = new double[]{0, 0.25, 0.50, 0.75, 1};
        double[] VACC_DURATION = new double[]{2 * 360};

        double[] vaccineSettingBase = new double[16];
        Arrays.fill(vaccineSettingBase, -1);

        vaccineSettingBase[VACCINE_SETTING_LENGTH] = 13;
        vaccineSettingBase[VACCINE_START_TIME] = 7920;

        int offset = 3;
        int VACC_COVERAGE_SETTING
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_PROPORTION_VACC_COVERAGE_SETTING + offset;
        int TRANMISSION_EFFICACY_G
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_G + offset;
        int TRANMISSION_EFFICACY_A
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_A + offset;
        int TRANMISSION_EFFICACY_R
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_R + offset;
        int SUSCEPTIBLE_EFFICACY_G
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_G + offset;
        int SUSCEPTIBLE_EFFICACY_A
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_A + offset;
        int SUSCEPTIBLE_EFFICACY_R
                = infection.vaccination.SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_R + offset;
        int VACCINE_DURATION_DEFAULT
                = infection.vaccination.SiteSpecificVaccination.OPTIONAL_EFFECT_VACCINE_DURATION_DEFAULT + offset;

        int REMOVE_SYM_RATE_DEFAULT
                = infection.vaccination.SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_RATE_DEFAULT + offset;
        int REMOVE_SYM_STATE_DEFAULT
                = infection.vaccination.SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_STATE_DEFAULT + offset;
        int REMOVE_SYM_INF_DUR_DEFAULT_MEDIAN
                = infection.vaccination.SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_MEDIAN + offset;
        int REMOVE_SYM_INF_DUR_DEFAULT_SD
                = infection.vaccination.SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_SD + offset;

        vaccineSettingBase[VACC_COVERAGE_SETTING] = -0.3;
        vaccineSettingBase[TRANMISSION_EFFICACY_G] = 1;
        vaccineSettingBase[TRANMISSION_EFFICACY_A] = 1;
        vaccineSettingBase[TRANMISSION_EFFICACY_R] = 1;
        vaccineSettingBase[SUSCEPTIBLE_EFFICACY_G] = 1;
        vaccineSettingBase[SUSCEPTIBLE_EFFICACY_A] = 1;
        vaccineSettingBase[SUSCEPTIBLE_EFFICACY_R] = 1;
        vaccineSettingBase[VACCINE_DURATION_DEFAULT] = 2 * 360;

        // Generate prop
        double[] vaccineSetting;
        String folderName;

        // Coverage, duration, protect efficiacy map
        for (double coverage : COVERAGE_RANGE) {
            vaccineSetting = Arrays.copyOf(vaccineSettingBase, vaccineSettingBase.length);
            vaccineSetting[VACC_COVERAGE_SETTING] = -coverage;
            for (double sus : SUSCEPTIBLE_ADJ) {
                vaccineSetting[SUSCEPTIBLE_EFFICACY_G] = sus;
                vaccineSetting[SUSCEPTIBLE_EFFICACY_A] = sus;
                vaccineSetting[SUSCEPTIBLE_EFFICACY_R] = sus;

                for (double tran : TRANMISSION_ADJ) {

                    vaccineSetting[TRANMISSION_EFFICACY_G] = tran;
                    vaccineSetting[TRANMISSION_EFFICACY_A] = tran;
                    vaccineSetting[TRANMISSION_EFFICACY_R] = tran;

                    if (tran == 0 || sus != 0) {

                        for (double dur : VACC_DURATION) {
                            vaccineSetting[VACCINE_DURATION_DEFAULT] = dur;
                            // Without Booster
                            vaccineSetting[VACCINE_END_TIME] = vaccineSetting[VACCINE_START_TIME];

                            // With Booster                            
                            //vaccineSetting[VACCINE_END_TIME] = -(dur + 360); // Booster
                            
                            boolean st_map = !true;
                            
                            double R_EFF = 0.5;                            
                            boolean extra_r_eff = R_EFF != 0;
                            
                            boolean extra_sym_r = false;
                            File genPropFile;
                            
                            if(st_map){
                            
                            if (vaccineSettingBase[VACCINE_START_TIME] == vaccineSetting[VACCINE_END_TIME]) {
                                folderName = String.format("Vacc_C%03d_S%03d_T%03d_D%04d",
                                        (int) (coverage * 100),
                                        (int) (sus * 100),
                                        (int) (tran * 100),
                                        (int) dur);
                            } else {
                                folderName = String.format("Vacc_C%03d_S%03d_T%03d_D%04d_B%04d",
                                        (int) (coverage * 100),
                                        (int) (sus * 100),
                                        (int) (tran * 100),
                                        (int) dur,
                                        (int) -vaccineSetting[VACCINE_END_TIME]);
                            }
                            genPropFile = new File(FILE_TARGET_DIR, folderName);
                            genPropFile.mkdirs();
                            genPropFile = new File(genPropFile, PROP_FILE_NAME);

                            src_vaccine.setTextContent(Arrays.deepToString(new double[][]{vaccineSetting}));
                            PropValUtils.replacePropEntryByDOM(xml_src, genPropFile,
                                    new Element[]{src_vaccine}, null);
                            
                            }

                            

                            if (extra_r_eff) {
                                double[] vaccineSettingR = Arrays.copyOf(vaccineSetting, vaccineSetting.length);
                          
                                // divide as it is stored as tranmission adjustment
                                vaccineSettingR[SUSCEPTIBLE_EFFICACY_R] = Math.min(1, vaccineSettingR[SUSCEPTIBLE_EFFICACY_R]/R_EFF); 
                                vaccineSettingR[TRANMISSION_EFFICACY_R] = Math.min(1, vaccineSettingR[TRANMISSION_EFFICACY_R]/R_EFF);
                                
                            

                                folderName = String.format("Vacc_C%03d_S%03d_T%03d_D%04d_R%03d",
                                        (int) (coverage * 100),
                                        (int) (sus * 100),
                                        (int) (tran * 100),
                                        (int) dur,
                                        (int) (R_EFF*100));

                                genPropFile = new File(FILE_TARGET_DIR, folderName);
                                genPropFile.mkdirs();
                                genPropFile = new File(genPropFile, PROP_FILE_NAME);

                                src_vaccine.setTextContent(Arrays.deepToString(new double[][]{vaccineSettingR}));
                                PropValUtils.replacePropEntryByDOM(xml_src, genPropFile,
                                        new Element[]{src_vaccine}, null);

                            }
                            if (extra_sym_r) {

                                double[] vaccineSettingSym = Arrays.copyOf(vaccineSetting, vaccineSetting.length);
                                vaccineSettingSym[REMOVE_SYM_RATE_DEFAULT] = 1;
                                vaccineSettingSym[REMOVE_SYM_STATE_DEFAULT] = 1;
                                vaccineSettingSym[REMOVE_SYM_INF_DUR_DEFAULT_MEDIAN] = 185;
                                vaccineSettingSym[REMOVE_SYM_INF_DUR_DEFAULT_SD] = 35;

                                folderName = String.format("Vacc_C%03d_S%03d_T%03d_D%04d_SYMR100",
                                        (int) (coverage * 100),
                                        (int) (sus * 100),
                                        (int) (tran * 100),
                                        (int) dur);

                                genPropFile = new File(FILE_TARGET_DIR, folderName);
                                genPropFile.mkdirs();
                                genPropFile = new File(genPropFile, PROP_FILE_NAME);

                                src_vaccine.setTextContent(Arrays.deepToString(new double[][]{vaccineSettingSym}));
                                PropValUtils.replacePropEntryByDOM(xml_src, genPropFile,
                                        new Element[]{src_vaccine}, null);

                            }
                            
                        }
                    }

                }
            }
        }

    }

}
