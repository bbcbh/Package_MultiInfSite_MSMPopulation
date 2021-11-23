package population.person;

import infection.AbstractInfection;
import java.util.logging.Level;
import java.util.logging.Logger;
import person.AbstractIndividualInterface;

/**
 *
 * @author Ben Hui
 */
public class RelationshipPerson_MSM extends RelationshipPerson {

    /**
	 * 
	 */
	private static final long serialVersionUID = 3606660625193387151L;
	public static final int BEHAV_REG_ONLY = 0;
    public static final int BEHAV_CAS_ONLY = 1;
    public static final int BEHAV_BOTH = 2;
    public static final int SITE_G = 0;
    public static final int SITE_A = 1;
    public static final int SITE_R = 2;
    
    // e.g. if purely receptive to anal sex only then the person's 
    // genital cannot be infected through anal sex 
    // i.e. PARAM_IMMUNE_ACT_SITE_G = IMMUNE_ACT_ANAL
    // A person receptive to both anal and oral sex will be
    // PARAM_IMMUNE_ACT_SITE_G = IMMUNE_ACT_ANAL & IMMMUN_ACT_ORAL
    
    public static final int IMMUNE_ACT_ANAL = 0b1;  
    public static final int IMMMUN_ACT_ORAL = 0b10;        
    
    protected static final String[] PARA_NANE = {"PARAM_MSM_BEHAV_TYPE", "PARAM_NUM_CASUAL_IN_LAST_6_MONTHS", 
        "PARAM_IMMUNE_ACT_SITE_G", "PARAM_IMMUNE_ACT_SITE_A", "PARAM_IMMUNE_ACT_SITE_R",
        "PARAM_LAST_UNPROTECTED_ANAL_SEX_AT_AGE", "PARAM_LAST_UNPROTECTED_ORAL_SEX_AT_AGE",
        "PARAM_LAST_SCREEN_AT_AGE", "PARAM_SCHEDULED_SCREEN_AT_AGE"
            
    };
    public static final int PARAM_BEHAV_TYPE_INDEX = RelationshipPerson.LENGTH_PARAM_TOTAL; // 0 = Reg only, 1 = Casual only, 2 = Mixed     
    public static final int PARAM_NUM_CASUAL_IN_LAST_6_MONTHS = PARAM_BEHAV_TYPE_INDEX + 1;
    public static final int PARAM_IMMUNE_ACT_SITE_G = PARAM_NUM_CASUAL_IN_LAST_6_MONTHS + 1;
    public static final int PARAM_IMMUNE_ACT_SITE_A = PARAM_IMMUNE_ACT_SITE_G + 1;
    public static final int PARAM_IMMUNE_ACT_SITE_R = PARAM_IMMUNE_ACT_SITE_A + 1;
    public static final int PARAM_LAST_UNPROTECTED_ANAL_SEX_AT_AGE = PARAM_IMMUNE_ACT_SITE_R + 1;
    public static final int PARAM_LAST_UNPROTECTED_ORAL_SEX_AT_AGE = PARAM_LAST_UNPROTECTED_ANAL_SEX_AT_AGE + 1;
    public static final int PARAM_LAST_SCREEN_AT_AGE = PARAM_LAST_UNPROTECTED_ORAL_SEX_AT_AGE + 1;
    public static final int PARAM_SCHEDULED_SCREEN_AT_AGE = PARAM_LAST_SCREEN_AT_AGE + 1;
    
    protected int[] param = new int[PARA_NANE.length];
    // Casual encounter record 
    protected int[] casualRecord = new int[6 * 30];
    protected int casualRecordIndex = 0;
    // Logger
    private static final Logger LOG = Logger.getLogger(RelationshipPerson_MSM.class.getName());

    public RelationshipPerson_MSM(int id, boolean isMale, double age, int numInf) {
        super(id, isMale, age);
        this.initalisedInfections(numInf);  // 3 sites                
    }

    public void addCasualPartner(AbstractIndividualInterface p) {
        casualRecord[casualRecordIndex] = p.getId();       
        param[PARAM_NUM_CASUAL_IN_LAST_6_MONTHS-RelationshipPerson.LENGTH_PARAM_TOTAL]++;
    }

    @Override
    public int incrementTime(int deltaT, AbstractInfection[] infectionList) {
        int res = super.incrementTime(deltaT, infectionList);
        for (int t = 0; t < deltaT; t++) {
            casualRecordIndex = (casualRecordIndex + 1) % casualRecord.length;
            if (casualRecord[casualRecordIndex] != 0) {
                // Remove casual from record
                param[PARAM_NUM_CASUAL_IN_LAST_6_MONTHS - RelationshipPerson.LENGTH_PARAM_TOTAL]--;
                casualRecord[casualRecordIndex] = 0;
            }
        }
        return res;
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
    public Comparable getParameter(String id) {
        for (int i = 0; i < PARA_NANE.length; i++) {
            if (PARA_NANE[i].equals(id)) {
                return param[i];
            }
        }
        return super.getParameter(id);
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
    @Override
    public Comparable setParameter(String id, Comparable val) {
        //Comparable ret;
        for (int i = 0; i < PARA_NANE.length; i++) {
            if (PARA_NANE[i].equals(id)) {
                //ret = param[i];
                try {
                    param[i] = ((Number) val).intValue();
                } catch (ClassCastException ex) {
                    LOG.log(Level.WARNING,
                            "Parameter [" + val.toString() + "] cannot be cast to fit " + PARA_NANE[i] + "'s type",
                            ex);
                }

            }
        }
        return super.setParameter(id, val);
    }

    @Override
    public String indexToParamName(int index) {
        if (index < RelationshipPerson.LENGTH_PARAM_TOTAL) {
            return super.indexToParamName(index);
        } else {
            return PARA_NANE[index - LENGTH_PARAM_TOTAL];
        }

    }

    public int[] getCasualRecord() {
        return casualRecord;
    }    

    public void setCasualRecord(int[] casualRecord) {
        this.casualRecord = casualRecord;
    }    
    
    
}
