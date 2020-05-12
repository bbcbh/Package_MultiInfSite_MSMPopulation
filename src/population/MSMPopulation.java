package population;

import availability.AbstractAvailability;
import infection.AbstractInfection;
import infection.GonorrhoeaSiteInfection;
import infection.vaccination.AbstractVaccination;
import infection.vaccination.SiteSpecificVaccination;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import person.AbstractIndividualInterface;
import static population.AbstractRegCasRelMapPopulation.*;

import population.availability.MSMAvailablity;

import population.person.RelationshipPerson_MSM;
import population.relationshipMap.RegCasRelationship;
import random.MersenneTwisterRandomGenerator;
import relationship.RelationshipMap;
import relationship.SingleRelationship;
import util.PersonClassifier;

import util.StaticMethods;
import population.person.MultiSiteMultiStrainPersonInterface;
import static infection.vaccination.SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_STATE_DEFAULT;
import java.util.ArrayList;
import sim.SinglePopRunnable;

/**
 * @author Ben Hui
 * @version 20190620
 *
 * <p>
 * History:</p>
 * <p>
 * 20140829 - Added treatPerson method</p>
 * <p>
 * 20140829 - Added support for re-screening and target screening</p>
 * <p>
 * 20140829 - Added support for screening </p>
 * <p>
 * 20140829 - Added support for insertive/receptive only partner and support function </p>
 * <p>
 * 20140827 - Additional index for genital infection for unprotected sex </p>
 * <p>
 * 20140807 - Add return value for perform act for unprotected sex (previously unused) </p>
 * <p>
 * 20140409 - Add support for Rimming (oral - anal tranmission) </p>
 * <p>
 * 20140328 - Add support for MSM_SYM_DURATION_OVERWRITE_PROB </p>
 * <p>
 * 20140212 - Add support for delayed condom usage adjustment MSM_REG_CONDOM_USAGE_ADJ, MSM_CAS_CONDOM_USAGE_ADJ </p>
 * <p>
 * 20140212 - Add support for delayed act adjustment MSM_REG_FREQ_ADJ, MSM_CAS_FREQ_ADJ </p>
 * <p>
 * 20140212 - Bug fix - Act frequency field are now fixed </p>
 * <p>
 * 20140211 - Add initial age structure </p>
 * <p>
 * 20140211 - Add alterative implementation of act frequency </p>
 * <p>
 * 20140211 - Added incident measure </p>
 * <p>
 * 20140130 - Added MSM_NUM_IN_GRP, MSM_CAS_PART_PROB and MSM_CAS_PART_SPREAD field </p>
 * <p>
 * 20140129 - Use product for transmission probability. Shared RNG for all infection </p>
 * <p>
 * 20130128 - Added getIndivdualInfectionList support</p>
 * <p>
 * 20180521 - Rework to fit in new package description</p>
 * <p>
 * 20180522 - Add treatment efficiency by site</p>
 * <p>
 * 20180528 - Add multiple strain support </p>
 * <p>
 * 20180920 - Add support for tranmission for kissing, slight change to the output definition of performAct to include tranmission from all site
 *
 * </p>
 * <p>
 * 20181011 - Add support for alterative format for freq. of act based on number of partner in last 12 months
 * </p>
 * <p>
 * 20190110 - Change the definition of symptomatic infection to self-sought treatment.
 * </p>
 * <p>
 * 20190111 - Adding input for % symptomatic by site.
 * </p>
 * <p>
 * 20190115 - Add support for population based distribution for casual partnership
 * </p>
 * <p>
 * 20190115 - Add alterate screening set format for target screening. Remove requirement on last unprotected anal sex
 * </p>
 * <p>
 * 20190523 - Add support for vaccination
 * </p>
 * <p>
 * 20190620 - Add support for adjusting population size.
 * </p>
 */
public class MSMPopulation extends AbstractRegCasRelMapPopulation {

    public int numInPop = 10000;

    public void setInitNumInPop(int numInPop) {

        if (getPop() != null) {
            System.err.println("Number in population already set as " + getPop().length
                    + ". New population size not setted.");

        } else {
            this.numInPop = numInPop;
        }
    }
    public static final int ACT_ANAL = 0;  // Coresponds to length of MSM_POP_REG_ACT_FREQ and MSM_POP_CAS_ACT_FREQ
    public static final int ACT_ORAL = 1;
    public static final int ACT_RIMMING = 2;
    public static final int ACT_KISSING = 3;

    public static final int TRAN_SUSC_INDEX_RIMMING_ANAL = 3; // After G,A,R
    public static final int TRAN_SUSC_INDEX_RIMMING_ORAL = TRAN_SUSC_INDEX_RIMMING_ANAL + 1;
    public static final int TRAN_SUSC_INDEX_KISSING = TRAN_SUSC_INDEX_RIMMING_ORAL + 1;

    public static final int MSM_POP_INDEX_OFFSET = AbstractRegCasRelMapPopulation.LENGTH_FIELDS;
    public static final int MSM_POP_REG_ACT_FREQ = MSM_POP_INDEX_OFFSET;
    public static final int MSM_POP_CAS_ACT_FREQ = MSM_POP_REG_ACT_FREQ + 1;
    public static final int MSM_POP_REG_CONDOM_USAGE = MSM_POP_CAS_ACT_FREQ + 1;
    public static final int MSM_POP_CAS_CONDOM_USAGE = MSM_POP_REG_CONDOM_USAGE + 1;
    public static final int MSM_REG_LENTH_AVE = MSM_POP_CAS_CONDOM_USAGE + 1;
    public static final int MSM_NUM_IN_GRP = MSM_REG_LENTH_AVE + 1;
    public static final int MSM_CAS_PART_PROB = MSM_NUM_IN_GRP + 1;
    public static final int MSM_CAS_PART_SPREAD = MSM_CAS_PART_PROB + 1;
    public static final int MSM_REG_FREQ_ADJ = MSM_CAS_PART_SPREAD + 1;
    public static final int MSM_CAS_FREQ_ADJ = MSM_REG_FREQ_ADJ + 1;
    public static final int MSM_REG_CONDOM_USAGE_ADJ = MSM_CAS_FREQ_ADJ + 1;
    public static final int MSM_CAS_CONDOM_USAGE_ADJ = MSM_REG_CONDOM_USAGE_ADJ + 1;
    public static final int MSM_SYM_DURATION_OVERWRITE_PROB = MSM_CAS_CONDOM_USAGE_ADJ + 1;
    public static final int MSM_INSERTIVE_RECEPTIVE_ACCUM_DIST = MSM_SYM_DURATION_OVERWRITE_PROB + 1;
    public static final int MSM_SCREENING_SETTING_GENERAL = MSM_INSERTIVE_RECEPTIVE_ACCUM_DIST + 1;
    public static final int MSM_SCREENING_SETTING_TARGETED = MSM_SCREENING_SETTING_GENERAL + 1;
    public static final int MSM_TREATMENT_EFFICIENCY = MSM_SCREENING_SETTING_TARGETED + 1;
    public static final int MSM_SYM_TREATMENT_PROB = MSM_TREATMENT_EFFICIENCY + 1;
    public static final int MSM_USE_GLOBAL_CASUAL_LIMIT = MSM_SYM_TREATMENT_PROB + 1;
    public static final int MSM_SITE_SPECIFIC_VACCINATION = MSM_USE_GLOBAL_CASUAL_LIMIT + 1;
    public static final int MSM_SITE_CURRENTLY_VACCINATED = MSM_SITE_SPECIFIC_VACCINATION + 1;
    public static final int MSM_SITE_VACC_BOOSTER_SCHEDULE = MSM_SITE_CURRENTLY_VACCINATED + 1;
    public static final int LENGTH_FIELDS_MSM_POP = MSM_SITE_VACC_BOOSTER_SCHEDULE + 1;

    public static final Object[] DEFAULT_MSM_FIELDS = {
        // Min/max for anal and oral sex for reg (per day) and causal relationship (per partnership), from 
        // Crawford et al. Number of risk acts by relationship status
        //  and partner serostatus: Findings from the HIM cohort of homosexually active men in Sydney,
        //   Australia. AIDS and behavior. 2006;10(3):325-31.
        // Alterative format : {[min, max, prob...], [min, max, prob...],}             
        // Alterative format: {[-1, number partners limit over 6 months ], 
        //                      [prob of anal under limit, prob of anal over limit],
        //                      [prob of oral under limit, prob of oral over limit]... }         
        new float[][]{{0.2286f, 0.3429f}, {0.2286f, 0.3429f}}, // {{1.6f / 7, 2.4f / 7}, {1.6f / 7, 2.4f / 7}},
        new float[][]{{1, 8, 0.415f}, {1, 10, 0.825f}},
        // Condom usage by relation type and act
        // KI Sureillance report 2013 table 5.1.1
        // Unprotected anal intercourse at Mel (R,C) =  (28.6%, 23.3% ) 2013 (34.8% 26.3%) 2012
        new float[]{1 - 0.286f, 0},
        new float[]{1 - 0.233f, 0},
        // MSM_REG_LENTH_AVE
        (double) 4 * 360,
        // GCPS Melbourne 2013: Reg partnership 33.0% only regular, 24.2 % only causal, 26.9% has both with 15.9 has no sex
        new float[]{0.33f / (0.33f + 0.242f + 0.269f), 0.242f / (0.33f + 0.242f + 0.269f), 0.269f / (0.33f + 0.242f + 0.269f)},
        /*
         * Fogarty A, et al. viour of HIV-negative and HIV-positive gay men, 2002–2005,2006.
         *
         * 1-1.5 26% 2-5 21% 6-10 16% 11-50 30% 51-60 7%
         *
         * For every 6 months
         */
        new float[]{0.26f, 0.21f, 0.16f, 0.30f, 0.07f},
        new int[][]{{1, 2}, {2, 5}, {6, 10}, {11, 50}, {51, 60}},
        // Regular act adjustment Format:[globalTime actType_adj...]
        new float[]{},
        // Casual act adjustment  Format: [globalTime actType_adj...]
        new float[]{},
        // Regular protection usage adjustment Format:[globalTime actType_adj...]
        new float[]{},
        // Casual protection usage adjustment  Format: [globalTime actType_adj...]
        new float[]{},
        // Symptom overwrite probability - probability of treatment can reduce duration on other sites
        // upon appearance of symtoms
        new float[]{1f, 1f, 1f},
        // Insertive/ receptive by acts, by cummulative probability       
        // Format: float{insertive only to all, insertive only to anal, insertive only to oral, 
        //               receptive only to all, receptive only to anal, receptive only to oral, }
        // i.e. index = 
        // 0b00 = all, 0b01 = 1<< ACT_ANAL, 0b10 = 1<< ACT_ORAL, 0b100 = 1<< ACT_RIMMING (if implemented)
        new float[]{0f, 0f, 0f, 0f, 0f, 0f},
        // General screening setting 
        // Format: float[]{startTime, totalScreened, frequency, 
        //    cumulative probability of screening at site (2^3 -1) in total}      
        new float[]{Float.POSITIVE_INFINITY, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        // Target screening
        // Format: float[]{startTime, 
        //   prob of targeted screen per day (e.g. 70% within 3 months := (1- 0.7) =  (1-p)^(3*30), or p = 0.0133)              
        //   prob of repeat screen after 3 months         
        new float[]{Float.POSITIVE_INFINITY, 0, 0},
        // Treatment efficiency by site
        // Format: float[site]{eff_for_s1, eff_for_s2 ...}
        new float[][]{{1, 1, 1}},
        // Probablility of syptomatic/sought treatment at each site
        // Urethral: Assumption
        // Anal : Bissessor 2011: 7 out 47 of has proctitis
        // Phary : Assumption
        new float[]{1, 7 / 4.7f, 0},
        // Global casual limit        
        true,
        // MSM_SITE_SPECIFIC_VACCINATION
        // A SiteSpecificVaccination object if applied
        null,
        // MSM_SITE_CURRENTLY_VACCINATED
        // HashMap<Integer, int[]> currentlyVaccinated of ID, VACC_SETTING if applied
        new HashMap<Integer, int[]>(),
        // MSM_SITE_VACC_BOOSTER_SCHEDULE
        // HashMap<Integer, ArrayList<Integer>> global_time, list of id who receive booseter 
        new HashMap<Integer, ArrayList<Integer>>(),};

    public static final int[] AGE_RANGE = {(int) (16 * AbstractRegCasRelMapPopulation.ONE_YEAR_INT),
        (int) (80 * AbstractRegCasRelMapPopulation.ONE_YEAR_INT)
    };
    public static final int MAPPING_REG = 0;
    public static final int MAPPING_CAS = 1;

    // Pop Snap Count (Acumulative)    
    private transient final int[] cumulativeIncidencesBySites = new int[3]; // 3 Sites, in accordance to Relationship_MSM.Site_X

    // Global casual limit    
    private transient float[] CASUAL_PARTNER_PROB;

    // Screening (General)
    private transient int[] screeningPersonByIndex = null; // Length of total screened
    private transient int[] screeningSite = null; // Length of total screened
    private transient int screenTarPt = 0;
    private transient int screenDayPt = 0;
    private transient int[] screeningToday = null; // Length of screening freq

    // Vaccine    
    public static final int VACC_SETTING_AGE_EXPIRY = 0;
    public static final int VACC_SETTING_LENGTH = VACC_SETTING_AGE_EXPIRY + 1;

    private transient AbstractRealDistribution vaccine_duration_dist = null;
    private transient AbstractRealDistribution vaccine_remove_sym_infect_duration = null;

    // Screening (Targeted)
    private PersonClassifier targetedScreenClassifier = new PersonClassifier() {
        @Override
        public int classifyPerson(AbstractIndividualInterface p) {
            if (p instanceof RelationshipPerson_MSM) {
                // From STIGMA guideline
                //int lastUnprotectAnalSex = (int) p.getParameter("PARAM_LAST_UNPROTECTED_ANAL_SEX_AT_AGE");
                //int lastScreenAtAge = (int) p.getParameter("PARAM_LAST_SCREEN_AT_AGE");
                //int numCasualPartnerInLast6Months = (int) p.getParameter("PARAM_NUM_CASUAL_IN_LAST_6_MONTHS");

                int[] casualRec = ((RelationshipPerson_MSM) p).getCasualRecord();
                int numCasual = 0;
                for (int i = 0; i < casualRec.length; i++) {
                    numCasual += (casualRec[i] != 0) ? 1 : 0;
                }

                return (numCasual >= 10) ? 0 : -1;

            } else {
                return -1;
            }

        }

        @Override
        public int numClass() {
            return 1;
        }
    };

    public MSMPopulation(long seed) {
        setSeed(seed);
        setRNG(new random.MersenneTwisterRandomGenerator(seed));

        Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_MSM_POP);
        for (int i = MSM_POP_INDEX_OFFSET; i < newFields.length; i++) {
            newFields[i] = DEFAULT_MSM_FIELDS[i - MSM_POP_INDEX_OFFSET];
        }

        // Mostly set for class definition.
        // Alterative - { g, a, r, tranmission for rimming (a), tranmission for rimming (r)},
        newFields[FIELDS_TRANSMIT] = new double[][]{{1, 0.0}, {0.024258237711416882, 0.0}, {0.08654020531624781, 0.0}};
        newFields[FIELDS_SUSCEPT] = new double[][]{{1, 0.0}, {0.8402332453833801, 0.0}, {0.6288277879484502, 0.0}};

        super.setFields(newFields);

    }

    public static MSMPopulation importMSMPopulation(java.io.ObjectInputStream objStr)
            throws IOException, ClassNotFoundException {
        Object[] importFields = (Object[]) objStr.readObject();
        Long seed = (Long) importFields[MSMPopulation.FIELDS_SEED];
        MSMPopulation newPop = new MSMPopulation(seed);
        newPop.setFields(importFields);
        return newPop;
    }

    public void setTargetedScreenClassifier(PersonClassifier targetedScreenClassifier) {
        this.targetedScreenClassifier = targetedScreenClassifier;
    }

    public int[] cumulativeIncidencesBySitesCount() {
        return Arrays.copyOf(cumulativeIncidencesBySites, cumulativeIncidencesBySites.length);
    }

    @Override
    public void initialise() {
        initialiseInfection(getSeed());
        /*
         * GCPS Melbourne 2013: Reg partnership 33.0% only regular, 24.2 % only causal, 26.9% has both with 15.9 has no sex
         */
        float[] numByGrp = (float[]) getFields()[MSM_NUM_IN_GRP];

        int[] numByGroupAccum = StaticMethods.accumulativeArray(numInPop, numByGrp);

        float[] casualPartProb = (float[]) getFields()[MSM_CAS_PART_PROB];
        int[][] casualPartLimit = (int[][]) getFields()[MSM_CAS_PART_SPREAD];
        float[] casualPartProbAccum = StaticMethods.accumulativeArray(casualPartProb);

        define_CASUAL_PARTNER_PROB(casualPartLimit, casualPartProb);

        setPop(new RelationshipPerson_MSM[numInPop]);
        setRelMap(new RelationshipMap[]{
            new RelationshipMap(),
            new RelationshipMap(),});
        setAvailablity(new AbstractAvailability[]{
            new MSMAvailablity(getRNG()),
            new MSMAvailablity(getRNG()),});

        int grpIndex = 0;

        for (int p = 0; p < numInPop; p++) {
            RelationshipPerson_MSM person = new RelationshipPerson_MSM(p + 1, true,
                    AGE_RANGE[0] + getRNG().nextInt(AGE_RANGE[1] - AGE_RANGE[0]), 3); // 3 sites
            getPop()[p] = person;
            person.setEnterPopulationAt(getGlobalTime());

            while (p >= numByGroupAccum[grpIndex] && (grpIndex + 1) < numByGroupAccum.length) {
                grpIndex++;
            }

            person.setParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX), grpIndex);

            if (grpIndex != RelationshipPerson_MSM.BEHAV_REG_ONLY) {
                float pC = getRNG().nextFloat();
                int f = 0;
                while (pC >= casualPartProbAccum[f] && (f + 1) < casualPartProbAccum.length) {
                    f++;
                }
                int numCas = casualPartLimit[f][0];
                if (casualPartLimit[f][1] - casualPartLimit[f][0] > 0) {
                    numCas += getRNG().nextInt(casualPartLimit[f][1] - casualPartLimit[f][0]);
                }

                //if(this instanceof HetroPopulation_MSMBehaviour){                    
                person.setTimeUntilNextRelationship(getRNG().nextInt(1 * 360)); // One offset for sexual debut                  
                //}
                person.setParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_MAX_PARTNER), Math.max(numCas, 1));

            }

            int behavType = ((Number) person.getParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue();

            if (person.getTimeUntilNextRelationship() < 0) {
                if (behavType != RelationshipPerson_MSM.BEHAV_CAS_ONLY) { // i.e. can have regular relationship
                    getRelMap()[MAPPING_REG].addAvailablePerson(person);
                }
                if (behavType != RelationshipPerson_MSM.BEHAV_REG_ONLY) { // i.e. can have casual relationship
                    getRelMap()[MAPPING_CAS].addAvailablePerson(person);
                }
            }

            // Check for insert
            String paramStrSus = AbstractInfection.PARAM_DIST_NEXT_VALUE_REGEX.replaceFirst("\\d+",
                    Integer.toString(GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX));
            for (int siteId = 0; siteId < person.getInfectionStatus().length; siteId++) {
                person.getProbSusBySite()[siteId] = (double) getInfList()[siteId].getParameter(paramStrSus);
            }

            if (getFields()[MSM_INSERTIVE_RECEPTIVE_ACCUM_DIST] != null) {
                float[] distByType = (float[]) getFields()[MSM_INSERTIVE_RECEPTIVE_ACCUM_DIST];

                if (distByType.length > 0) {
                    float pType = getRNG().nextFloat();
                    int pt = 0;

                    while (pt < distByType.length && pType >= distByType[pt]) {
                        pt++;
                    }

                    if (pt < distByType.length) { // insertive or receptive                          
                        boolean isReceptive = !(pt < distByType.length / 2);

                        // = 2 without rimming, or 3 with rimming
                        int actIncluded = (distByType.length / 2) - 1; // Not include all                          

                        // 0b01 = 1<< ACT_ANAL, 0b10 = 1<< ACT_ORAL, 
                        // 0b100 = 1<< ACT_RIMMING (if implemented)
                        int toActType = pt % (distByType.length / 2);
                        if (toActType == 0) { // Special case for all
                            toActType = (1 << actIncluded) - 1;
                        }

                        for (int actId = 0; actId < actIncluded; actId++) {
                            String immuneTarSite;
                            if (isReceptive) {
                                immuneTarSite = actId == ACT_RIMMING
                                        ? "PARAM_IMMUNE_ACT_SITE_R" : "PARAM_IMMUNE_ACT_SITE_G";
                            } else {
                                // Insertive
                                switch (actId) {
                                    case ACT_ANAL:
                                    case ACT_RIMMING:
                                        immuneTarSite = "PARAM_IMMUNE_ACT_SITE_A";
                                        break;
                                    case ACT_ORAL:
                                        immuneTarSite = "PARAM_IMMUNE_ACT_SITE_R";
                                        break;
                                    default:
                                        throw new UnsupportedOperationException("Insertive/receptive not defined");
                                }
                            }

                            if ((toActType & (1 << actId)) != 0) {
                                int immuneAct = 1 << actId;
                                int currentImmune = (int) person.getParameter(immuneTarSite);
                                currentImmune += immuneAct;
                                person.setParameter(immuneTarSite, currentImmune);
                            }
                        }

                    }
                }
            }

            snapshotCount(person);

        }

        getFields()[AbstractRegCasRelMapPopulation.FIELDS_NEXT_ID] = new Integer(getPop().length + 1);
        updatePairs();

    }

    protected void define_CASUAL_PARTNER_PROB(int[][] casualPartLimit, float[] casualPartProb) {
        CASUAL_PARTNER_PROB = new float[casualPartLimit[casualPartLimit.length - 1][1] + 1];

        for (int i = 0; i < casualPartLimit.length; i++) {
            int[] limit = Arrays.copyOf(casualPartLimit[i], casualPartLimit[i].length);

            if (i == 0) {
                limit[0] = 0;
            }

            if (i + 1 < casualPartLimit.length && limit[1] == casualPartLimit[i + 1][0]) {
                limit[1]--;
            }

            for (int k = limit[0]; k <= limit[1]; k++) {
                CASUAL_PARTNER_PROB[k] = casualPartProb[i] / (limit[1] - limit[0] + 1);
            }
        }
    }

    @Override
    public void advanceTimeStep(int deltaT) {

        incrementTime(deltaT);

        float[] scrSetting = (float[]) getFields()[MSM_SCREENING_SETTING_GENERAL];
        int[] todayScreenTarget = new int[0];
        int[] todayScreenSite = new int[0];
        int todayScreenPt = 0;

        if (getGlobalTime() >= scrSetting[0]) {
            defineScreeningSchedule(scrSetting);
        }

        if (CASUAL_PARTNER_PROB == null) {
            float[] casualPartProb = (float[]) getFields()[MSM_CAS_PART_PROB];
            int[][] casualPartLimit = (int[][]) getFields()[MSM_CAS_PART_SPREAD];
            define_CASUAL_PARTNER_PROB(casualPartLimit, casualPartProb);
        }

        if (screeningToday != null && screenDayPt < screeningToday.length) {
            todayScreenTarget = Arrays.copyOfRange(screeningPersonByIndex, screenTarPt,
                    screenTarPt + screeningToday[screenDayPt]);
            todayScreenSite = Arrays.copyOfRange(screeningSite, screenTarPt,
                    screenTarPt + screeningToday[screenDayPt]);
            screenTarPt += screeningToday[screenDayPt];
            screenDayPt++;
        }

        HashSet<RelationshipPerson_MSM> toBeRemoved = new HashSet<>();

        // Casual partnership stat 
        boolean useGlobalPopLimit = (Boolean) getFields()[MSM_USE_GLOBAL_CASUAL_LIMIT];
        int[] casualPart_Stat = new int[CASUAL_PARTNER_PROB.length];
        RelationshipPerson_MSM[][] casualPart_candidateCollection = new RelationshipPerson_MSM[CASUAL_PARTNER_PROB.length][getPop().length];
        int[] casualPart_candidateCollectionPt = new int[CASUAL_PARTNER_PROB.length];
        int casualPart_Total = 0;
        int regOnlyTotal = 0;

        AbstractVaccination vacc = ((AbstractVaccination) getFields()[MSM_SITE_SPECIFIC_VACCINATION]);

        for (int p = 0; p < getPop().length; p++) {
            RelationshipPerson_MSM person = (RelationshipPerson_MSM) getPop()[p];
            incrementPersonStat(person, deltaT);

            if (person.getAge() > AGE_RANGE[1]) {
                toBeRemoved.add(person);

                for (RelationshipMap relMap : getRelMap()) {
                    relMap.removeAvailablePerson(person);
                }

                // Remove person from vaccine record
                if (getFields()[MSM_SITE_SPECIFIC_VACCINATION] != null) {
                    ((AbstractVaccination) getFields()[MSM_SITE_SPECIFIC_VACCINATION]).getVaccinationRecord().remove(person.getId());
                }
                if (getFields()[MSM_SITE_CURRENTLY_VACCINATED] != null) {
                    ((HashMap<Integer, int[]>) getFields()[MSM_SITE_CURRENTLY_VACCINATED]).remove(person.getId());
                }

                int nextId = ((Number) getFields()[FIELDS_NEXT_ID]).intValue();

                RelationshipPerson_MSM newP = (RelationshipPerson_MSM) generateNewPerson(nextId, person, AGE_RANGE[0]);

                getLocalData().put(p, newP);
                person = newP;
                getFields()[FIELDS_NEXT_ID] = nextId + 1;

                if (vacc != null) {

                    if (vacc instanceof SiteSpecificVaccination) {
                        SiteSpecificVaccination ssv = (SiteSpecificVaccination) vacc;
                        double propVacc = ssv.getParameters()[SiteSpecificVaccination.EFFECT_INDEX_PROPORTION_VACC_COVERAGE_SETTING];
                        if (propVacc > 0) {
                            if (propVacc < 1) {
                                propVacc = getRNG().nextDouble() < propVacc ? 0 : 1;
                            }
                            if (propVacc >= 1) {
                                vaccinatePerson(vacc, person);
                            }
                        }
                    } else {
                        System.err.println("Vaccination format for " + vacc.getClass() + " not defined");
                    }

                }

            }

            int behaviour = ((Number) person.getParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX))).intValue();

            if (useGlobalPopLimit) {
                if (behaviour != RelationshipPerson_MSM.BEHAV_REG_ONLY) {
                    int[] casualRec = person.getCasualRecord();
                    int numCasual = 0;
                    for (int i = 0; i < casualRec.length; i++) {
                        numCasual += (casualRec[i] != 0) ? 1 : 0;
                    }

                    numCasual = Math.min(numCasual, casualPart_Stat.length - 1);
                    casualPart_Stat[numCasual]++;

                    if (person.getTimeUntilNextRelationship() < 0) {
                        casualPart_candidateCollection[numCasual][casualPart_candidateCollectionPt[numCasual]] = person;
                        casualPart_candidateCollectionPt[numCasual]++;
                    }
                    casualPart_Total++;
                } else {
                    regOnlyTotal++;
                }

            }

            if (person.getTimeUntilNextRelationship() < 0) {

                if (behaviour != RelationshipPerson_MSM.BEHAV_CAS_ONLY) { // i.e. can have regular partnership                
                    if (!getRelMap()[MAPPING_REG].containsVertex(person.getId())
                            || getRelMap()[MAPPING_REG].degreeOf(person.getId()) < 1) {

                        getRelMap()[MAPPING_REG].addAvailablePerson(person);
                    } else {
                        getRelMap()[MAPPING_REG].removeAvailablePerson(person);
                    }
                }
                if (!useGlobalPopLimit && behaviour != RelationshipPerson_MSM.BEHAV_REG_ONLY) { // i.e. can have casual 
                    boolean seekingCasual = person.getParameter(
                            person.indexToParamName(RelationshipPerson_MSM.PARAM_NUM_CASUAL_IN_LAST_6_MONTHS)).compareTo(
                            person.getParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_MAX_PARTNER))) < 0;

                    seekingCasual &= !getRelMap()[MAPPING_CAS].containsVertex(person.getId())
                            || getRelMap()[MAPPING_CAS].degreeOf(person.getId()) < 1;

                    if (seekingCasual) {
                        getRelMap()[MAPPING_CAS].addAvailablePerson(person);
                    } else {
                        getRelMap()[MAPPING_CAS].removeAvailablePerson(person);
                    }
                }
            }

            boolean screenedToday = false;
            if (todayScreenPt < todayScreenTarget.length
                    && p == todayScreenTarget[todayScreenPt]) {
                screenPerson((RelationshipPerson_MSM) person, todayScreenSite[todayScreenPt]);
                todayScreenPt++;
                screenedToday = true;
            }

            if (!screenedToday && getGlobalTime() >= ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[0]) {

                // Rescreening 
                screenedToday = person.getParameter("PARAM_SCHEDULED_SCREEN_AT_AGE").equals(person.getAge());

                int screenTypeStartIndex = 3;

                if (!screenedToday) {
                    // Targeted screening                    

                    if (targetedScreenClassifier.classifyPerson(person) >= 0
                            && ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[1] > 0) {

                        float srnTodayProb = ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[1];

                        if (srnTodayProb > 1) {
                            // Multi class option
                            int numClass = (int) srnTodayProb;
                            screenTypeStartIndex += numClass * 2;
                            float classFloat = getRNG().nextFloat();

                            int pt = Arrays.binarySearch((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED], 2, 2 + numClass, classFloat);

                            if (pt <= 0) {
                                pt = -(pt + 1);
                            }

                            if (pt < 2 + numClass) {
                                int probPt = pt + numClass;

                                srnTodayProb = ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[probPt];
                            } else {
                                srnTodayProb = 0;
                            }

                        }

                        if (srnTodayProb > 0) {
                            screenedToday = getRNG().nextFloat() < srnTodayProb;
                        }

                    }
                }

                if (screenedToday) {
                    int screeningType = 0b111;

                    if (((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED]).length > screenTypeStartIndex) {
                        screeningType = 0;
                        float pScrType = getRNG().nextFloat();
                        while (pScrType >= ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[screeningType + screenTypeStartIndex]) {
                            screeningType++;
                        }
                    }
                    screenPerson((RelationshipPerson_MSM) person, screeningType + 1);
                }
            }

            snapshotCount(person);
        }

        if (useGlobalPopLimit) {
            int[] casualPart_Target = new int[CASUAL_PARTNER_PROB.length];
            int targetSum = 0;
            for (int i = 0; i < casualPart_Target.length; i++) {
                casualPart_Target[i] = (int) Math.round(CASUAL_PARTNER_PROB[i] * casualPart_Total);
                targetSum += casualPart_Target[i];
            }
            // Remove randomly until exact match
            while (targetSum != casualPart_Total) {
                if (targetSum > casualPart_Total) {
                    casualPart_Target[getRNG().nextInt(casualPart_Target.length)]--;
                    targetSum--;
                } else {
                    casualPart_Target[getRNG().nextInt(casualPart_Target.length)]++;
                    targetSum++;
                }
            }

            int diff;

            //System.out.println(getGlobalTime() + ":" + Arrays.toString(casualPart_Stat));            
            for (int g = 0; g < casualPart_Stat.length; g++) {

                RelationshipPerson_MSM[] changeCasualBehaviorCandidates;

                int excessAbove = 0; // Need to form new partnership
                int excessBelow = 0;

                for (int k = 0; k < g; k++) {
                    excessBelow += Math.max(casualPart_Target[k] - casualPart_Stat[k], 0);
                }
                for (int k = g + 1; k < casualPart_Target.length; k++) {
                    excessAbove += Math.max(casualPart_Target[k] - casualPart_Stat[k], 0);
                }

                diff = Math.max(casualPart_Stat[g] - casualPart_Target[g], 0);

                changeCasualBehaviorCandidates = new RelationshipPerson_MSM[Math.min(diff, casualPart_candidateCollectionPt[g])];

                int cPt = 0;
                for (int c = 0; c < casualPart_candidateCollectionPt[g]; c++) {
                    if (getRNG().nextInt(casualPart_candidateCollectionPt[g] - c) < diff) {
                        changeCasualBehaviorCandidates[cPt] = casualPart_candidateCollection[g][c];
                        cPt++;
                        diff--;
                    } else {
                        getRelMap()[MAPPING_CAS].removeAvailablePerson(casualPart_candidateCollection[g][c]);
                    }
                }

                if (changeCasualBehaviorCandidates.length > 0) {

                    for (RelationshipPerson_MSM candidate : changeCasualBehaviorCandidates) {
                        if (getRNG().nextInt(excessAbove + excessBelow) < excessAbove) {
                            // Can seek new partner 
                            if (getRelMap()[MAPPING_CAS].containsVertex(candidate.getId())
                                    && getRelMap()[MAPPING_CAS].degreeOf(candidate.getId()) > 1) {
                                java.util.Set<SingleRelationship> edgeSet = getRelMap()[MAPPING_CAS].edgesOf(candidate.getId());
                                SingleRelationship[] existingCasual
                                        = edgeSet.toArray(new SingleRelationship[edgeSet.size()]);
                                getRelMap()[MAPPING_CAS].removeEdge(existingCasual[0]);
                            }
                            getRelMap()[MAPPING_CAS].addAvailablePerson(candidate);
                            excessAbove--;
                        } else {
                            getRelMap()[MAPPING_CAS].removeAvailablePerson(candidate);
                            candidate.setTimeUntilNextRelationship(getRNG().nextInt(12 * 30));
                            excessBelow--;
                        }
                    }

                }

            }

        }

        HashMap<Integer, ArrayList<Integer>> booster
                = (HashMap<Integer, ArrayList<Integer>>) getFields()[MSM_SITE_VACC_BOOSTER_SCHEDULE];

        ArrayList<Integer> boosterSchedule = booster.remove(getGlobalTime());
        for (int boosterId : boosterSchedule) {
            vaccinatePerson(vacc, getLocalData().get(boosterId));
        }

        for (RelationshipPerson_MSM removePerson : toBeRemoved) {
            for (RelationshipMap relMap : getRelMap()) {
                if (relMap.containsVertex(removePerson.getId())) {
                    relMap.removeVertex(removePerson);
                }
            }
        }

        updatePairs();

    }

    protected void vaccinatePerson(AbstractVaccination vaccine, AbstractIndividualInterface person) {
        if (vaccine != null) {
            if (vaccine.vaccinationApplicableAt(getGlobalTime())) {
                if (getFields()[MSM_SITE_CURRENTLY_VACCINATED] == null) {
                    HashMap<Integer, int[]> currentVaccinated = new HashMap<>();
                    getFields()[MSM_SITE_CURRENTLY_VACCINATED] = currentVaccinated;
                }

                vaccine.vaccinatePerson(person);

                // Test booster
                int[] vTime = vaccine.getValidTime();
                if (vTime[SinglePopRunnable.VACCINE_END_TIME] < 0) {
                    HashMap<Integer, ArrayList<Integer>> booster
                            = (HashMap<Integer, ArrayList<Integer>>) getFields()[MSM_SITE_VACC_BOOSTER_SCHEDULE];
                    int boosterTime = getGlobalTime() + -vTime[SinglePopRunnable.VACCINE_END_TIME];

                    ArrayList<Integer> ent = booster.get(boosterTime);

                    if (ent == null) {
                        ent = new ArrayList<>();
                        booster.put(boosterTime, ent);
                    }
                    ent.add(person.getId());
                }

                int[] vacSetting = new int[VACC_SETTING_LENGTH];
                vacSetting[VACC_SETTING_AGE_EXPIRY] = -1; // -1 = lifelong;

                if (vaccine instanceof SiteSpecificVaccination) {
                    if (vaccine.getParameters().length > SiteSpecificVaccination.OPTIONAL_EFFECT_VACCINE_DURATION_DEFAULT) {
                        double defaultDuration = ((SiteSpecificVaccination) vaccine).getParameters()[SiteSpecificVaccination.OPTIONAL_EFFECT_VACCINE_DURATION_DEFAULT];

                        if (defaultDuration > 0) {
                            if (vaccine_duration_dist == null) {
                                vaccine_duration_dist = new ExponentialDistribution(getRNG(), defaultDuration);
                            }
                            vacSetting[VACC_SETTING_AGE_EXPIRY] = (int) Math.round(person.getAge() + vaccine_duration_dist.sample());
                        }
                    }

                    // Immediate vaccine effect 
                    for (int site = 0; site < getInfList().length; site++) {
                        if (getInfList()[site].isInfectious(person)) {
                            vaccine_effect_adj_infection_duration((SiteSpecificVaccination) vaccine, person, site);
                        }
                        if (getInfList()[site].hasSymptoms(person)) {
                            vaccine_effect_remove_symptom((SiteSpecificVaccination) vaccine, person, site);
                        }
                    }

                }
                ((HashMap<Integer, int[]>) getFields()[MSM_SITE_CURRENTLY_VACCINATED]).put(person.getId(), vacSetting);
            }
        }

    }

    protected void vaccine_effect_adj_infection_duration(SiteSpecificVaccination vaccine, AbstractIndividualInterface person, int site) {
        if (vaccine.getParameters().length > SiteSpecificVaccination.OPTIONAL_EFFECT_ADJ_INF_DUR_DEFAULT) {
            double durAdj = vaccine.getParameters()[SiteSpecificVaccination.OPTIONAL_EFFECT_ADJ_INF_DUR_DEFAULT];
            if (durAdj >= 0) {
                person.setTimeUntilNextStage(site, person.getTimeUntilNextStage(site) * durAdj);
            }
        }
    }

    private void screenPerson(RelationshipPerson_MSM msm, int screenType) {
        boolean foundInfection = false;

        msm.setParameter("PARAM_LAST_SCREEN_AT_AGE", msm.getAge());

        for (int site = 0; site < msm.getInfectionStatus().length && !foundInfection; site++) {
            foundInfection |= (((1 << site) & screenType) != 0)
                    && (msm.getInfectionStatus()[site] != AbstractIndividualInterface.INFECT_S);
        }
        if (foundInfection) {

            int reIndex = 2;
            float pRepeat;
            if (((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[1] > 1) {
                reIndex += ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[1] * 2;
            }

            pRepeat = ((float[]) getFields()[MSM_SCREENING_SETTING_TARGETED])[reIndex];

            if (pRepeat > 0 && getRNG().nextFloat() < pRepeat) {
                msm.setParameter("PARAM_SCHEDULED_SCREEN_AT_AGE", msm.getAge() + 3 * 30);
            }
            treatPerson(msm);
        }

        // Vaccination by screening
        if (getFields()[MSM_SITE_SPECIFIC_VACCINATION] != null) {
            AbstractVaccination vacc = ((AbstractVaccination) getFields()[MSM_SITE_SPECIFIC_VACCINATION]);
            HashMap<Integer, int[]> currentlyVaccinated = (HashMap<Integer, int[]>) getFields()[MSM_SITE_CURRENTLY_VACCINATED];

            if (currentlyVaccinated == null || !currentlyVaccinated.containsKey(msm.getId())) {
                if (vacc instanceof SiteSpecificVaccination) {
                    SiteSpecificVaccination ssv = (SiteSpecificVaccination) vacc;
                    double propVaccThruScreening = -ssv.getParameters()[SiteSpecificVaccination.EFFECT_INDEX_PROPORTION_VACC_COVERAGE_SETTING];
                    if (propVaccThruScreening > 0) {
                        if (getRNG().nextDouble() < propVaccThruScreening) {
                            vaccinatePerson(vacc, msm);
                        }
                    }
                }
            }
        }

    }

    protected void treatPerson(RelationshipPerson_MSM msm) {
        for (int site = 0; site < msm.getInfectionStatus().length; site++) {
            if (msm.getInfectionStatus()[site] != AbstractIndividualInterface.INFECT_S) {
                float[] effByStrain = ((float[][]) getFields()[MSM_TREATMENT_EFFICIENCY])[site];

                for (int strain = 0; strain < effByStrain.length; strain++) {
                    int strainAtSite = msm.getCurrentStrainsAtSite()[site];

                    if ((strainAtSite & 1 << strain) != 0) {
                        float eff = effByStrain[strain];

                        if (eff < 1) {
                            eff = getRNG().nextFloat() < eff ? 1 : 0;
                        }
                        if (eff >= 1) {
                            msm.getCurrentStrainsAtSite()[site] -= (1 << strain);
                        }
                    }
                }

            }

            if (msm.getCurrentStrainsAtSite()[site] == MultiSiteMultiStrainPersonInterface.STRAIN_NONE) {
                msm.getInfectionStatus()[site] = AbstractIndividualInterface.INFECT_S;
                msm.setTimeUntilNextStage(site, Double.POSITIVE_INFINITY);
            }

        }
    }

    private void defineScreeningSchedule(float[] scrSetting) {

        // Set screening if needed
        // Format: float[]{startTime, totalScreened, frequency, 
        //    cummuative probabiliy of screening at site (2^3 -1) in total}    
        //  e.g. 0b001 = screen at genital only, 0b011 = screen at genital and anal
        // private transient int[] screeningPerson // Length of total screened
        // private transient int[] screeningSite // Length of total screened
        // private transient int screenPt = 0;
        // private transient int[] screeningToday // Length of screening freq
        if ((getGlobalTime() - scrSetting[0]) % scrSetting[2] == 0) {
            if (screeningPersonByIndex == null) {
                screeningPersonByIndex = new int[(int) scrSetting[1]];
                screeningSite = new int[(int) scrSetting[1]];
                screeningToday = new int[(int) scrSetting[2]];
            }
            screenDayPt = 0;
            screenTarPt = 0;
            float[] cummulaSiteProb = Arrays.copyOfRange(scrSetting, 3, scrSetting.length);
            int num2SrcTotal = screeningPersonByIndex.length;
            int num2SrcPt = 0;

            for (int pt = 0; pt < getPop().length; pt++) {
                if (getRNG().nextInt(getPop().length - pt) < num2SrcTotal) {
                    screeningPersonByIndex[num2SrcPt] = pt;

                    if (cummulaSiteProb.length > 0) {
                        int screeningType = 0;
                        float pScrType = getRNG().nextFloat();
                        while (pScrType >= cummulaSiteProb[screeningType]) {
                            screeningType++;
                        }
                        screeningSite[num2SrcPt] = screeningType + 1;
                    } else {
                        // All sites
                        screeningSite[num2SrcPt] = 0b111;
                    }

                    num2SrcPt++;
                    num2SrcTotal--;
                }
            }
            // Randomise screening schedule                
            util.ArrayUtilsRandomGenerator.shuffleArray(screeningPersonByIndex, getRNG());

            int num2SrcPerDay = screeningPersonByIndex.length / screeningToday.length;
            int numExtra = screeningPersonByIndex.length - num2SrcPerDay * screeningToday.length;

            Arrays.fill(screeningToday, num2SrcPerDay);

            int start = 0;
            for (int d = 0; d < screeningToday.length; d++) {
                if (getRNG().nextInt(screeningToday.length - d) < numExtra) {
                    screeningToday[d]++;
                    numExtra--;
                }
                // Sort the id for screening at the same day
                Arrays.sort(screeningPersonByIndex, start, start + screeningToday[d]);
                start += screeningToday[d];
            }
        }
    }

    @Override
    protected AbstractIndividualInterface generateNewPerson(int nextId,
            AbstractIndividualInterface p, double newAge) {
        RelationshipPerson_MSM person = (RelationshipPerson_MSM) p;
        RelationshipPerson_MSM newP = new RelationshipPerson_MSM(nextId,
                true, newAge, p.getInfectionStatus().length);
        newP.setEnterPopulationAt(getGlobalTime());
        newP.setParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX),
                person.getParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_BEHAV_TYPE_INDEX)));
        newP.setMaxPartners(person.getMaxPartners());

        // Act specific immunity
        newP.setParameter("PARAM_IMMUNE_ACT_SITE_G", p.getParameter("PARAM_IMMUNE_ACT_SITE_G"));
        newP.setParameter("PARAM_IMMUNE_ACT_SITE_A", p.getParameter("PARAM_IMMUNE_ACT_SITE_A"));
        newP.setParameter("PARAM_IMMUNE_ACT_SITE_R", p.getParameter("PARAM_IMMUNE_ACT_SITE_R"));

        return newP;
    }

    protected AbstractInfection[] incrementPersonStat(AbstractIndividualInterface person, int deltaT) {
        AbstractInfection[] infList = getIndivdualInfectionList(person);
        int[] preInfectStat = Arrays.copyOf(person.getInfectionStatus(), person.getInfectionStatus().length);
        person.incrementTime(deltaT, infList);
        int[] afterInfectStat = person.getInfectionStatus();

        boolean[] justBecomeInfectious = new boolean[infList.length];
        Arrays.fill(justBecomeInfectious, false);

        boolean justDevelopedSym = false;
        double symDur = Double.POSITIVE_INFINITY;

        boolean symSubside = false;

        for (int i = 0; i < preInfectStat.length; i++) {

            justBecomeInfectious[i]
                    = (preInfectStat[i] != GonorrhoeaSiteInfection.STATUS_SYM
                    && preInfectStat[i] != GonorrhoeaSiteInfection.STATUS_ASY)
                    && infList[i].isInfectious(person);

            justDevelopedSym |= preInfectStat[i] != GonorrhoeaSiteInfection.STATUS_SYM
                    && afterInfectStat[i] == GonorrhoeaSiteInfection.STATUS_SYM;

            symSubside |= preInfectStat[i] == GonorrhoeaSiteInfection.STATUS_SYM
                    && afterInfectStat[i] != GonorrhoeaSiteInfection.STATUS_SYM;

            if (justDevelopedSym) {
                symDur = Math.min(symDur, person.getTimeUntilNextStage(i));
            }

        }
        if (justDevelopedSym) {
            // Possibly reducing duration of other sites as well
            for (int i = 0; i < preInfectStat.length; i++) {
                if (afterInfectStat[i] == GonorrhoeaSiteInfection.STATUS_ASY) {
                    boolean hasReducedDur = ((float[]) getFields()[MSM_SYM_DURATION_OVERWRITE_PROB])[i] == 1f;
                    if (!hasReducedDur && ((float[]) getFields()[MSM_SYM_DURATION_OVERWRITE_PROB])[i] > 0) {
                        hasReducedDur = infList[i].getRNG().nextFloat()
                                < ((float[]) getFields()[MSM_SYM_DURATION_OVERWRITE_PROB])[i];
                    }
                    if (hasReducedDur) {
                        person.setTimeUntilNextStage(i, symDur);
                    }
                }
            }
        }
        // Effect of Vaccination
        if (getFields()[MSM_SITE_CURRENTLY_VACCINATED] != null
                && getFields()[MSM_SITE_SPECIFIC_VACCINATION] != null) {
            HashMap<Integer, int[]> currentVaccinated = (HashMap<Integer, int[]>) getFields()[MSM_SITE_CURRENTLY_VACCINATED];
            SiteSpecificVaccination vaccine = (SiteSpecificVaccination) getFields()[MSM_SITE_SPECIFIC_VACCINATION];
            int[] vacRes = currentVaccinated.get(person.getId());

            if (vacRes != null) {
                if (person.getAge() < vacRes[VACC_SETTING_AGE_EXPIRY]) { // Active vaccine
                    // Just become infectious 
                    for (int site = 0; site < getInfList().length; site++) {
                        if (justBecomeInfectious[site]) {
                            vaccine_effect_adj_infection_duration(vaccine, person, site);
                        }
                        if (justDevelopedSym && getInfList()[site].hasSymptoms(person)) {
                            vaccine_effect_remove_symptom(vaccine, person, site);
                        }

                    }
                }
            }
        }

        if (symSubside) {
            treatPerson((RelationshipPerson_MSM) person);
        }

        return infList;

    }

    protected boolean vaccine_effect_remove_symptom(SiteSpecificVaccination vaccine, AbstractIndividualInterface person, int site) {
        boolean removeSym = false;
        if (vaccine.getParameters().length > SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_SD) {
            double sym_adj = vaccine.getParameters()[SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_RATE_DEFAULT];
            if (sym_adj > 0) {
                removeSym = sym_adj >= 1;
                if (sym_adj < 1) {
                    removeSym = getInfList()[site].getRNG().nextDouble() < sym_adj;
                }
                if (removeSym) {
                    person.getInfectionStatus()[site] = (int) vaccine.getParameters()[OPTIONAL_EFFECT_REMOVE_SYM_STATE_DEFAULT];
                    double newDuration = vaccine.getParameters()[SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_MEDIAN];
                    if (vaccine.getParameters()[SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_SD] > 0) {
                        if (vaccine_remove_sym_infect_duration == null) {
                            double[] ent = StaticMethods.generatedGammaParam(new double[]{newDuration,
                                vaccine.getParameters()[SiteSpecificVaccination.OPTIONAL_EFFECT_REMOVE_SYM_INF_DUR_DEFAULT_SD]});
                            vaccine_remove_sym_infect_duration = new GammaDistribution(ent[0], 1 / ent[1]);
                        }
                        newDuration = (int) Math.round(vaccine_remove_sym_infect_duration.sample());

                        person.setTimeUntilNextStage(site, newDuration);

                    }
                }
            }
        }
        return removeSym;
    }

    @Override
    protected boolean[][] performAct(RegCasRelationship rel) {

        boolean[] hasActed = rel.hasActToday();
        boolean[][] res = new boolean[hasActed.length][3]; // hasUnprotectedSex[actType]{occured, from_person_1, from_person_2}

        int[][] infectStat = new int[rel.getLinks().length][];
        int[][] strainStat = new int[rel.getLinks().length][];

        double[][] vaccineImpact = new double[rel.getLinks().length][];

        RelationshipPerson_MSM[] person = new RelationshipPerson_MSM[rel.getLinks().length];

        for (int p = 0; p < rel.getLinks().length; p++) {
            person[p] = (RelationshipPerson_MSM) getLocalData().get(rel.getLinks()[p].intValue());
            infectStat[p] = person[p].getInfectionStatus();
            strainStat[p] = person[p].getCurrentStrainsAtSite();

            if (getFields()[MSM_SITE_SPECIFIC_VACCINATION] != null) {
                SiteSpecificVaccination vacc = (SiteSpecificVaccination) getFields()[MSM_SITE_SPECIFIC_VACCINATION];

                int[] vacState = null;

                if (getFields()[MSM_SITE_CURRENTLY_VACCINATED] != null) {
                    vacState = ((HashMap<Integer, int[]>) getFields()[MSM_SITE_CURRENTLY_VACCINATED]).get(person[p].getId());
                }

                if (vacState != null
                        && (vacState[VACC_SETTING_AGE_EXPIRY] < 0
                        || vacState[VACC_SETTING_AGE_EXPIRY] > person[p].getAge())) {
                    vaccineImpact[p] = vacc.vaccineImpact(person[p], null);
                } else {
                    vaccineImpact[p] = null;
                }
            }

        }

        float[] protectAdj = rel.getType() == RegCasRelationship.REL_TYPE_REG
                ? (float[]) getFields()[MSM_REG_CONDOM_USAGE_ADJ]
                : (float[]) getFields()[MSM_CAS_CONDOM_USAGE_ADJ];

        for (int a = 0; a < res.length; a++) {
            if (hasActed[a]) {
                switch (a) {
                    case ACT_KISSING:
                        res[a][0] = true;
                        res[a][1] = false;
                        res[a][2] = false;

                        double r2r_def = ((double[][]) getFields()[FIELDS_TRANSMIT])[TRAN_SUSC_INDEX_KISSING][0]
                                * ((double[][]) getFields()[FIELDS_SUSCEPT])[TRAN_SUSC_INDEX_KISSING][0];

                        boolean tranR2R = r2r_def > 0
                                && strainStat[0][RelationshipPerson_MSM.SITE_R] != strainStat[1][RelationshipPerson_MSM.SITE_R];

                        if (tranR2R) {
                            for (int s = 0; s < rel.getLinks().length; s++) {
                                if (strainStat[s][RelationshipPerson_MSM.SITE_R] != strainStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_R]
                                        && (infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_R] == GonorrhoeaSiteInfection.STATUS_ASY
                                        || infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_R] == GonorrhoeaSiteInfection.STATUS_SYM)) {

                                    double r2r = r2r_def;

                                    // Impact of vaccine
                                    if (vaccineImpact[s] != null) {
                                        r2r *= vaccineImpact[s][SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_R];
                                    }
                                    if (vaccineImpact[(s + 1) % 2] != null) {
                                        r2r *= vaccineImpact[(s + 1) % 2][SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_R];
                                    }

                                    if (r2r < 1) {
                                        tranR2R = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_R]).getRNG().nextDouble() < r2r;
                                    }

                                    if ((((int) person[s].getParameter("PARAM_IMMUNE_ACT_SITE_R")) & (1 << ACT_KISSING)) != 0) {
                                        tranR2R = false;
                                    }

                                    if (tranR2R) {
                                        res[a][(s + 1) % 2] = true;
                                        cumulativeIncidencesBySites[RelationshipPerson_MSM.SITE_R]++;
                                        person[s].setLastActInfectious(RelationshipPerson_MSM.SITE_R, true);
                                        if (person[s] instanceof MultiSiteMultiStrainPersonInterface) {
                                            ((MultiSiteMultiStrainPersonInterface) person[s]).getLastActStainsAtSite()[RelationshipPerson_MSM.SITE_R]
                                                    = ((MultiSiteMultiStrainPersonInterface) person[(s + 1) % 2]).getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_R];
                                        }
                                    }
                                }
                            }
                        }
                        break;
                    case ACT_RIMMING:
                        res[a][0] = true;
                        res[a][1] = false;
                        res[a][2] = false;

                        double a2r_def = ((double[][]) getFields()[FIELDS_TRANSMIT])[TRAN_SUSC_INDEX_RIMMING_ANAL][0]
                                * ((double[][]) getFields()[FIELDS_SUSCEPT])[TRAN_SUSC_INDEX_RIMMING_ORAL][0];
                        double r2a_def = ((double[][]) getFields()[FIELDS_SUSCEPT])[TRAN_SUSC_INDEX_RIMMING_ANAL][0]
                                * ((double[][]) getFields()[FIELDS_TRANSMIT])[TRAN_SUSC_INDEX_RIMMING_ORAL][0];

                        boolean tranA2R = a2r_def > 0;
                        boolean tranR2A = r2a_def > 0;

                        for (int s = 0; s < rel.getLinks().length; s++) {
                            boolean susR = strainStat[s][RelationshipPerson_MSM.SITE_R] != strainStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_A];
                            boolean susA = strainStat[s][RelationshipPerson_MSM.SITE_A] != strainStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_R];

                            double a2r = a2r_def;
                            double r2a = r2a_def;

                            // Impact of vaccine
                            if (vaccineImpact[s] != null) {
                                a2r *= vaccineImpact[s][SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_R];
                                r2a *= vaccineImpact[s][SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_A];
                            }
                            if (vaccineImpact[(s + 1) % 2] != null) {
                                a2r *= vaccineImpact[(s + 1) % 2][SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_A];
                                r2a *= vaccineImpact[(s + 1) % 2][SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_R];
                            }

                            // Anal to oral
                            if (tranA2R && susR
                                    && (infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_A] == GonorrhoeaSiteInfection.STATUS_ASY
                                    || infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_A] == GonorrhoeaSiteInfection.STATUS_SYM)) {

                                if (a2r < 1) {
                                    tranA2R = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_A]).getRNG().nextDouble() < a2r;
                                }
                                if ((((int) person[s].getParameter("PARAM_IMMUNE_ACT_SITE_R")) & (1 << ACT_RIMMING)) != 0) {
                                    tranA2R = false;
                                }
                                if (tranA2R && !person[s].getLastActInfectious()[RelationshipPerson_MSM.SITE_R]) {
                                    cumulativeIncidencesBySites[RelationshipPerson_MSM.SITE_R]++;
                                    res[a][(s + 1) % 2] = true;
                                    person[s].setLastActInfectious(RelationshipPerson_MSM.SITE_R, true);

                                    if (person[s] instanceof MultiSiteMultiStrainPersonInterface) {
                                        ((MultiSiteMultiStrainPersonInterface) person[s]).getLastActStainsAtSite()[RelationshipPerson_MSM.SITE_R] = ((MultiSiteMultiStrainPersonInterface) person[(s + 1) % 2]).getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_A];
                                    }
                                }
                            }

                            // Oral to anal
                            if (tranR2A && susA
                                    && (infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_R] == GonorrhoeaSiteInfection.STATUS_ASY
                                    || infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_R] == GonorrhoeaSiteInfection.STATUS_SYM)) {

                                if (r2a < 1) {
                                    tranR2A = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_R]).getRNG().nextDouble() < r2a;
                                }
                                if ((((int) person[s].getParameter("PARAM_IMMUNE_ACT_SITE_A")) & (1 << ACT_RIMMING)) != 0) {
                                    tranR2A = false;
                                }
                                if (tranR2A && !person[s].getLastActInfectious()[RelationshipPerson_MSM.SITE_A]) {
                                    res[a][(s + 1) % 2] = true;
                                    cumulativeIncidencesBySites[RelationshipPerson_MSM.SITE_A]++;
                                    person[s].setLastActInfectious(RelationshipPerson_MSM.SITE_A, true);
                                    if (person[s] instanceof MultiSiteMultiStrainPersonInterface) {
                                        ((MultiSiteMultiStrainPersonInterface) person[s]).getLastActStainsAtSite()[RelationshipPerson_MSM.SITE_A] = ((MultiSiteMultiStrainPersonInterface) person[(s + 1) % 2]).getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_R];
                                    }
                                }

                            }

                        }
                        break;
                    case ACT_ANAL:
                    case ACT_ORAL:
                        int nonGTarget = a == ACT_ANAL ? RelationshipPerson_MSM.SITE_A : RelationshipPerson_MSM.SITE_R;
                        String nonGTargetImmune = a == ACT_ANAL ? "PARAM_IMMUNE_ACT_SITE_A" : "PARAM_IMMUNE_ACT_SITE_R";
                        float probCondomUse = rel.getType() == RegCasRelationship.REL_TYPE_REG
                                ? ((float[]) getFields()[MSM_POP_REG_CONDOM_USAGE])[a]
                                : ((float[]) getFields()[MSM_POP_CAS_CONDOM_USAGE])[a];

                        // Condom usage adjust
                        if (protectAdj != null && protectAdj.length > 0) {
                            if (getGlobalTime() >= protectAdj[0]) {
                                if (protectAdj[a + 1] >= 0) {
                                    probCondomUse *= protectAdj[a + 1];
                                } else {
                                    probCondomUse = -protectAdj[a + 1]; // If < 0, replacement instead
                                }
                            }
                        }

                        // Determine if condom is use for the act
                        boolean unprotectedAct = probCondomUse == 0;
                        if (!unprotectedAct && probCondomUse < 1) {
                            unprotectedAct = getRNG().nextFloat() > probCondomUse;
                        }

                        for (int s = 0; s < rel.getLinks().length; s++) {

                            if (unprotectedAct) {
                                if (a == ACT_ANAL) {
                                    person[s].setParameter("PARAM_LAST_UNPROTECTED_ANAL_SEX_AT_AGE", person[s].getAge());
                                } else if (a == ACT_ORAL) {
                                    person[s].setParameter("PARAM_LAST_UNPROTECTED_ORAL_SEX_AT_AGE", person[s].getAge());
                                }
                                res[a][0] = true;

                                double tranSucProb;
                                boolean susNonG = strainStat[s][nonGTarget] != strainStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G];
                                boolean susG = strainStat[s][RelationshipPerson_MSM.SITE_G] != strainStat[(s + 1) % 2][nonGTarget];

                                if (susNonG && (infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G] == GonorrhoeaSiteInfection.STATUS_ASY
                                        || infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G] == GonorrhoeaSiteInfection.STATUS_SYM)) {

                                    // From G to nonG - set to global version if available
                                    if (getFields()[FIELDS_SUSCEPT] != null) {
                                        if (person[s].getProbSusBySite()[nonGTarget] != ((double[][]) getFields()[FIELDS_SUSCEPT])[nonGTarget][0]) {
                                            person[s].getProbSusBySite()[nonGTarget] = ((double[][]) getFields()[FIELDS_SUSCEPT])[nonGTarget][0];
                                        }
                                    }

                                    tranSucProb = person[(s + 1) % 2].getProbTransBySite()[RelationshipPerson_MSM.SITE_G]
                                            * person[s].getProbSusBySite()[nonGTarget];

                                    if ((((int) person[s].getParameter(nonGTargetImmune)) & (1 << a)) != 0) {
                                        // E.g. Insertive to anal or oral sex
                                        tranSucProb = 0;
                                    }

                                    // Impact of vaccine
                                    if (vaccineImpact[s] != null) {
                                        int eff_tar = nonGTarget == RelationshipPerson_MSM.SITE_A
                                                ? SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_A
                                                : SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_R;
                                        tranSucProb *= vaccineImpact[s][eff_tar];
                                    }
                                    if (vaccineImpact[(s + 1) % 2] != null) {
                                        tranSucProb *= vaccineImpact[(s + 1) % 2][SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_G];
                                    }

                                    boolean g2NG = tranSucProb > 0;
                                    if (g2NG && tranSucProb < 1) {
                                        g2NG = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_G]).getRNG().nextDouble()
                                                < tranSucProb;
                                    }

                                    if (g2NG) {
                                        res[a][(s + 1) % 2] = true;
                                        cumulativeIncidencesBySites[nonGTarget]++;
                                        person[s].setLastActInfectious(nonGTarget, true);
                                        if (person[s] instanceof MultiSiteMultiStrainPersonInterface) {
                                            ((MultiSiteMultiStrainPersonInterface) person[s]).getLastActStainsAtSite()[nonGTarget] = ((MultiSiteMultiStrainPersonInterface) person[(s + 1) % 2]).getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_G];
                                        }

                                    }

                                }
                                if (susG && (infectStat[(s + 1) % 2][nonGTarget] == GonorrhoeaSiteInfection.STATUS_ASY
                                        || infectStat[(s + 1) % 2][nonGTarget] == GonorrhoeaSiteInfection.STATUS_SYM)) {

                                    // From nonG to G - set to global version if available
                                    if (getFields()[FIELDS_SUSCEPT] != null) {
                                        if (person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G]
                                                != ((double[][]) getFields()[FIELDS_SUSCEPT])[RelationshipPerson_MSM.SITE_G][0]) {
                                            person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G] = ((double[][]) getFields()[FIELDS_SUSCEPT])[RelationshipPerson_MSM.SITE_G][0];
                                        }
                                    }

                                    tranSucProb = person[(s + 1) % 2].getProbTransBySite()[nonGTarget]
                                            * person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G];

                                    if ((((int) person[s].getParameter("PARAM_IMMUNE_ACT_SITE_G")) & (1 << a)) != 0) {
                                        // Eg. Receptive to anal or oral sex
                                        tranSucProb = 0;
                                    }

                                    // Impact of vaccine
                                    if (vaccineImpact[s] != null) {
                                        tranSucProb *= vaccineImpact[s][SiteSpecificVaccination.EFFECT_INDEX_SUSCEPTIBLE_EFFICACY_G];
                                    }
                                    if (vaccineImpact[(s + 1) % 2] != null) {
                                        int eff_tar = nonGTarget == RelationshipPerson_MSM.SITE_A
                                                ? SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_A
                                                : SiteSpecificVaccination.EFFECT_INDEX_TRANMISSION_EFFICACY_R;
                                        tranSucProb *= vaccineImpact[(s + 1) % 2][eff_tar];
                                    }

                                    boolean nG2G = tranSucProb > 0;

                                    if (nG2G && tranSucProb < 1) {
                                        nG2G = ((GonorrhoeaSiteInfection) getInfList()[nonGTarget]).getRNG().nextDouble()
                                                < tranSucProb;
                                    }
                                    if (nG2G) {
                                        res[a][(s + 1) % 2] = true;
                                        cumulativeIncidencesBySites[RelationshipPerson_MSM.SITE_G]++;
                                        person[s].setLastActInfectious(RelationshipPerson_MSM.SITE_G, true);
                                        if (person[s] instanceof MultiSiteMultiStrainPersonInterface) {
                                            ((MultiSiteMultiStrainPersonInterface) person[s]).getLastActStainsAtSite()[RelationshipPerson_MSM.SITE_G] = ((MultiSiteMultiStrainPersonInterface) person[(s + 1) % 2]).getCurrentStrainsAtSite()[nonGTarget];
                                        }
                                    }

                                }
                            }
                        }
                        break;
                }
            }
        }

        return res;
    }
    public int[] relTotal = new int[2];
    public int[] relLen = new int[2];

    @Override
    protected SingleRelationship formRelationship(AbstractIndividualInterface[] pair, RelationshipMap relMap, int d, int mapType) {
        RegCasRelationship rel;
        int dur;
        float[][] actFreq;
        if (mapType == MAPPING_REG) {

            double regPartLength = ((Number) getFields()[MSM_REG_LENTH_AVE]).doubleValue();
            AbstractRealDistribution regLength;
            regLength = new ExponentialDistribution(getRNG(), regPartLength);
            dur = (int) Math.max(Math.round(regLength.sample()), 1);
            actFreq = (float[][]) getFields()[MSM_POP_REG_ACT_FREQ];

        } else {
            dur = Integer.MAX_VALUE;
            for (AbstractIndividualInterface pair1 : pair) {
                RelationshipPerson_MSM person = (RelationshipPerson_MSM) pair1;
                int numCasualPast6Months = ((Number) person.getParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_NUM_CASUAL_IN_LAST_6_MONTHS))).intValue();
                int maxCasual = person.getMaxPartners();
                if (maxCasual > 0 && maxCasual > numCasualPast6Months) {
                    dur = Math.min(dur, 6 * 30 / (maxCasual - numCasualPast6Months));
                } else {
                    dur = 1; // Dummy relationship of 1 day
                }
            }
            // Min of 1 day
            dur = Math.max(1, dur);
            actFreq = (float[][]) getFields()[MSM_POP_CAS_ACT_FREQ];
        }

        actFreq = Arrays.copyOf(actFreq, actFreq.length); // So it won't change the global value 

        if (actFreq.length > 0 && actFreq[0][0] == -1) {
            // For Alterative format: {[-1, number partners limit], 
            //                      [prob of anal under limit, prob of anal over limit],
            //                      [prob of oral under limit, prob of oral over limit]... }                   

            float[][] actFreqAlt = new float[actFreq.length - 1][2];
            for (AbstractIndividualInterface pair1 : pair) {
                RelationshipPerson_MSM person = (RelationshipPerson_MSM) pair1;
                int numCasualPast6Months = ((Number) person.getParameter(person.indexToParamName(RelationshipPerson_MSM.PARAM_NUM_CASUAL_IN_LAST_6_MONTHS))).intValue();
                for (int a = 0; a < actFreqAlt.length; a++) {
                    actFreqAlt[a][0] += actFreq[a + 1][(numCasualPast6Months <= actFreq[0][1]) ? 0 : 1] / 2;
                    actFreqAlt[a][1] = actFreqAlt[a][0];
                }
            }

            actFreq = actFreqAlt;
        }

        relTotal[mapType]++;
        relLen[mapType] += dur;

        for (int a = 0; a < actFreq.length; a++) { // Check for 3 parameters option
            float[] freqForAct = actFreq[a];
            if (freqForAct.length > 2) {  // Alterative implementation , [min,max, prob of occurring (cumulative)]                
                actFreq[a] = new float[2];
                float pActInRel = getRNG().nextFloat();
                int pt = 0;
                while (pt + 2 < freqForAct.length
                        && pActInRel >= freqForAct[pt + 2]) {
                    pt += 2;
                }
                if (pt + 2 < freqForAct.length) {
                    actFreq[a] = Arrays.copyOfRange(freqForAct, pt, pt + 2);
                }
            }
        }

        rel = new RegCasRelationship(new Integer[]{pair[0].getId(), pair[1].getId()}, mapType, actFreq.length);
        rel.setDurations(Math.max((int) dur, 1)); // Relations of at least one day

        for (int a = 0; a < actFreq.length; a++) {
            float actFloat = actFreq[a][0];
            if (actFloat != actFreq[a][1]) {
                if (actFreq[a][1] >= 1) { // As whole number                    
                    actFloat += getRNG().nextInt((int) (actFreq[a][1] - actFreq[a][0]));
                } else if (actFreq[a][1] > 0) {
                    actFloat += getRNG().nextFloat() * (actFreq[a][1] - actFreq[a][0]);
                }

                float[] actAdj = (float[]) getFields()[mapType == MAPPING_REG ? MSM_REG_FREQ_ADJ : MSM_CAS_FREQ_ADJ];

                if (actAdj != null && actAdj.length > 0) {
                    if (getGlobalTime() >= actAdj[0] && a + 1 < actAdj.length) {
                        boolean isRatio = actFloat < 1;
                        actFloat *= actAdj[a + 1];

                        if (isRatio) {
                            actFloat = Math.min(actFloat, 0.9999f); // Cannot go beyond 1
                        }

                    }
                }
            }
            rel.setActSchedule(a, actFloat, getRNG());
        }
        for (AbstractIndividualInterface pair1 : pair) {
            if (!relMap.containsVertex(pair1.getId())) {
                relMap.addVertex(pair1.getId());
            }
        }
        boolean added = relMap.addEdge(pair[0].getId(), pair[1].getId(), rel);

        if (added) {

            for (AbstractIndividualInterface person : pair) {
                relMap.removeAvailablePerson(person);
            }

            if (mapType == MAPPING_CAS) {
                for (int p = 0; p < pair.length; p++) {
                    RelationshipPerson_MSM person = (RelationshipPerson_MSM) pair[p];
                    person.addCasualPartner(pair[(p + 1) % 2]);
                }
            }
        } else {
            rel = null;
        }

        return rel;
    }

    @Override
    public void initialiseInfection(long seed) {
        AbstractInfection[] infList = new AbstractInfection[3]; // 3 sites        
        random.RandomGenerator infRNG;
        float[] probSymBySite = (float[]) getFields()[MSM_SYM_TREATMENT_PROB];

        if (seed != 0) {
            getFields()[FIELDS_INF_RNG] = new random.RandomGenerator[infList.length];
            infRNG = new MersenneTwisterRandomGenerator(seed);
        } else {
            infRNG = null;
        }
        for (int i = 0; i < infList.length; i++) {
            double[][] siteParam = new double[GonorrhoeaSiteInfection.DIST_TOTAL][];
            AbstractRealDistribution[] distributions = new AbstractRealDistribution[GonorrhoeaSiteInfection.DIST_TOTAL];

            if (infRNG == null) {
                infRNG = ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[i]; // Imported
            }

            ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[i] = infRNG;
            siteParam[GonorrhoeaSiteInfection.DIST_EXPOSED_DUR_INDEX] = new double[]{4, 0};
            siteParam[GonorrhoeaSiteInfection.DIST_IMMUNE_DUR_INDEX] = new double[]{7, 0};
            siteParam[GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX] = new double[]{1, 0};
            double[] var;
            switch (i) {
                case RelationshipPerson_MSM.SITE_G:
                    // Assume all have sym i.e treated
                    siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{probSymBySite[i], 0};
                    // Kit email 20130911 for urethal
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new double[]{2.5714, 2.24343};
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new GammaDistribution(infRNG, var[0], 1 / var[1]);
                    // Assume duration same as Garnett 1999, Brunham 1991, SD same as Johnson 2010                                        
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{185, 5 * 7};
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new GammaDistribution(infRNG, var[0], 1 / var[1]);
                    // Transmission from G - 0.5 Chen 2010
                    siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = new double[]{0.5, 0};
                    break;
                case RelationshipPerson_MSM.SITE_A:
                    // Anal : Bissessor 2011: 7 out 47 of has proctitis
                    siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{probSymBySite[i], 0};
                    // Assume same as Urethal
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new double[]{2.5714, 2.24343};
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new GammaDistribution(infRNG, var[0], 1 / var[1]);
                    // Assume duration same as Garnett 1999, Brunham 1991, SD same as Johnson 2010
                    // siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{185, 5 * 7};
                    // From C. Fairley email at 20131120
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{12 * 30, 5 * 7};
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new GammaDistribution(infRNG, var[0], 1 / var[1]);
                    // Transmission from A - assumption
                    siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = new double[]{0.08, 0};
                    break;
                case RelationshipPerson_MSM.SITE_R:
                    // Assume none have sym
                    siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{probSymBySite[i], 0};
                    // Assume same as Johnson 2010
                    //siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{15 * 7, 5 * 7};
                    // From C. Fairley email at 20131121, Fairley et al 2011
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{12 * 7, 3 * 7};
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new GammaDistribution(infRNG, var[0], 1 / var[1]);
                    // Transmission from R - assumption
                    siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = new double[]{0.4, 0};
                    break;
            }
            if (getFields()[FIELDS_TRANSMIT] != null) {
                siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = ((double[][]) getFields()[FIELDS_TRANSMIT])[i];
            }
            if (getFields()[FIELDS_SUSCEPT] != null) {
                siteParam[GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX] = ((double[][]) getFields()[FIELDS_SUSCEPT])[i];
            }
            infList[i] = new GonorrhoeaSiteInfection(infRNG, i, siteParam, distributions);
        }
        setInfList(infList);
    }

    protected void snapshotCount(AbstractIndividualInterface person) {
        for (int i = 0; i < getIndivdualInfectionList(person).length; i++) {
            if (person.getInfectionStatus(i) != AbstractIndividualInterface.INFECT_S) {
                getNumInf()[i]++;
            }
        }
        if (getSnapshotClassifier() != null) {
            for (int i = 0; i < snapshotClassifier.length; i++) {
                if (snapshotClassifier[i] != null && snapshotCount[i] != null) {
                    int cI = snapshotClassifier[i].classifyPerson(person);
                    if (cI >= 0) {
                        snapshotCount[i][cI]++;
                    }
                }
            }
        }
    }

}
