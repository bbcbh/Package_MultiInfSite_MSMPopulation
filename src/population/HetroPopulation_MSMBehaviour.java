package population;


import infection.AbstractInfection;
import infection.GonorrhoeaSiteInfection;
import java.util.Arrays;
import org.apache.commons.math3.distribution.GammaDistribution;
import availability.AbstractAvailability;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import population.availability.MSMAvailbility_Bipartite;
import person.AbstractIndividualInterface;
import population.person.RelationshipPerson_MSM;
import population.person.SiteSpecificTransmissivity;
import population.relationshipMap.RegCasRelationship;
import random.RandomGenerator;

import util.StaticMethods;

/**
 * A heterosexual population with MSM behaviour.
 *
 * The difference are:
 * <ol>
 * <li>Bipartiate - for now assume female has ID of odd number</li>
 * <li>No little anal sex</li>
 * <li>Genital - Genital transmission</li>
 * </ol>
 *
 * <p>
 * History:</p>
 * <p>
 * 20140522 - Add support for HETRO_POP_SYM_DURATION_ADJ </p>
 * <p>
 * 20140328 - Add support for MSM_SYM_DURATION_OVERWRITE_PROB </p>
 * <p>
 * 20140313 - New per act tranmission prob. (see for example, <a herf=http://en.wikipedia.org/wiki/Sexually_transmitted_disease">here</a>) </p>
 * <p>
 * 20140228 - Move availability into proper subclass</p>
 * <p>
 * 20140128 - Add proportional anal sex to genital transmission</p>
 * <p>
 * 20140129 - Use product for transmission probability. Shared RNG for all infection </p>
 * <p>
 * 20140130 - Added MSM_CAS_PART_PROB and MSM_CAS_PART_SPREAD field </p>
 *
 * @author Ben Hui
 * @version 20140522
 *
 */
public class HetroPopulation_MSMBehaviour extends MSMPopulation {

    public static final int HETRO_POP_PROB_ANAL_INDEX = LENGTH_FIELDS_MSM_POP;
    public static final int HETRO_POP_TREAT_RATE_GENITAL = HETRO_POP_PROB_ANAL_INDEX + 1;
    public static final int HETRO_POP_SYM_DURATION_ADJ = HETRO_POP_TREAT_RATE_GENITAL + 1;
    public static final int LENGTH_FIELDS_HETRO_POP = HETRO_POP_SYM_DURATION_ADJ + 1;

    public static final Object[] DEFAULT_HETRO_FIELDS = {
        new Float(0.008), // Proportion of genital acts to anal - Maybe de Visser 2003. If negative, it is proportion of anal + genitial
        new double[][]{{1, 0}, {0.4, 0}},
        new float[][]{{1, 1, 1}, {1, 1, 1}},};

    private final double[][] DEFAULT_TRANSMIT = new double[][]{ // Genital tranmission based from Holmes 1970, Platt 1983 and others,
        //{0.50, 0.0}, {0.0, 0.0},                 { 0.08498397168684897, 0.0},
        //{0.22, 0.0}, {0.026149873567442243, 0.0}, { 0.08498397168684897, 0.0},};    
        {0.50, 0.0}, {0.0, 0.0}, {0.08654020531624781, 0.0},
        {0.22, 0.0}, {0.024258237711416882, 0.0}, {0.08654020531624781, 0.0},};

    private final double[][] DEFAULT_SUSCEPT = new double[][]{
        //{1, 0.0}, {0.0, 0.0},                       { 0.6268051255667103 / 0.22, 0.0},
        //{1, 0.0}, {0.8413473872398357 / 0.50, 0.0}, { 0.6268051255667103 / 0.50, 0.0},};
        {1, 0.0}, {0.0, 0.0}, {0.6288277879484502 / 0.22, 0.0},
        {1, 0.0}, {0.8402332453833801 / 0.50, 0.0}, {0.6288277879484502 / 0.50, 0.0},};

    private final RandomGenerator actTypeRNG;

    public HetroPopulation_MSMBehaviour(long seed) {
        super(seed);

        Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_HETRO_POP);
        for (int i = LENGTH_FIELDS_MSM_POP; i < newFields.length; i++) {
            newFields[i] = DEFAULT_HETRO_FIELDS[i - LENGTH_FIELDS_MSM_POP];
        }
        super.setFields(newFields);

        getFields()[FIELDS_TRANSMIT] = DEFAULT_TRANSMIT;
        getFields()[FIELDS_SUSCEPT] = DEFAULT_SUSCEPT;

        // Rissel 2000 
        // 77.5% has single partner only, and assume 4% can have both at the same time (Rissel 2003, ASHR)
        getFields()[MSM_NUM_IN_GRP] = new float[]{0.775f / (1 - 0.123f), 1 - (1.04f * 0.775f / (1 - 0.123f)), 0.775f / (1 - 0.123f) * 0.04f};

        // Tothill ASHR data - for those more than 1 partner -
        //getFields()[MSM_CAS_PART_PROB] = new float[]{1f};
        //getFields()[MSM_CAS_PART_SPREAD] = new int[][]{{1,2}};                           
        actTypeRNG = new random.MersenneTwisterFastEngine(seed);
    }
    public static final int ACT_GENITAL = ACT_ORAL + 1;
    public static final int ACT_GENITAL_ANAL = ACT_GENITAL + 1;

    @Override
    protected boolean[][] performAct(RegCasRelationship rel) {
        boolean[] hasActed = rel.hasActToday();
        boolean[][] res = new boolean[hasActed.length][2];

        int[][] infectStat = new int[rel.getLinks().length][];
        RelationshipPerson_MSM[] person = new RelationshipPerson_MSM[rel.getLinks().length];

        for (int p = 0; p < rel.getLinks().length; p++) {
            person[p] = (RelationshipPerson_MSM) getLocalData().get(rel.getLinks()[p].intValue());
            infectStat[p] = person[p].getInfectionStatus();

        }

        for (int a = 0; a < res.length; a++) {
            if (hasActed[a]) {
                int actType = a;
                if (actType == ACT_ANAL) {
                    float pAnal = ((Float) getFields()[HETRO_POP_PROB_ANAL_INDEX]).floatValue();
                    if (pAnal == 0) {
                        actType = ACT_GENITAL;
                    } else if (pAnal > 0) {
                        actType = actTypeRNG.nextFloat() < pAnal ? ACT_ANAL : ACT_GENITAL;
                    } else if (pAnal < 0) {
                        if (pAnal <= -1) {
                            actType = ACT_GENITAL_ANAL;
                        } else {
                            actType = actTypeRNG.nextFloat() < -pAnal ? ACT_GENITAL_ANAL : ACT_GENITAL;
                        }
                    }
                }

                int siteNonG; // Although in this case it can G as well
                switch (actType) {
                    case ACT_GENITAL:
                        siteNonG = RelationshipPerson_MSM.SITE_G;
                        break;
                    case ACT_GENITAL_ANAL:
                        siteNonG = RelationshipPerson_MSM.SITE_G;
                        break;
                    case ACT_ANAL:
                        siteNonG = RelationshipPerson_MSM.SITE_A;
                        break;
                    default:
                        siteNonG = RelationshipPerson_MSM.SITE_R;
                }

                float probCondomUseProb = rel.getType() == RegCasRelationship.REL_TYPE_REG
                        ? ((float[]) getFields()[MSM_POP_REG_CONDOM_USAGE])[a]
                        : ((float[]) getFields()[MSM_POP_CAS_CONDOM_USAGE])[a];

                boolean unprotectedAct = probCondomUseProb == 0;

                if (!unprotectedAct && probCondomUseProb < 1) {
                    unprotectedAct = getRNG().nextFloat() > probCondomUseProb;
                }

                for (int s = 0; s < rel.getLinks().length; s++) {

                    if (unprotectedAct) {
                        double tranSucProb;
                        if (actType == ACT_GENITAL || actType == ACT_GENITAL_ANAL) { // Special case to prevent double count
                            if (infectStat[s][RelationshipPerson_MSM.SITE_G] == AbstractIndividualInterface.INFECT_S
                                    && (infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G] == GonorrhoeaSiteInfection.STATUS_ASY
                                    || infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G] == GonorrhoeaSiteInfection.STATUS_SYM)) {

                                // From G to G  
                                if (person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G]
                                        != ((double[][]) getFields()[FIELDS_SUSCEPT])[RelationshipPerson_MSM.SITE_G
                                        + (s == 0 ? 0 : getInfList().length / 2)][0]) {
                                    person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G]
                                            = ((double[][]) getFields()[FIELDS_SUSCEPT])[RelationshipPerson_MSM.SITE_G
                                            + (s == 0 ? 0 : getInfList().length / 2)][0];
                                }

                                tranSucProb = person[(s + 1) % 2].getProbTransBySite()[RelationshipPerson_MSM.SITE_G]
                                        * person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G];

                                boolean g2g = tranSucProb >= 1;

                                if (!g2g && tranSucProb > 0) {
                                    g2g = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_G]).getRNG().nextDouble()
                                            < tranSucProb;
                                }

                                if (g2g) {
                                    person[s].setLastActInfectious(RelationshipPerson_MSM.SITE_G, true);
                                    if (person[s] instanceof SiteSpecificTransmissivity) {
                                        ((SiteSpecificTransmissivity) person[s]).getLastActStainsAtSite()[RelationshipPerson_MSM.SITE_G]
                                                = ((SiteSpecificTransmissivity) person[(s + 1) % 2]).getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_G];
                                    }
                                }
                            }
                            if (actType == ACT_GENITAL_ANAL) { // Addtional for anal
                                nonG2GTranmission(infectStat, s, RelationshipPerson_MSM.SITE_A, ACT_ANAL, person);
                            }
                        } else {
                            nonG2GTranmission(infectStat, s, siteNonG, actType, person);
                        }
                    }
                }
            }
        }

        return res;
    }

    private void nonG2GTranmission(int[][] infectStat, int s, int siteNonG, int actType, RelationshipPerson_MSM[] person) {
        double tranSucProb;
        if (infectStat[s][siteNonG] == AbstractIndividualInterface.INFECT_S
                && (infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G] == GonorrhoeaSiteInfection.STATUS_ASY
                || infectStat[(s + 1) % 2][RelationshipPerson_MSM.SITE_G] == GonorrhoeaSiteInfection.STATUS_SYM)
                && (actType != ACT_ANAL || s != 0)) { // Male's anal cannot be infected by female's genital

            // From G to altSite
            if (person[s].getProbSusBySite()[siteNonG]
                    != ((double[][]) getFields()[FIELDS_SUSCEPT])[siteNonG
                    + (s == 0 ? 0 : getInfList().length / 2)][0]) {
                person[s].getProbSusBySite()[siteNonG]
                        = ((double[][]) getFields()[FIELDS_SUSCEPT])[siteNonG + (s == 0 ? 0 : getInfList().length / 2)][0];
            }

            tranSucProb = person[(s + 1) % 2].getProbTransBySite()[RelationshipPerson_MSM.SITE_G]
                    * person[s].getProbSusBySite()[siteNonG];

            boolean g2AltSite = tranSucProb >= 1;

            if (!g2AltSite && tranSucProb > 0) {
                g2AltSite = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_G]).getRNG().nextDouble()
                        < tranSucProb;
            }

            if (g2AltSite) {
                person[s].setLastActInfectious(siteNonG, true);

                if (person[s] instanceof SiteSpecificTransmissivity) {
                    ((SiteSpecificTransmissivity) person[s]).getLastActStainsAtSite()[siteNonG]
                            = ((SiteSpecificTransmissivity) person[(s + 1) % 2]).getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_G];
                }

            }

        }
        if (infectStat[s][RelationshipPerson_MSM.SITE_G] == AbstractIndividualInterface.INFECT_S
                && (infectStat[(s + 1) % 2][siteNonG] == GonorrhoeaSiteInfection.STATUS_ASY
                || infectStat[(s + 1) % 2][siteNonG] == GonorrhoeaSiteInfection.STATUS_SYM)
                && (actType != ACT_ANAL || ((s + 1) % 2) != 0)) { // Male's anal cannot transmit to female's genitial

            // From altSite to G
            if (person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G]
                    != ((double[][]) getFields()[FIELDS_SUSCEPT])[RelationshipPerson_MSM.SITE_G
                    + (s == 0 ? 0 : getInfList().length / 2)][0]) {
                person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G]
                        = ((double[][]) getFields()[FIELDS_SUSCEPT])[RelationshipPerson_MSM.SITE_G + (s == 0 ? 0 : getInfList().length / 2)][0];
            }

            tranSucProb = person[(s + 1) % 2].getProbTransBySite()[siteNonG]
                    * person[s].getProbSusBySite()[RelationshipPerson_MSM.SITE_G];

            boolean altSite2G = tranSucProb >= 1;

            if (!altSite2G && tranSucProb > 0) {
                altSite2G = ((GonorrhoeaSiteInfection) getInfList()[RelationshipPerson_MSM.SITE_G]).getRNG().nextDouble()
                        < tranSucProb;
            }

            if (altSite2G) {
                person[s].setLastActInfectious(RelationshipPerson_MSM.SITE_G, true);
                if (person[s] instanceof SiteSpecificTransmissivity) {
                    ((SiteSpecificTransmissivity) person[s]).getLastActStainsAtSite()[RelationshipPerson_MSM.SITE_G]
                            = ((SiteSpecificTransmissivity) person[(s + 1) % 2]).getCurrentStrainsAtSite()[siteNonG];
                }
            }

        }
    }

    @Override
    protected AbstractInfection[] getIndivdualInfectionList(AbstractIndividualInterface person) {
        AbstractInfection[] res;

        int[] range = new int[]{0, getInfList().length / 2};

        if (person.getId() % 2 == 1) {
            range[0] += getInfList().length / 2;
            range[1] += getInfList().length / 2;
        }

        res = Arrays.copyOfRange(getInfList(), range[0], range[1]);

        return res;
    }

    @Override
    public void setAvailablity(AbstractAvailability[] availablity) {
        // Overwrite availability with a Bipartite version
        for (int i = 0; i < availablity.length; i++) {
            availablity[i] = new MSMAvailbility_Bipartite(getRNG());
        }
        super.setAvailablity(availablity);
    }

    @Override
    public void initialiseInfection(long seed) {
        AbstractInfection[] infList = new AbstractInfection[3 * 2]; // 3 sites, 2 genders   
        random.RandomGenerator infRNG;
        if (seed != 0) {
            getFields()[FIELDS_INF_RNG] = new random.RandomEngine[infList.length];
            infRNG = new random.MersenneTwisterFastEngine(seed);
        } else {
            infRNG = null;
        }

        for (int i = 0; i < infList.length; i++) {
            boolean male = i < infList.length / 2;
            double[][] siteParam = new double[GonorrhoeaSiteInfection.DIST_TOTAL][];
            AbstractRealDistribution[] distributions = new AbstractRealDistribution[GonorrhoeaSiteInfection.DIST_TOTAL];

            if (infRNG == null) {
                infRNG = ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[i];
            }
            ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[i] = infRNG;

            siteParam[GonorrhoeaSiteInfection.DIST_EXPOSED_DUR_INDEX] = new double[]{4, 0};
            siteParam[GonorrhoeaSiteInfection.DIST_IMMUNE_DUR_INDEX] = new double[]{7, 0};
            siteParam[GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX] = new double[]{1, 0};

            if (getFields()[FIELDS_TRANSMIT] != null) {
                siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = ((double[][]) getFields()[FIELDS_TRANSMIT])[i];
            } else {
                siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = DEFAULT_TRANSMIT[ i - (male ? 0 : DEFAULT_TRANSMIT.length / 2)];
            }
            if (getFields()[FIELDS_SUSCEPT] != null) {
                siteParam[GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX] = ((double[][]) getFields()[FIELDS_SUSCEPT])[i];
            } else {
                siteParam[GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX] = DEFAULT_SUSCEPT[ i - (male ? 0 : DEFAULT_SUSCEPT.length / 2)];
            }

            double[] var;

            switch (i % (infList.length / 2)) {
                case RelationshipPerson_MSM.SITE_G:
                    // Assume all have sym for male (assumption), and 0.4 for female (Johnson)                    

                    if (getFields()[HETRO_POP_TREAT_RATE_GENITAL] != null) {
                        siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX]
                                = ((double[][]) getFields()[HETRO_POP_TREAT_RATE_GENITAL])[male ? 0 : 1];

                    } else {
                        if (male) {
                            siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{1, 0};
                        } else {
                            siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{0.4, 0};
                        }
                    }

                    // Kit email 20130911 for urethal
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new double[]{2.5714, 2.24343};

                    if (getFields()[HETRO_POP_SYM_DURATION_ADJ] != null
                            && ((float[][]) getFields()[HETRO_POP_SYM_DURATION_ADJ])[male ? 0 : 1][0] != 1) {
                        siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX][0]
                                *= ((float[][]) getFields()[HETRO_POP_SYM_DURATION_ADJ])[male ? 0 : 1][0];
                    }

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]
                            = new GammaDistribution(infRNG,var[0], 1/var[1]);

                    // Johnson                  
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{15 * 7, 5 * 7};
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]
                            = new GammaDistribution(infRNG,var[0], 1/var[1]);
                    break;

                case RelationshipPerson_MSM.SITE_A:
                    // Anal : Bissessor 2011: 7 out 47 of has proctitis
                    siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{7.0 / 47, 0};

                    // Assume same as Urethal
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new double[]{2.5714, 2.24343};

                    if (getFields()[HETRO_POP_SYM_DURATION_ADJ] != null
                            && ((float[][]) getFields()[HETRO_POP_SYM_DURATION_ADJ])[male ? 0 : 1][1] != 1) {
                        siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX][0]
                                *= ((float[][]) getFields()[HETRO_POP_SYM_DURATION_ADJ])[male ? 0 : 1][1];
                    }

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]
                            = new GammaDistribution(infRNG,var[0], 1/var[1]);

                    // From C. Fairley email at 20131120
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{12 * 30, 5 * 7};

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]
                            = new GammaDistribution(infRNG,var[0], 1/var[1]);

                    break;
                case RelationshipPerson_MSM.SITE_R:
                    // Assume none have sym
                    siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{0, 0};

                    /*
                     if(getFields()[HETRO_POP_SYM_DURATION_ADJ] != null 
                     && ((float[][]) getFields()[HETRO_POP_SYM_DURATION_ADJ])[male? 0:1][2] != 1){                        
                     siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX][0] 
                     *= ((float[][]) getFields()[HETRO_POP_SYM_DURATION_ADJ])[male? 0:1][2];                        
                     }    
                     */
                    // From C. Fairley email at 20131121, Fairley et al 2011
                    siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{12 * 7, 3 * 7};

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]
                            = new GammaDistribution(infRNG,var[0], 1/var[1]);
                    break;
            }

            infList[i] = new GonorrhoeaSiteInfection(infRNG, i % (infList.length / 2),
                    siteParam, distributions);
        }

        setInfList(infList);
    }

}
