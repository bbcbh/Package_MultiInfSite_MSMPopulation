package population;


import infection.AbstractInfection;
import infection.GonorrhoeaSiteInfection;
import java.util.Arrays;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import static population.AbstractRegCasRelMapPopulation.*;
import person.AbstractIndividualInterface;
import population.person.RelationshipPerson_MSM;
import population.person.SiteSpecificTransmissivity;
import population.relationshipMap.RegCasRelationship;
import random.MersenneTwisterFastEngine;
import util.PersonClassifier;
import util.StaticMethods;

/**
 * A MSM population with strains markers for multiple strain infection. Marker will be identified by binary representation of (strain marker)
 *
 * e.g. Strain type of
 * <ul>
 * <li>-1 or 0 = No strain </li><li>1 = Background strains</li><li>2 = New strains only </li><li>3 = Co-infection</li>
 * </ul>
 *
 *
 * @author Ben Hui
 * @version 20150414
 *
 * History:
 * <p>
 * 20150414: Add negative option for strain mixing.</p>
 * <p>
 * 20150402: Added treatment modification setting and associated method </p>
 * <p>
 * 20140922: Updated setInstantInfection method </p>
 * <p>
 * 20140219: Add setting for strain mixing </p>
 *
 *
 */
public class MSMPopulation_MultiType_Strains extends MSMPopulation {

    public static final int MSM_POP_MTS_INDEX_OFFSET = LENGTH_FIELDS_MSM_POP;
    public static final int MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS = MSM_POP_MTS_INDEX_OFFSET;
    public static final int MSM_POP_MTS_STRAIN_PRIORITY_DEFAULT = MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS + 1;
    public static final int MSM_POP_MTS_STRAIN_PROP = MSM_POP_MTS_STRAIN_PRIORITY_DEFAULT + 1;
    public static final int MSM_POP_MTS_STRAIN_MIXING = MSM_POP_MTS_STRAIN_PROP + 1;
    public static final int MSM_POP_MTS_STRAIN_TREATMENT_CONVERT = MSM_POP_MTS_STRAIN_MIXING + 1;
    public static final int LENGTH_FIELDS_MSM_POP_MTS = MSM_POP_MTS_STRAIN_TREATMENT_CONVERT + 1;

    private final int[] incidentCount; // { by transmission : G01, G10, G11, A01, A10, ... , by mixing: G01, G10, G11, A01, A10, ...}

    public static final Object[] DEFAULT_MSM_POP_MTS_FIELDS = {
        // MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS
        // Number of strains included, with number of strains type = 2^number of strain -1        
        2,
        // MSM_POP_MTS_STRAIN_PRIORITY_DEFAULT
        // Bit position by increasing dominant. e.g.   new int[]{0, 1} = 2nd (strainPt.e. last) bit dominant  
        new int[]{0, 1},
        // MSM_POP_MTS_STRAIN_PROP        
        // double[strain][site]{siteParam_id, param_length, parameter values...}
        /* 
         Example:
         new double[][]{
         // RelationshipPerson_MSM.SITE_G
         {},
         // RelationshipPerson_MSM.SITE_A:
         {GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX, 2, 6 * 30, 5 * 7}, // Only half as long than the main strans
         // RelationshipPerson_MSM.SITE_R:
         {GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX, 2, 6 * 7, 3 * 7}, // Only half as long than the main strans
         }
         */
        new double[][][]{
            {}, // Default for strain 0            
            {}, //            
        },
        // MSM_POP_MTS_STRAIN_MIXING        
        // float[combined_strain][siteId][ transfromed_strain, likelihood ...] 
        new float[][][]{},
        // MSM_POP_MTS_STRAIN_TREATMENT_CONVERT
        // float[siteId][current_strain_type, transformed_strain, likelihood...]
        // e.g. 100% failure for last bit strain (strainPt.e. A10 -> A10, A11 -> A10) will 
        // have entry float[MSM.SITE_ANAL] = {0b10, 0b10, 1, 0b11, 0b10, 1}
        new float[][]{}

    };

    @Override
    protected void treatPerson(RelationshipPerson_MSM msm) {
        // Implementation treatment failure

        float[][] treatConvertOptionAll = (float[][]) getFields()[MSM_POP_MTS_STRAIN_TREATMENT_CONVERT];

        if (treatConvertOptionAll == null || treatConvertOptionAll.length == 0) {
            super.treatPerson(msm);
        } else {
            for (int site = 0; site < msm.getInfectionStatus().length; site++) {

                // Sucuess treatment (default)                
                int newInfectStatus = AbstractIndividualInterface.INFECT_S;
                int newStrain = SiteSpecificTransmissivity.STRAIN_NONE;
                int currentStrain = msm.getCurrentStrainsAtSite()[site];

                double timeUntilNextStage = Double.POSITIVE_INFINITY;

                if (site < treatConvertOptionAll.length) {
                    float[] treatConvertOptions = treatConvertOptionAll[site];

                    if (treatConvertOptions != null && treatConvertOptions.length > 0) {
                        int cSPt = 0;
                        boolean converted = false;

                        while (cSPt < treatConvertOptions.length && !converted) {
                            if (treatConvertOptions[cSPt] == currentStrain) {
                                float likeihood = treatConvertOptions[cSPt + 2];
                                converted = getInfList()[site].getRNG().nextFloat() < likeihood;
                                if (converted) {
                                    newInfectStatus = msm.getInfectionStatus()[site];
                                    newStrain = (int) treatConvertOptions[cSPt + 1];                                   
                                    int[] strainPriority = msm.getStrainPriorityBitAtSite()[site];
                                    
                                    for(int strainPt = 0; strainPt < strainPriority.length; strainPt++){
                                        double strainLasted = msm.getCurrentStrainLastUntilOfAgeAtSite()[site][strainPriority[strainPt]];
                                        
                                        if(!Double.isNaN(strainLasted)){
                                            timeUntilNextStage = Math.min(timeUntilNextStage,
                                                strainLasted - msm.getAge());
                                        }                                        
                                    }                                                                       
                                 
                                }
                            }
                            cSPt += 3;
                        }
                    }
                }

                msm.getInfectionStatus()[site] = newInfectStatus;
                msm.getCurrentStrainsAtSite()[site] = newStrain;
                msm.setTimeUntilNextStage(site, timeUntilNextStage);

                if (newInfectStatus == AbstractIndividualInterface.INFECT_S) {                    // Sucess treatment 
                    Arrays.fill(msm.getCurrentStrainLastUntilOfAgeAtSite()[site], Double.NaN);
                }
            }
        }

    }

    public MSMPopulation_MultiType_Strains(long seed) {
        super(seed);
        Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_MSM_POP_MTS);
        setFields(newFields);
        for (int f = MSM_POP_MTS_INDEX_OFFSET; f < newFields.length; f++) {
            getFields()[f] = DEFAULT_MSM_POP_MTS_FIELDS[f - MSM_POP_MTS_INDEX_OFFSET];

        }

        incidentCount = new int[3 * 2 * ((1 << (int) getFields()[MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS]) - 1)];
    }

    // StrainType = 1 << strainBit 
    @Override
    public void introduceStrains(int siteId, PersonClassifier infClassifer, float[] numIntroduced, int strainBit) {

        int[] numForEachStrainByClassifier = new int[infClassifer.numClass()];
        int[] numInfTotalByClassifier = new int[infClassifer.numClass()];

        int infId = (getInfList().length / ((Integer) getFields()[MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS])) * (strainBit) + siteId;

        for (AbstractIndividualInterface indiv : getPop()) {
            int cI = infClassifer.classifyPerson(indiv);
            if (cI >= 0) {
                numInfTotalByClassifier[cI]++;
            }
        }

        for (int c = 0; c < infClassifer.numClass(); c++) {
            numForEachStrainByClassifier[c] = Math.round(numIntroduced[c] >= 1
                    ? numIntroduced[c] : numInfTotalByClassifier[c] * numIntroduced[c]);
        }

        for (AbstractIndividualInterface person : getPop()) {
            int cI = infClassifer.classifyPerson(person);
            if (cI >= 0) {
                RelationshipPerson_MSM msm = (RelationshipPerson_MSM) person;
                random.RandomGenerator infRNG = getInfList()[infId].getRNG();
                if (infRNG.nextInt(numInfTotalByClassifier[cI]) < numForEachStrainByClassifier[cI]) {
                    if (person.getInfectionStatus()[siteId] == AbstractIndividualInterface.INFECT_S) {
                        getInfList()[infId].infecting(msm);
                        msm.getCurrentStrainsAtSite()[siteId] = 1 << (strainBit);
                        msm.getCurrentStrainLastUntilOfAgeAtSite()[siteId][strainBit] = msm.getAge() + msm.getTimeUntilNextStage(siteId);
                    } else {
                        strainMixing(msm, siteId, (1 << strainBit));
                    }

                    numForEachStrainByClassifier[cI]--;
                }
                numInfTotalByClassifier[cI]--;
            }
        }
    }

    @Override
    protected boolean[][] performAct(RegCasRelationship rel) {
        RelationshipPerson_MSM[] person = new RelationshipPerson_MSM[rel.getLinks().length];
        int[][] preStrainType = new int[2][];
        for (int p = 0; p < rel.getLinks().length; p++) {
            person[p] = (RelationshipPerson_MSM) getLocalData().get(rel.getLinks()[p]);
            preStrainType[p] = Arrays.copyOf(person[p].getCurrentStrainsAtSite(), person[p].getCurrentStrainsAtSite().length);
        }

        boolean[][] hasUnprotectedSex = super.performAct(rel); // hasUnprotectedSex[actType]{occured, from_genital_person_1, from_genital_person_2}

        for (int a = 0; a < hasUnprotectedSex.length; a++) {

            if (hasUnprotectedSex[a][0]) {
                int nonG = -1;
                switch (a) {
                    case ACT_ANAL:
                        nonG = RelationshipPerson_MSM.SITE_A;
                        break;
                    case ACT_ORAL:
                        nonG = RelationshipPerson_MSM.SITE_R;
                        break;
                    default:
                        throw new UnsupportedOperationException("Non anal or oral acts not supported yet.");
                }

                for (int p = 0; p < person.length; p++) {
                    // Full transfer with no mulation at this stage 
                    int target;
                    if (hasUnprotectedSex[a][p + 1]) {
                        target = RelationshipPerson_MSM.SITE_G;
                        strainMixing(person[p], target, person[(p + 1) % person.length].getCurrentStrainsAtSite()[nonG]);
                    } else {
                        target = nonG;
                        strainMixing(person[p], target, person[(p + 1) % person.length].getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_G]);
                    }

                    int newStrain = person[p].getCurrentStrainsAtSite()[target];

                    // Check for gaining of strains
                    if (person[p].getLastActInfectious()[target]) {
                        newStrain = person[p].getLastActStainsAtSite()[target];
                    }

                    if (newStrain > 0 && newStrain != preStrainType[p][target]) {
                        if (person[p].getLastActInfectious()[target]) {
                            // Transmittion to a susceptible
                            incidentCount[3 * target + (newStrain - 1)]++;
                        } else {
                            // Transmission to an already infected
                            incidentCount[3 * ((1 << (int) getFields()[MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS]) - 1)
                                    + 3 * target + (newStrain - 1)]++;
                        }

                    }

                }

            }
        }

        return hasUnprotectedSex;
    }

    @Override
    public int[] generatedPopSnapCount() {
        int[] incidentTotal = super.generatedPopSnapCount();
        int orgLen = incidentTotal.length;
        incidentTotal = Arrays.copyOf(incidentTotal, orgLen + incidentCount.length);
        System.arraycopy(incidentCount, 0, incidentTotal, orgLen, incidentCount.length);

        return incidentTotal;
    }

    private void strainMixing(RelationshipPerson_MSM msm, int siteIndex, int newStrain) {
        if (newStrain > 0 && newStrain != msm.getCurrentStrainsAtSite()[siteIndex]
                && msm.getCurrentStrainsAtSite()[siteIndex] > 0) {
            
            int orgStrain = msm.getCurrentStrainsAtSite()[siteIndex];
            // Combine all   
            msm.getCurrentStrainsAtSite()[siteIndex] |= newStrain;

            float[][][] strainMixingByCombineStrain = (float[][][]) getFields()[MSM_POP_MTS_STRAIN_MIXING];
            if (strainMixingByCombineStrain != null && msm.getCurrentStrainsAtSite()[siteIndex] < strainMixingByCombineStrain.length) {
                float[][] strainMixingBySite = strainMixingByCombineStrain[msm.getCurrentStrainsAtSite()[siteIndex]];
                if (siteIndex < strainMixingBySite.length) {
                    float[] strainMixingEntry = strainMixingBySite[siteIndex];
                    if (strainMixingEntry.length > 0) {
                        float pConvert = getIndivdualInfectionList(msm)[siteIndex].getRNG().nextFloat();
                        int pt = 0;
                        while ((pt + 1) < strainMixingEntry.length && pConvert >= strainMixingEntry[pt + 1]) {
                            pt += 2;
                        }
                        int strainChoice = (int) strainMixingEntry[pt];                        
                        
                        if(strainChoice < 0) { // Unmodified
                            strainChoice = orgStrain;
                        }                                                                     
                        msm.getCurrentStrainsAtSite()[siteIndex] = strainChoice;

                    }
                }
            }

            int domainantStrainBit = msm.getCurrentDomainantStrainAtSite(siteIndex);

            double dominantStrainLasted = msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][domainantStrainBit];

            for (int stPt = 0; stPt < msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex].length; stPt++) {
                if (stPt != domainantStrainBit) {
                    long strainAloneLength;
                    String paramStr;

                    // Assume to be the same as dominant strain initially
                    msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][stPt] = dominantStrainLasted;
                    // Determine possible length for non-domainant strian                                         
                    AbstractInfection strain_inf = getInfList()[stPt * msm.getInfectionStatus().length + siteIndex];

                    // Check if new mix strain will lead to symptom 
                    String paramStrSym = AbstractInfection.PARAM_DIST_NEXT_VALUE_REGEX.replaceFirst("\\d+",
                            Integer.toString(GonorrhoeaSiteInfection.DIST_SYM_INDEX));

                    double pSym = (Double) strain_inf.getParameter(paramStrSym);

                    boolean nonDomainantSym = pSym >= 1;

                    if (!nonDomainantSym && pSym > 0) {
                        nonDomainantSym = strain_inf.getRNG().nextDouble() < pSym;
                    }

                    if (nonDomainantSym) {
                        paramStr = AbstractInfection.PARAM_DIST_NEXT_VALUE_REGEX.replaceFirst("\\d+",
                                Integer.toString(GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX));
                    } else {
                        paramStr = AbstractInfection.PARAM_DIST_NEXT_VALUE_REGEX.replaceFirst("\\d+",
                                Integer.toString(GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX));
                    }

                    strainAloneLength = Math.round((Double) strain_inf.getParameter(paramStr));

                    // Non domainant strain lead to symptom 
                    if (msm.getInfectionStatus()[siteIndex] == GonorrhoeaSiteInfection.STATUS_ASY
                            && nonDomainantSym) {
                        msm.getInfectionStatus()[siteIndex] = GonorrhoeaSiteInfection.STATUS_SYM;

                        if (msm.getTimeUntilNextStage(siteIndex) > strainAloneLength) {

                            msm.setTimeUntilNextStage(siteIndex, strainAloneLength);
                            msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][domainantStrainBit] = strainAloneLength + msm.getAge();
                            dominantStrainLasted = Math.min(strainAloneLength + msm.getAge(), dominantStrainLasted);
                        }
                    }
                    msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][stPt] = Math.min(dominantStrainLasted, msm.getAge() + strainAloneLength);
                }
            }

            // Trim dominant by site
            for (int stPt = 0; stPt < msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex].length; stPt++) {
                msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][stPt] = Math.min(dominantStrainLasted,
                        msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][stPt]);

            }

        }

    }

    @Override
    public int[] setInstantInfection(int infId, PersonClassifier infClassifer, float[] prevalByClass, int preExposeMax, float[][] strainDecompositionByClassStrains) {

        int[] res = super.setInstantInfection(infId, infClassifer, prevalByClass, preExposeMax, strainDecompositionByClassStrains);
        // Initialised strainLastUntilOfAgeAtSite fields
        int numStrainType = 1 << ((Integer) getFields()[MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS]) - 1;
        for (AbstractIndividualInterface person : getPop()) {
            RelationshipPerson_MSM msm = (RelationshipPerson_MSM) person;
            for (int site = 0; site < msm.getInfectionStatus().length; site++) {
                msm.getCurrentStrainLastUntilOfAgeAtSite()[site] = new double[numStrainType];
                Arrays.fill(msm.getCurrentStrainLastUntilOfAgeAtSite()[site], Double.NaN);
                msm.getStrainPriorityBitAtSite()[site] = (int[]) getFields()[MSM_POP_MTS_STRAIN_PRIORITY_DEFAULT];

                int strainType = msm.getCurrentStrainsAtSite()[site];

                if (msm.getInfectionStatus()[site] != AbstractIndividualInterface.INFECT_S) {
                    // Was infected at some stage                    
                    if (strainType == 0) { // Error
                        System.err.println("Warning : Converting strain type from 0 to 1 for infected");
                        msm.getCurrentStrainsAtSite()[site] = 1;
                    }

                    int matcher = 1;
                    int strainPt = 0;
                    while (strainPt < numStrainType) {
                        if ((matcher & strainType) > 0) {
                            msm.getCurrentStrainLastUntilOfAgeAtSite()[site][strainPt]
                                    = msm.getAge() + msm.getTimeUntilNextStage(site);
                        }
                        strainPt++;
                        matcher = matcher << 1;
                    }
                } else if (strainType > 0) {
                    System.err.println("Warning : Converting strain type from "
                            + strainType + " to 0 for non infected");
                    msm.getCurrentStrainsAtSite()[site] = 0;
                }
            }
        }

        return res;
    }

    @Override
    protected AbstractInfection[] incrementPersonStat(AbstractIndividualInterface person, int deltaT) {
        AbstractInfection[] dominatedInfectionList = super.incrementPersonStat(person, deltaT);

        // Check for duration for all non-domainate infection and remove if necessary
        RelationshipPerson_MSM msm = (RelationshipPerson_MSM) person;

        for (int siteIndex = 0; siteIndex < dominatedInfectionList.length; siteIndex++) {
            int currentStrainAtSite = msm.getCurrentStrainsAtSite()[siteIndex];
            if (currentStrainAtSite > 0) {

                int domainantStrainPointer = dominatedInfectionList[siteIndex].getInfectionIndex() / msm.getInfectionStatus().length;
                int nonDominatantStrains = currentStrainAtSite - (1 << domainantStrainPointer);

                int nonDominantStraintPointer = 0;
                while (nonDominatantStrains > 0) {
                    // The removal of lesser strains
                    if (nonDominatantStrains % 2 != 0
                            && msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][nonDominantStraintPointer] <= person.getAge()) {

                        msm.getCurrentStrainLastUntilOfAgeAtSite()[siteIndex][nonDominantStraintPointer] = Double.NaN;
                        msm.getCurrentStrainsAtSite()[siteIndex] = currentStrainAtSite - (1 << nonDominantStraintPointer);
                    }

                    nonDominatantStrains = nonDominatantStrains >> 1;
                    nonDominantStraintPointer++;
                }
            }
        }

        return dominatedInfectionList;
    }

    @Override
    protected AbstractInfection[] getIndivdualInfectionList(AbstractIndividualInterface person) {
        RelationshipPerson_MSM msm = (RelationshipPerson_MSM) person;
        AbstractInfection[] inf = new AbstractInfection[msm.getInfectionStatus().length];
        for (int site = 0; site < msm.getInfectionStatus().length; site++) {
            if (msm.getInfectionStatus()[site] != AbstractIndividualInterface.INFECT_S
                    && msm.getCurrentStrainsAtSite()[site] <= 0) {
                System.err.println("getIndivdualInfectionList: Warning, Strain type <= 0 and infectious! Switch to 1 instead");
                msm.getCurrentStrainsAtSite()[site] = 1;
            }
            inf[site] = getInfList()[inf.length * msm.getCurrentDomainantStrainAtSite(site) + site];
        }

        return inf;

    }

    @Override
    public void initialiseInfection(long seed) {
        int numSite = 3;
        AbstractInfection[] infList = new AbstractInfection[numSite
                * (Integer) getFields()[MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS]]; // 3 sites  * Number of strain 

        random.RandomGenerator infRNG;

        if (seed != 0) {
            getFields()[FIELDS_INF_RNG] = new random.RandomGenerator[infList.length];
            infRNG = new MersenneTwisterFastEngine(seed);
        } else {
            infRNG = null;
        }

        for (int infListId = 0; infListId < infList.length; infListId++) {

            double[][] siteParam = new double[GonorrhoeaSiteInfection.DIST_TOTAL][];
            AbstractRealDistribution[] distributions = new AbstractRealDistribution[GonorrhoeaSiteInfection.DIST_TOTAL];

            if (infRNG == null) {
                infRNG = ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[infListId]; // Imported
            }

            ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[infListId] = infRNG;

            if (getFields()[MSM_POP_MTS_STRAIN_PROP] != null) {
                loadPredefinedStrainProp(infListId, numSite, siteParam);
            }

            // Backward comp.
            if (getFields()[FIELDS_TRANSMIT] != null && infListId < ((double[][]) getFields()[FIELDS_TRANSMIT]).length) {
                siteParam[GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX] = ((double[][]) getFields()[FIELDS_TRANSMIT])[infListId];
            }
            if (getFields()[FIELDS_SUSCEPT] != null && infListId < ((double[][]) getFields()[FIELDS_SUSCEPT]).length) {
                siteParam[GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX] = ((double[][]) getFields()[FIELDS_SUSCEPT])[infListId];
            }

            //siteParam[GonorrhoeaSiteInfection.DIST_EXPOSED_DUR_INDEX] = new double[]{4, 0};
            //siteParam[GonorrhoeaSiteInfection.DIST_IMMUNE_DUR_INDEX] = new double[]{7, 0};
            replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_EXPOSED_DUR_INDEX, new double[]{4, 0});
            replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_IMMUNE_DUR_INDEX, new double[]{7, 0});

            double[] var;
            // Same for every strain atm
            switch (infListId % numSite) {
                case RelationshipPerson_MSM.SITE_G:

                    // Assume all have sym
                    // siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{1, 0};
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_SYM_INDEX, new double[]{0.92, 0}); // Johnson = 0.9

                    // Kit email 20130911 for urethal
                    // siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new double[]{2.5714, 2.24343};                    
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX, new double[]{2.5714, 2.24343});

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new GammaDistribution(infRNG, var[0], 1/var[1]);

                    // Assume duration same as Garnett 1999, Brunham 1991, SD same as Johnson 2010
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX, new double[]{185, 5 * 7});
                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new GammaDistribution(infRNG, var[0], 1/var[1]);

                    // Transmission and susceptiblity - from MSM paper                   
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX, new double[]{1, 0});
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX, new double[]{1, 0});

                    break;
                case RelationshipPerson_MSM.SITE_A:

                    // Anal : Bissessor 2011: 7 out 47 of has proctitis
                    //siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{7.0 / 47, 0};
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_SYM_INDEX, new double[]{7.0 / 47, 0});

                    // Assume same as Urethal
                    //siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new double[]{2.5714, 2.24343};
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX, new double[]{2.5714, 2.24343});

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_SYM_INDEX] = new GammaDistribution(infRNG, var[0], 1/var[1]);

                    // From C. Fairley email at 20131120                    
                    //siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{12 * 30, 5 * 7};                    
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX, new double[]{12 * 30, 5 * 7});

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new GammaDistribution(infRNG, var[0], 1/var[1]);

                    // Transmission and susceptiblity - from MSM paper                                       
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX, new double[]{0.024258237711416882, 0});
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX, new double[]{0.8402332453833801, 0});

                    break;

                case RelationshipPerson_MSM.SITE_R:
                    // Assume none have sym
                    //siteParam[GonorrhoeaSiteInfection.DIST_SYM_INDEX] = new double[]{0, 0};
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_SYM_INDEX, new double[]{0, 0});

                    // From C. Fairley email at 20131121, Fairley et al 2011
                    //siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new double[]{12 * 7, 3 * 7};                    
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX, new double[]{12 * 7, 3 * 7});

                    var = StaticMethods.generatedGammaParam(siteParam[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX]);
                    distributions[GonorrhoeaSiteInfection.DIST_INFECT_DUR_ASY_INDEX] = new GammaDistribution(infRNG, var[0], 1/var[1]);

                    // Transmission and susceptiblity - from MSM paper                                       
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_TRANS_PROB_INDEX, new double[]{0.08654020531624781, 0});
                    replaceIfNull(siteParam, GonorrhoeaSiteInfection.DIST_SUS_PROB_INDEX, new double[]{0.6288277879484502, 0});

                    break;
            }

            infList[infListId] = new GonorrhoeaSiteInfection(infRNG, infListId, infListId % numSite, siteParam, distributions);
        }
        setInfList(infList);
    }

    private void replaceIfNull(double[][] siteParam, int site_param_id, double[] ent) {
        if (siteParam[site_param_id] == null) {
            siteParam[site_param_id] = ent;
        }
    }

    private void loadPredefinedStrainProp(int infListId, int numSite, double[][] siteParam) {
        int strainId = infListId / numSite;
        int siteId = infListId % numSite;
        double[][] byStrain = ((double[][][]) getFields()[MSM_POP_MTS_STRAIN_PROP])[strainId];
        if (byStrain != null && byStrain.length == numSite) { // 3 site
            double[] bySite = byStrain[siteId];
            if (bySite != null) {
                int pt = 0;
                while (pt < bySite.length) {
                    int siteParam_id = (int) bySite[pt];
                    pt++;
                    int siteParam_length = (int) bySite[pt];
                    pt++;
                    double[] param = Arrays.copyOfRange(bySite, pt, pt + siteParam_length);
                    if (siteParam_id > 0 && siteParam_length != 0) {
                        siteParam[siteParam_id] = param;
                    }
                    pt += siteParam_length;
                }
            }
        }
    }

    @Override
    protected AbstractIndividualInterface generateNewPerson(int nextId, AbstractIndividualInterface p, double newAge) {
        RelationshipPerson_MSM msm = (RelationshipPerson_MSM) super.generateNewPerson(nextId, p, newAge);
        int numStrainType = 1 << ((Integer) getFields()[MSM_POP_MTS_STRAIN_MAX_NUM_STRAINS]) - 1;
        for (int site = 0; site < msm.getInfectionStatus().length; site++) {
            msm.getCurrentStrainLastUntilOfAgeAtSite()[site] = new double[numStrainType];
            Arrays.fill(msm.getCurrentStrainLastUntilOfAgeAtSite()[site], Double.NaN);
            msm.getStrainPriorityBitAtSite()[site] = (int[]) getFields()[MSM_POP_MTS_STRAIN_PRIORITY_DEFAULT];
        }

        return msm;
    }

}
