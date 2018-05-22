package population;

import java.util.Arrays;
import person.AbstractIndividualInterface;
import population.person.RelationshipPerson_MSM;
import population.relationshipMap.RegCasRelationship;

/**
 * A MSM population with strains markers
 *
 * @author Ben Hui
 * @version 20140807
 */
public class MSMPopulation_SingleType_Strains extends MSMPopulation {

    public static final int MSM_POP_STS_INDEX_OFFSET = LENGTH_FIELDS_MSM_POP;
    public static final int MSM_POP_STS_STRAIN_ADJ = MSM_POP_STS_INDEX_OFFSET;
    public static final int LENGTH_FIELDS_MSM_POP_STS = MSM_POP_STS_STRAIN_ADJ + 1;

    public static final int INEXT_TO_MIN_STRAIN = 0;
    public static final int INDEX_TO_MAX_STRAIN = 1;

    public static final Object[] DEFAULT_MSM_POP_STS_FIELDS = {
        new double[][]{{0, 0}, {0, 0}, {0, 0}} // Cumulative probability of switching to lower, higher strain (i.e.  x < ent[0] = lower, ent[0]< x < ent[1] = higher)                                 
    };

    public MSMPopulation_SingleType_Strains(long seed) {
        super(seed);
        Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_MSM_POP_STS);
        for (int i = LENGTH_FIELDS_MSM_POP; i < newFields.length; i++) {
            newFields[i] = DEFAULT_MSM_POP_STS_FIELDS[i - LENGTH_FIELDS_MSM_POP];
        }
        super.setFields(newFields);
    }

    @Override
    protected boolean[][] performAct(RegCasRelationship rel) {
        boolean[][] hasUnprotectedSex = super.performAct(rel); // hasUnprotectedSex       
        RelationshipPerson_MSM[] person = new RelationshipPerson_MSM[rel.getLinks().length];
        int[][] infectStat = new int[rel.getLinks().length][];
        int[][] infectStrainType = new int[rel.getLinks().length][];

        for (int p = 0; p < rel.getLinks().length; p++) {
            person[p] = (RelationshipPerson_MSM) getLocalData().get(rel.getLinks()[p]);
            infectStat[p] = person[p].getInfectionStatus();
            infectStrainType[p] = person[p].getCurrentStrainsAtSite();

        }
        // Check for tranmission across double infection 
        for (int a = 0; a < hasUnprotectedSex.length; a++) {
            if (hasUnprotectedSex[a][0]) {
                boolean infectedWithdifferentStrains;
                int nonGTarget = -1;

                switch (a) {
                    case ACT_ANAL:
                        nonGTarget = RelationshipPerson_MSM.SITE_A;
                        break;
                    case ACT_ORAL:
                        nonGTarget = RelationshipPerson_MSM.SITE_R;
                        break;
                    case ACT_RIMMING:
                        throw new UnsupportedOperationException(getClass().getName() + ".performAct(): Rimming not supported in the current version");
                    default:
                        throw new UnsupportedOperationException(getClass().getName() + ".performAct(): Act type not supported in the current version");
                }

                // Infection found at already infected sites               
                // with different strains
                infectedWithdifferentStrains
                        = (infectStat[0][RelationshipPerson_MSM.SITE_G] != AbstractIndividualInterface.INFECT_S
                        && infectStat[1][nonGTarget] != AbstractIndividualInterface.INFECT_S
                        && infectStrainType[0][RelationshipPerson_MSM.SITE_G] != infectStrainType[1][nonGTarget])
                        || (infectStat[0][nonGTarget] != AbstractIndividualInterface.INFECT_S
                        && infectStat[1][RelationshipPerson_MSM.SITE_G] != AbstractIndividualInterface.INFECT_S
                        && infectStrainType[0][nonGTarget] != infectStrainType[1][RelationshipPerson_MSM.SITE_G]);

                if (infectedWithdifferentStrains) {

                    switchStrains(person, infectStat, infectStrainType, nonGTarget);

                }

            }

        }

        return hasUnprotectedSex;
    }
    
    /**
     * Handle switching of strains 
     * @param person person involved
     * @param infectStat infection type matrix
     * @param infectStrainType strain type matrix
     * @param nonGTarget Non-genital target
     */

    protected void switchStrains(RelationshipPerson_MSM[] person, int[][] infectStat, int[][] infectStrainType, int nonGTarget) {
        random.RandomGenerator infRNG;
        double pSwitch;
        double[][] switchProb = (double[][]) getFields()[MSM_POP_STS_STRAIN_ADJ];
        
        int[] strainId = new int[2]; // lowest, highest
        strainId[0] = Integer.MAX_VALUE;
        strainId[1] = Integer.MIN_VALUE;
        
        for (int p = 0; p < person.length; p++) {
            if (infectStat[p][RelationshipPerson_MSM.SITE_G] != AbstractIndividualInterface.INFECT_S) {
                strainId[0] = Math.min(strainId[0], infectStrainType[p][RelationshipPerson_MSM.SITE_G]);
                strainId[1] = Math.max(strainId[1], infectStrainType[p][RelationshipPerson_MSM.SITE_G]);
            }
            if (infectStat[p][nonGTarget] != AbstractIndividualInterface.INFECT_S) {
                strainId[0] = Math.min(strainId[0], infectStrainType[p][nonGTarget]);
                strainId[1] = Math.max(strainId[1], infectStrainType[p][nonGTarget]);
            }
        }
        
        if (switchProb[RelationshipPerson_MSM.SITE_G][MSMPopulation_SingleType_Strains.INDEX_TO_MAX_STRAIN] > 0) {
            infRNG = ((random.RandomGenerator[]) getFields()[FIELDS_INF_RNG])[RelationshipPerson_MSM.SITE_G];
            
            pSwitch = infRNG.nextDouble();
            
            if (pSwitch < switchProb[RelationshipPerson_MSM.SITE_G][MSMPopulation_SingleType_Strains.INEXT_TO_MIN_STRAIN]) {
                for (RelationshipPerson_MSM per : person) {
                    if (per.getInfectionStatus()[RelationshipPerson_MSM.SITE_G] != AbstractIndividualInterface.INFECT_S) {
                        per.getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_G] = strainId[0];
                    }
                }
            } else if (pSwitch < switchProb[RelationshipPerson_MSM.SITE_G][MSMPopulation_SingleType_Strains.INDEX_TO_MAX_STRAIN]) {
                for (RelationshipPerson_MSM per : person) {
                    if (per.getInfectionStatus()[RelationshipPerson_MSM.SITE_G] != AbstractIndividualInterface.INFECT_S) {
                        per.getCurrentStrainsAtSite()[RelationshipPerson_MSM.SITE_G] = strainId[1];
                    }
                }
            }
        }
        
        if (nonGTarget != -1 && switchProb[nonGTarget][MSMPopulation_SingleType_Strains.INDEX_TO_MAX_STRAIN] > 0) {
            infRNG = ((random.RandomEngine[]) getFields()[FIELDS_INF_RNG])[nonGTarget];
            pSwitch = infRNG.nextDouble();
            if (pSwitch < switchProb[nonGTarget][MSMPopulation_SingleType_Strains.INEXT_TO_MIN_STRAIN]) {
                for (RelationshipPerson_MSM per : person) {
                    if (per.getInfectionStatus()[nonGTarget] != AbstractIndividualInterface.INFECT_S) {
                        per.getCurrentStrainsAtSite()[nonGTarget] = strainId[0];
                    }
                }
            } else if (pSwitch < switchProb[nonGTarget][MSMPopulation_SingleType_Strains.INDEX_TO_MAX_STRAIN]) {
                for (RelationshipPerson_MSM per : person) {
                    if (per.getInfectionStatus()[nonGTarget] != AbstractIndividualInterface.INFECT_S) {
                        per.getCurrentStrainsAtSite()[nonGTarget] = strainId[1];
                    }
                }
            }
        }
    }

}
