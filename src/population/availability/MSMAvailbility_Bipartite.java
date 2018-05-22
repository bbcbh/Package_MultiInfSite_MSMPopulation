package population.availability;

import java.util.Arrays;
import person.AbstractIndividualInterface;
import random.RandomGenerator;


/**
 *
 * @author Ben Hui
 * @version 20180521
 */
public class MSMAvailbility_Bipartite extends MSMAvailablity {

    public MSMAvailbility_Bipartite(RandomGenerator RNG) {
        super(RNG);
    }

    // Bipartite by "gender"
    @Override
    public void setAvailablePopulation(AbstractIndividualInterface[][] available) {
        AbstractIndividualInterface[][] bipart = new AbstractIndividualInterface[2][available[0].length];
        int[] counter = new int[2];
        for (AbstractIndividualInterface item : available[0]) {
            int gI = (item.getId() % 2 == 0) ? 0 : 1;
            bipart[gI][counter[gI]] = item;
            counter[gI]++;
        }
        for (int i = 0; i < bipart.length; i++) {
            bipart[i] = Arrays.copyOf(bipart[i], counter[i]);
        }
        super.setAvailablePopulation(bipart);
    }

    @Override
    public int generatePairing() {
        int numPairing = Integer.MAX_VALUE;
        
        for (AbstractIndividualInterface[] available1 : available) {
            numPairing = Math.min(numPairing, available1.length);
            if (available1.length > 1) {
                util.ArrayUtilsRandomGenerator.shuffleArray(available1, getRNG());
            }
        }

        pairing = new AbstractIndividualInterface[numPairing][2];

        for (int p = 0; p < numPairing; p++) {
            pairing[p][0] = available[0][p];
            pairing[p][1] = available[1][p];
        }
        return numPairing;
    }
}
