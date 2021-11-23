package population.availability;

import availability.AbstractAvailability;
import person.AbstractIndividualInterface;
import random.RandomGenerator;

/**
 *
 * @author Ben Hui 
 */
public class MSMAvailablity extends AbstractAvailability {
    
    /**
	 * 
	 */
	private static final long serialVersionUID = 475269122305141031L;
	protected AbstractIndividualInterface[][] available;    
    protected AbstractIndividualInterface[][] pairing = null;             
    
    public MSMAvailablity(RandomGenerator RNG) {
        super(RNG);        
    }        

    @Override
    public int generatePairing() {        
        int numAvailable  = available[0].length; // MSM only        
        int numPairing = numAvailable/2;   
        
        util.ArrayUtilsRandomGenerator.shuffleArray(available[0], getRNG());        
        pairing = new AbstractIndividualInterface[numPairing][2];
        
        int r = 0, pt = 0;        
        while(r < numPairing){
            pairing[r][0] = available[0][pt];
            pt++;
            pairing[r][1] = available[0][pt];
            pt++;                        
            r++;
        }                                       
        return numPairing;       
    }    
  

    @Override
    public void setAvailablePopulation(AbstractIndividualInterface[][] available) {        
       this.available = available;
    }
    
    @Override
    public AbstractIndividualInterface[][] getPairing() {
        return pairing;
    }

    @Override
    public boolean setParameter(String id, Object value) {
        throw new UnsupportedOperationException("Not supported in this verison."); 
    }

    @Override
    public Object getParameter(String id) {
        throw new UnsupportedOperationException("Not supported in this verison."); 
    }
    
    @Override
    public boolean removeMemberAvailability(AbstractIndividualInterface p) {
        throw new UnsupportedOperationException("Not supported in this verison."); 
    }

    @Override
    public boolean memberAvailable(AbstractIndividualInterface p) {
       throw new UnsupportedOperationException("Not supported in this verison."); 
    }
    
}
