package gov.pppl.spec;

import de.labathome.FortranNamelist;
import gov.pppl.spec.SpecInput.InputPhysics;

/**
 * Toolbox for the SPEC MRxMHD equilibrium code
 * @author Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
 */
public class jSPEC {

	/**
	 * read the input data from a namelist
	 * @param input_nmlists
	 * @return
	 */
	public static SpecInput parseInputNamelists(String inputNamelists) {
		SpecInput specInput = new SpecInput();
		
		// parse physicslist
		FortranNamelist inputParser = new FortranNamelist(inputNamelists, "physicslist", specInput.physics);
		specInput.physics = (InputPhysics)inputParser.getParsed();
		
        
		
		
		
		return specInput;
	}
	
	
}
