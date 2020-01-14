package gov.pppl.spec;

import de.labathome.FortranNamelist;

/**
 * Toolbox for the SPEC MRxMHD equilibrium code
 * @author Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
 */
public class jSPEC {

	public static SpecInput read_input(String input_nmlists) {
		SpecInput physicslist = new SpecInput();
		
		FortranNamelist input_parser = new FortranNamelist(input_nmlists, "physicslist", physicslist);
		physicslist = (SpecInput)input_parser.getParsed();
		
		return physicslist;
	}
	
	
}
