package gov.pppl.spec;

import java.util.Arrays;
import java.util.Locale;

import de.labathome.namelist_variable;

/**
 * Input data for the SPEC MrxMHD equilibrium code
 * @author Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
 */
public class SpecInput {
	
	/**
	 * maximum resolution, i.e. array dimensions
	 */
	public final int MNvol = 256; // maximum value of \c Nvol
	public final int MMpol =  32; // maximum value of \c Mpol
	public final int MNtor =  16; // maximum value of \c Ntor
	
	/**
	 * input namelist 'physicslist'
	 */
	public class InputPhysics {
		/**
		 * selects Cartesian, cylindrical or toroidal geometry;
		 * <ul>
		 * <li> \c Igeometry=1 : Cartesian; geometry determined by \f$R\f$; </li>
		 * <li> \c Igeometry=2 : cylindrical; geometry determined by \f$R\f$; </li>
		 * <li> \c Igeometry=3 : toroidal; geometry determined by \f$R\f$ *and* \f$Z\f$; </li>
		 * </ul>
		 */
		@namelist_variable
		public int Igeometry;
		
		/**
		 * number of volumes
		 * <ul>
		 * <li> each volume \f${\cal V}_l\f$ is bounded by the \f${\cal I}_{l-1}\f$ and \f${\cal I}_{l}\f$ interfaces </li>
		 * <li> note that in cylindrical or toroidal geometry, \f${\cal I}_{0}\f$ is the degenerate coordinate axis </li>
		 * <li> constraint: \c Nvol<=MNvol </li>
		 * </ul>
		 */
		@namelist_variable
		public int Nvol;
		
		/**
		 * Chebyshev resolution in each volume
		 * <ul>
		 * <li> constraint : \c Lrad(1:Mvol) >= 2 </li>
		 * </ul>
		 */
		@namelist_variable
		public int[] Lrad;
		
		/**
		 * total enclosed toroidal magnetic flux
		 * unit: Vs
		 */
		@namelist_variable
		public double phiedge;
		
		/**
		 * Instantiate the SPEC input data object and initialize the contents to their default values
		 * as in the original Fortran source code.
		 */
		public InputPhysics() {
			Igeometry = 3;
			Nvol = 1;
			Lrad = new int[MNvol];
			Arrays.fill(Lrad, 4);
			phiedge = 1.0;
		}
		
		/**
		 * Print the contents of this class to the screen.
		 */
		public void dump() {
			System.out.println(String.format(Locale.ENGLISH, "Igeometry = %d", Igeometry));
			System.out.println(String.format(Locale.ENGLISH, "  phiedge = %g", phiedge));
			System.out.println(String.format(Locale.ENGLISH, "     Nvol = %d", Nvol));
			for (int i=0; i<Nvol; ++i) { // TODO: Nvol+Lfreeboundary
				System.out.println(String.format(Locale.ENGLISH, "  Lrad[%d] = %d", i, Lrad[i]));
			}
		}
	}
	
	/**
	 * input namelist 'physicslist'
	 */
	public InputPhysics physics;
	
	public SpecInput() {
		physics = new InputPhysics();
	}
	
	/**
	 * Print the contents of this class to the screen.
	 */
	public void dump() {
		if (physics != null) physics.dump();
	}
}
