package gov.pppl.spec;

import java.util.Arrays;

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
	 */
	@namelist_variable
	public double phiedge;
			
	public SpecInput() {
		Igeometry = 3;
		Nvol = 1;
		Lrad = new int[MNvol];
		Arrays.fill(Lrad, 4);
		phiedge = 1.0;
	}
}
