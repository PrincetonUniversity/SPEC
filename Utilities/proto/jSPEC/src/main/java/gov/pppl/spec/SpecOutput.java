package gov.pppl.spec;

import java.io.IOException;
import java.util.Locale;

import ucar.ma2.DataType;
import ucar.nc2.NetcdfFile;

/**
 * Output data for the SPEC MrxMHD equilibrium code
 * @author Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
 */
public class SpecOutput {

	/**
	 * mirror image of input data
	 */
	public SpecInput input;
	
	/** Initialize complex datatypes. */
	public SpecOutput() {
		input = new SpecInput();
	}
	
	/**
	 * Initalize complex datatypes and load SpecOutput contents from a HDF5 file identified by {@code filename}.
	 * @param filename path to the HDF5 file to load
	 */
	public SpecOutput(String filename) {
		this();
		try {
			NetcdfFile file = NetcdfFile.open(filename);
			loadFrom(file);
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Initalize complex datatypes and load SpecOutput contents from an already-open NetCDF file identified by {@code file}.
	 * @param file open file to load the data from
	 */
	public SpecOutput(NetcdfFile file) {
		this();
		try {
			loadFrom(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Load SpecOutput contents from an already-open NetCDF file identified by {@code file}.
	 * @param file open file to load the data from
	 * @return initialized SpecOutput object
	 */
	public SpecOutput loadFrom(NetcdfFile file) throws IOException {
		input.physics.Igeometry = file.findVariable("/input/physics/Igeometry").readScalarInt();
		input.physics.phiedge = file.findVariable("/input/physics/phiedge").readScalarDouble();
		input.physics.Nvol = file.findVariable("/input/physics/Nvol").readScalarInt();
		input.physics.Lrad = (int[])file.findVariable("/input/physics/Lrad").read().get1DJavaArray(DataType.INT);
		return this;
	}
	
	/**
	 * Print the contents of this class to the screen.
	 */
	public void dump() {
		input.dump();
	}
}
