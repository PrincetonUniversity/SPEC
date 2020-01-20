package gov.pppl.spec;

import java.nio.file.Files;
import java.nio.file.Paths;

public class TestSpecReader {

	public static void main(String[] args) {
		System.out.println("readInput() {");
		readInput();
		System.out.println("}\n");
		System.out.println("readOutput() {");
		readOutput();
		System.out.println("}\n");
	}
	
	public static void readInput() {
		String filename = "/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/Utilities/proto/demo_nml.sp";
		String inputFile = "";
		try {
			inputFile = new String(Files.readAllBytes(Paths.get(filename)));
		} catch (Exception e) {
			e.printStackTrace();
		}
		SpecInput input = jSPEC.parseInputNamelists(inputFile);
		input.dump();
	}
	
	public static void readOutput() {
		SpecOutput s = new SpecOutput("/home/jonathan/Uni/04_PhD/00_programs/SPEC/SPEC/Utilities/proto/demo_nml.sp.h5");
		s.dump();		
	}

}
