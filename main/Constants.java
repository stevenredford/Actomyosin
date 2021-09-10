package main;

import java.io.*;
import parameters.*;

// this class contains all the physical constants used in the simulation
public class Constants {
	
	/********* PHYSICAL CONSTANTS  ***********/
	
	/** Avogadro's number = #molecules/Mole */
	static double avagadro = 6.022*10e23;
	
	/** Number of molecules in 1 nm3 for a 1µM solution. */
	static double microMolePerSquareNM = avagadro*1e-30;
	
	/** Viscosity of water in (N*sec)/m^2). */
	public static double viscosityOfWater = 0.001;
	public static double xWater = 100;
	
	/** Conversion of viscosity to standard units (pN * sec) / nm^2) */
	static public double viscosity = (xWater*viscosityOfWater * Math.pow(10,12)) / (Math.pow(10,9)*Math.pow(10,9));
	
	/** Unit of thermal energy (pN * nm) */
	static public double kT = 4.13;
	
	static String className = new String("main.Constants");
	
	static public void setXWater(double x) {
		xWater = x;
		viscosity = (xWater*viscosityOfWater * Math.pow(10,12)) / (Math.pow(10,9)*Math.pow(10,9));
	}
		
	static void staticInit() {
		Parameters.addParameter(className, "ViscosityMultiplier",100,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					setXWater(p.getDoubleValue());
				}
			}
		);
	}
}


