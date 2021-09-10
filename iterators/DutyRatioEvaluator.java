package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;

public class DutyRatioEvaluator extends Evaluator {
	static String actinMyoSepString = " ";
	static String sepString = ";";
	boolean firstEval = true;
	static String dutyRatioAssayFileName;
	static boolean dutyRatioFileMade = false;
	static File dutyRatioFile;
	static FileWriter dRFW;
	static PrintWriter dRPW;
	double nameOfInput;
	static double curInterval;
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
		curInterval = this.interval;
	}
	
	public void reset()  {
		super.reset();
		firstEval = true;
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dutyRatioAssayFileName = dataPath + File.separator + "dutyratio.dat";
	}

	
	/** Override this method to specify whether this Evaluator writes data files.  These
	 * files MUST be written to the path defined by the variable dataPata to be compiled
	 * by the iterator.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}

	public boolean stop(double tm)  {
		if(tm > nextEval) {
			if (firstEval){
				firstEval = false;
			} else {
				writeActinData();
				nextEval += interval;
			}
		}
			
		if(super.stop(tm)) {
			// end of evaluation actions here, such as writing assembled statistics/histograms to file
			return true;
		}
		
		return false;
	}
	
	public static void mkDutyRatioDataFile () {
		try {
			dutyRatioFile = new File (dutyRatioAssayFileName);
			dRFW = new FileWriter(dutyRatioFile);
			dRPW = new PrintWriter(dRFW,true);
			
			//first line
			dRPW.print("ThreeBeadDutyRatioOutput: adpReleaseRateMult=" + Myosin.adpReleaseRateMult);
			dRPW.print(" Visc=" + Constants.xWater + "xWater");
			dRPW.print(" ResistiveLoad=" + Actin.resistiveLoad + "pN");
			dRPW.println (" NoNameRatio=" + Actin.getNoNameRatio());
			
			//header lines
			dRPW.print("simTime ");
			dRPW.print("dutyRatio");
			dRPW.println();
			
		} catch (IOException ioe) { System.out.println("An error creating mkActinDataFile"); }
		dutyRatioFileMade = true;
	}
	
	public static void writeActinData () {
		if (!dutyRatioFileMade) mkDutyRatioDataFile();
		
		printValue(Sim2D.simulationTime,dRPW,actinMyoSepString);
		
		printValue(Myosin.getAveDutyRatio(),dRPW,actinMyoSepString);

		dRPW.println("");
	}
	
	public static double getDutyRatioInterval () {
		return curInterval;
	}
	
	
	public void cleanUp()  {}
				

}
