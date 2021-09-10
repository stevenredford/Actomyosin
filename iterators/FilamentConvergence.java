package iterators;


import java.util.*;
import java.io.*;

import main.*;
import util.*;
import io.*;

public class FilamentConvergence extends Evaluator implements ActinListener  {
	
	
	static int checkConvCntr = 0;
	static double convergenceCheckInterval = 1;
	static double nextConvergenceCheck = 1;
	boolean converged = false;
	double [] aveFilLength;
	double [] stdvFilLength;
	double [] aveLifetime;
	double [] stdvLifetime;
	double steadyStateThreshold;
	
	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus lifetimeHistogram;
	
	/** A histogram that keeps running track of distribution of filament lengths.. */
	HistogramPlus lengthHistogram;
	
	public FilamentConvergence() {}

	public void init(String path, String name) throws Exception  {
		Actin.setMonomerListener(this);
		aveFilLength = new double [1000];
		stdvFilLength = new double [1000];
		aveLifetime = new double [1000];
		stdvLifetime = new double[1000];
		nextConvergenceCheck = 0;
		checkConvCntr = 0;
		converged = false;
		super.init(path,name);
	}
	
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		lifetimeHistogram = new HistogramPlus(20,0,100,dataPath,"actinLifetimes",null,false,true,false,false);
		lengthHistogram = new HistogramPlus(50,0,0.5*Sim2D.xDimension,dataPath,"actinLengths",null,false,true,false,false);
	}

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("SteadyStateThreshold"))  {
			steadyStateThreshold = in.nextDouble();
		}
		else if(tag.equals("ConvergenceCheckInterval"))  {
			convergenceCheckInterval = in.nextDouble();
		}
		else super.loadParameter(tag,in);
	}

	public void monomerReleased(Monomer m) {
		double lifetime = Sim2D.simulationTime - m.birthTime;
		if(m.birthTime > 0) {
			lifetimeHistogram.addValue(lifetime);
		}
	}
		
	void updateLengthHistogram() {
		lengthHistogram.clearBins();
		for (int i=0;i<Actin.actinCt;i++) {
			lengthHistogram.addValue(Actin.theActins[i].physicalLength);
		}
	}

	
	void writeHistograms() {
		lengthHistogram.writeToFile();
		lifetimeHistogram.writeToFile();
	}
	
	public boolean stop(double tm)  {
		if(super.stop(tm)) return true;
		if(tm >= nextConvergenceCheck) {
			nextConvergenceCheck += convergenceCheckInterval;
			return checkConvergence(tm);
		}
		return false;
	}
	
	public boolean gotHit()  {
		return converged;
	}

	public String getHeaderString() {
		return new String("AverageFilamentLength" + "\t" + "StdvFilamentLength" + "\t" + "AverageFilamentLifetime" + "\t" + "StdvFilamentLifetime");
	}
	
	public String getDataString() {
		return new String(aveFilLength[checkConvCntr] + "\t" + stdvFilLength[checkConvCntr] + "\t" + aveLifetime[checkConvCntr] + "\t" + stdvLifetime[checkConvCntr]);
	}
	
	/** Override this method to specify whether this Evaluator writes data files.  These
	 * files MUST be written to the path defined by the variable dataPata to be compiled
	 * by the iterator.
	 *
	 * */
	public boolean hasData() {
		return false;
	}
	
	public String writeDataFile() throws Exception {
		File f = new File(dataPath + File.separator + "histograms.txt");
		AMOutputStream out = new AMOutputStream(f);
		lengthHistogram.writeHistogram(out);
		lifetimeHistogram.writeHistogram(out);
		return f.getAbsolutePath();
	}
		
	private boolean checkConvergence(double tm) {
		storeData(checkConvCntr);
		//writeHistograms();
		if(checkConvCntr > 0) {
			boolean [] singleConvergence = new boolean[4];
			singleConvergence[0] = checkSingleConvergence(steadyStateThreshold,aveFilLength[checkConvCntr-1],aveFilLength[checkConvCntr]);
			singleConvergence[1] = checkSingleConvergence(steadyStateThreshold,stdvFilLength[checkConvCntr-1],stdvFilLength[checkConvCntr]);
			singleConvergence[2] = checkSingleConvergence(steadyStateThreshold,aveLifetime[checkConvCntr-1],aveLifetime[checkConvCntr]);
			singleConvergence[3] = checkSingleConvergence(steadyStateThreshold,stdvLifetime[checkConvCntr-1],stdvLifetime[checkConvCntr]);
			if(checkAllConvergence(singleConvergence[0],singleConvergence[1],singleConvergence[2],singleConvergence[3])) {
				converged = true;
				return true;
			}
		}
		checkConvCntr++;
		return false;
	}
			
	public static boolean checkSingleConvergence(double stablizedSlope, double valueA, double valueB) {
		return (Math.abs(valueA-valueB)/Math.abs(valueA))<=stablizedSlope;
	}
	
	static public boolean checkAllConvergence(boolean aveFilLengthConv, boolean stdvFilLengthConv, boolean aveLifetimeConv, boolean stdvLifetimeConv){
		if(aveFilLengthConv && stdvFilLengthConv && aveLifetimeConv && stdvLifetimeConv){
			return true;
		}
		return false;
	}
	

	private void storeData(int i){
		aveFilLength[i] = lengthHistogram.getMean();
		stdvFilLength[i] = lengthHistogram.getStdev();
		aveLifetime[i] = lifetimeHistogram.getMean();
		stdvLifetime[i] = lifetimeHistogram.getStdev();
	}

	
}
