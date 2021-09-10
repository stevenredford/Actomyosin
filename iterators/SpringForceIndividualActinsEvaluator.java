package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class SpringForceIndividualActinsEvaluator extends Evaluator {

	PrintWriter dataPW2;
	String dataFileName2;
	double binSize;
	double estimatedMaxForce;
	int bin;
	int [] forceHistogram;
	int binCt;
	int maxBin=0;
	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}
		

	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		binCt= (int) (estimatedMaxForce/binSize);
		forceHistogram = new int [binCt];

	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equals("binSize"))  {
			binSize= in.nextDouble();
		}
		else if(tag.equals("estimatedMaxForce"))  {
			estimatedMaxForce= in.nextDouble();
		}

		else super.loadParameter(tag, in);
	}
	
	public boolean decide(double tm) {
		
		if(tm > stopTime) {
			if(dataPW != null) dataPW.close();
			if(dataPW2 != null) dataPW2.close();
			return true;
			
		}
		return false;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			dataPW2 = new PrintWriter(new FileWriter(new File (dataFileName2)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("SpringForceIndividualActinsEvaluator.mkDataFile(): error creating File"); }
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "filamentSpringForceData.dat");
		dataFileName2 = new String(dataPath + File.separator + "springForceHistogram.dat");
	}

	public void evaluate(double tm)  {
		trackSpringForces();
		if(tm>(stopTime-interval)){
		printHistogram();
		}
	}
	
	/**
	 * track spring forces of actin filaments
	 */
	FilamentAttachment pEnd, bEnd;
	double springForce;
	public void trackSpringForces(){
		for (int i=0; i<Actin.actinCt; i++){
			pEnd = Actin.theActins[i].getPEndAttachment();
			bEnd = Actin.theActins[i].getBEndAttachment();
			springForce = Math.abs(pEnd.getPinForce().x)+ Math.abs(bEnd.getPinForce().x);
			bin = (int) (springForce/binSize);
			forceHistogram[bin]++;
			if (bin>maxBin){
				maxBin=bin;
			}
			dataPW.print(springForce+"\t");
		}
		dataPW.println();
	}
	
	public void printHistogram () {
		for (int i=0; i<=maxBin; i++){
			dataPW2.print(forceHistogram[i] + "\t");
		}
		dataPW2.println();
	}
	
	/** Override this method to write out a datafile containing Evaluator subtype-specific
	 * data in whatever format you choose.  The method should return true if the dataFie was successfuly written,
	 * false otherwise.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}


}

