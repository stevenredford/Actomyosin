package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;
import collision.*;

public class SpringForceManyActinsEvaluator extends Evaluator {

	//Sum of magnitudes of all spring forces
	double [] avgMagnitudeSum;

	//Average of the squares of the magnitude sums for RMSD calculation
	double [] avgMagnitudeSumSquared;

	//root mean squared deviation of the spring force magnitudes
	double rMSD;

	//time at which averaging of force magnitudes begins--can find force average at final steady state
	double timeToStartAveragingForces;
	
	boolean startAveragingForces=false;
	
	int evaluationCt;
	
	int numActins;
	
	double maxPEndMagnitude=0;
	
	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}
		

	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		numActins = Actin.actinCt;
		avgMagnitudeSum=new double[numActins];
		avgMagnitudeSumSquared=new double[numActins];
		evaluationCt = 0;
		maxPEndMagnitude=0;
		startAveragingForces=false;
		rMSD=0;
		j=0;
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "pinForceData.dat");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("SpringForceEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		trackPinForcesAtEnds();
		if (tm >= timeToStartAveragingForces){
			startAveragingForces=true;
		}
	}



	/**
	 * calculate sum of pin forces at each end of bundle
	 */
	int j=0;
	public void trackPinForcesAtEnds(){
		double [] pF = new double[Actin.actinCt];
		FilamentAttachment pEnd, bEnd;
		
		double numAHeads = MyosinMiniFilament.theMiniFilaments[0].getNumAHeadsBound();
		double numBHeads = MyosinMiniFilament.theMiniFilaments[0].getNumBHeadsBound();
		dataPW.print(numAHeads + "\t" + numBHeads + "\t");

		//Sum of all spring forces
		for (int i=0; i<Actin.actinCt; i++){
			if (Actin.theActins[i].isPPinned()){
				pEnd = Actin.theActins[i].getPEndAttachment();
				pF[i]+=pEnd.getPinForce().x;
			}
			if (Actin.theActins[i].isBPinned()){
				bEnd = Actin.theActins[i].getBEndAttachment();
				pF[i]+=bEnd.getPinForce().x;
			}
			dataPW.print("Actin " + i + "\t" + pF[i] + "\t");
		}
		dataPW.println();
		if (startAveragingForces){
			for (int i=0; i<Actin.actinCt; i++){
				avgMagnitudeSum[i] += pF[i];
				avgMagnitudeSumSquared[i] += Math.pow(pF[i], 2);
			}
			evaluationCt++;
		}
	}

	/** Override this method to write out a datafile containing Evaluator subtype-specific
	 * data in whatever format you choose.  The method should return true if the dataFie was successfuly written,
	 * false otherwise.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equals("StartAverageTime"))  {
			timeToStartAveragingForces= in.nextDouble();
		}
		else super.loadParameter(tag, in);
	}


	public String getHeaderString() {
		String s = new String(name + "\t");
		for(int i = 0; i < numActins; i++) {
			s = s + "Actin" + i + "\t" + "avgMagnitudeSum\t"  + "RMSD magnitude sum\t";
		}
		return s;
	}

	public String getDataString() {
		String s = new String(name + "\t");
		for(int i = 0; i < numActins; i++) {
			avgMagnitudeSum[i]/=evaluationCt;
			avgMagnitudeSumSquared[i]/=evaluationCt;
			rMSD = Math.sqrt(avgMagnitudeSumSquared[i] - Math.pow(avgMagnitudeSum[i], 2));
			s = s + "Actin" + i + "\t" + avgMagnitudeSum[i] + "\t" + rMSD + "\t";
		}
		return s;
	}
}

