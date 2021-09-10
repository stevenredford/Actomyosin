package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;
import collision.*;

public class SpringForceEvaluatorWithPulse extends Evaluator {

	double minX;
	double maxX;

	//Cumulative sum of spring forces for pointed ends on left of bundle
	double avgLeftPEndSum=0;

	//Cumulative sum of spring forces for pointed ends on right of bundle
	double avgRightPEndSum=0;

	double printForceInterval;

	double printForceCounter;

	//Cumulative sum of spring forces for barbed ends on left of bundle
	double avgLeftBEndSum=0;

	//Cumulative sum of spring forces for barbed ends on left of bundle
	double avgRightBEndSum=0;

	//Sum of magnitudes of all spring forces
	double avgMagnitudeSum=0;

	double startAveragingCounter;

	//Average of the squares of the magnitude sums for RMSD calculation
	double avgMagnitudeSumSquared=0;

	//root mean squared deviation of the spring force magnitudes
	double rMSD=0;

	//time at which averaging of force magnitudes begins--can find force average at final steady state
	double timeToStartAveragingForces;

	//indicates successful switching into state that generates 70% Fmax if 1; otherwise 0
	int switched = 0;

	int nHeadsEquil = 0;
	double equilStallF = 0;

	boolean stopAt50PFmax=false;
	boolean startAveragingForces=false;

	int averageCt;

	int evaluationCt;

	double maxPEndMagnitude=0;
	
	double magnitudeSum;


	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		minX = 0;
		maxX = Sim2D.xDimension;
		avgLeftPEndSum=0;
		avgRightPEndSum=0;
		avgLeftBEndSum=0;
		avgRightBEndSum=0;
		avgMagnitudeSum=0;
		avgMagnitudeSumSquared=0;
		averageCt = 0;
		maxPEndMagnitude=0;
		startAveragingForces=false;
		startAveragingCounter=0.0;
		rMSD=0;
		magnitudeSum=0;
		evaluationCt=0;
		printForceCounter=0;
		switched=0;
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

	public void writeDataFileHeader() {
		dataPW.println("numAHeads\t" + "numBHeads\t" + "leftPEndSum\t"+ "leftBEndSum\t"+ "rightPEndSum\t"+ "rightBEndSum\t"+ "magnitudeSums\t");
	}

	public boolean decide(double tm) {
		double maxF = MyosinMiniFilament.predictStallForce();
		if(stopAt50PFmax){
			if(evaluationCt > 0 && magnitudeSum > 0.5*maxF) {
				switched=1;
				return true;
			}
		}

		return super.decide(tm);
	}

	public void evaluate(double tm)  {

		if (tm < Actin.pulseTime){
			nHeadsEquil = MyosinMiniFilament.predictNumHeadsAtStall();
			equilStallF = nHeadsEquil*Myosin.getCrossbridgeForce();
			Actin.resistiveLoad = equilStallF*Actin.fractionStallForceApplied;


		}
		else{
			Actin.resistiveLoad=0;
		}

		trackPinForcesAtEnds();
		startAveragingCounter+=interval;
		if (startAveragingCounter >= timeToStartAveragingForces){
			startAveragingForces=true;
		}
	}



	/**
	 * calculate sum of pin forces at each end of bundle
	 */
	public void trackPinForcesAtEnds(){
		double leftPEndSum=0;									//Sum of spring forces for pointed ends on left of bundle
		double rightPEndSum=0;									//Sum of spring forces for pointed ends on right of bundle
		double leftBEndSum=0;									//Sum of spring forces for barbed ends on left of bundle
		double rightBEndSum=0;									//Sum of spring forces for barbed ends on left of bundle
		double pEndMagnitude=0;
		FilamentAttachment pEnd, bEnd;
		Actin a;

		//Sum of all spring forces
		for (int i=0; i<Actin.actinCt; i++){
			a = Actin.theActins[i];
			if(a.uVect.x > 0) {  // barbed end points to the right, i.e. towards positive X
				if (a.isPPinned()){
					pEnd = a.getPEndAttachment();
					rightPEndSum+=pEnd.getPinForce().x;
				}
				if (a.isBPinned()){
					bEnd = a.getBEndAttachment();
					rightBEndSum+=bEnd.getPinForce().x;
				}
			}
			else {
				if (a.isPPinned()){
					pEnd = a.getPEndAttachment();
					leftPEndSum+=pEnd.getPinForce().x;
				}
				if (a.isBPinned()){
					bEnd = a.getBEndAttachment();
					leftBEndSum+=bEnd.getPinForce().x;
				}
			}

		}

		double numAHeads = MyosinMiniFilament.theMiniFilaments[0].getNumAHeadsBound();
		double numBHeads = MyosinMiniFilament.theMiniFilaments[0].getNumBHeadsBound();

		magnitudeSum=(Math.abs(leftPEndSum)+ Math.abs(rightPEndSum)+ Math.abs(leftBEndSum)+ Math.abs(rightBEndSum))/Actin.actinCt;
		printForceCounter+=interval;
		if (printForceCounter>=printForceInterval){
			dataPW.println(numAHeads + "\t" + numBHeads + "\t" + leftPEndSum+"\t"+leftBEndSum+"\t"+rightPEndSum+"\t"+rightBEndSum+"\t"+ magnitudeSum);
			printForceCounter=0;
		}
		pEndMagnitude=leftPEndSum-rightPEndSum;
		if (pEndMagnitude>maxPEndMagnitude){
			maxPEndMagnitude=pEndMagnitude;
		}
		if (startAveragingForces){
			avgLeftPEndSum += leftPEndSum;
			avgRightPEndSum += rightPEndSum;
			avgLeftBEndSum += leftBEndSum;
			avgRightBEndSum += rightBEndSum;
			avgMagnitudeSum += magnitudeSum;
			avgMagnitudeSumSquared += Math.pow(magnitudeSum, 2);
			averageCt++;
		}
		evaluationCt++;
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
		if(tag.equalsIgnoreCase("startAverageTime"))  {
			timeToStartAveragingForces= in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("printForceInterval"))  {
			printForceInterval= in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("stopAt50PFmax"))  {
			stopAt50PFmax= in.nextBoolean();
		}

		else super.loadParameter(tag, in);
	}


	public String getHeaderString() {
		return new String(name + "\t" + "maxPEndMagnitude" + "\t" + "avgMagnitudeSum"  + "\t" + "predictedMaxF" + "\t" + "RMSD magnitude sum" + "\t"+"switched?"+"\t");
	}

	public String getDataString() {

		double maxF =  MyosinMiniFilament.predictStallForce();
		avgLeftPEndSum/=averageCt;
		avgRightPEndSum/=averageCt;
		avgLeftBEndSum/=averageCt;
		avgRightBEndSum/=averageCt;
		avgMagnitudeSum/=averageCt;
		avgMagnitudeSumSquared/=averageCt;
		rMSD = Math.sqrt(avgMagnitudeSumSquared - Math.pow(avgMagnitudeSum, 2));
	
	

		return new String(name + "\t" + maxPEndMagnitude + "\t" + avgMagnitudeSum + "\t" + maxF + "\t" + rMSD + "\t"+ switched + "\t");
	}

}

