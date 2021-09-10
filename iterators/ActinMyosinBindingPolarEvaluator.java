package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class ActinMyosinBindingPolarEvaluator extends Evaluator {
	int [] timeUnboundIntervalCt;
	double startTime=0;
	double startTimeInterval;
	int currentStartTimeInterval;
	int numStartTimeIntervals;
	String dataFileName2 = null;
	PrintWriter dataPW2;
	String dataFileName3 = null;
	PrintWriter dataPW3;
	String dataFileName4 = null;
	PrintWriter dataPW4;
	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();

	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "headsBoundLeft.dat");
		dataFileName2 = new String(dataPath + File.separator + "headsBoundRight.dat");
		dataFileName3 = new String(dataPath + File.separator + "bindingSitesLeft.dat");
		dataFileName4 = new String(dataPath + File.separator + "bindingSitesRight.dat");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			dataPW2 = new PrintWriter(new FileWriter(new File (dataFileName2)),true);
			dataPW3 = new PrintWriter(new FileWriter(new File (dataFileName3)),true);
			dataPW4 = new PrintWriter(new FileWriter(new File (dataFileName4)),true);
			writeDataFileHeader();

		} catch (IOException ioe) { System.out.println("ActinMyosinBindingEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		trackPolarMyosinBinding();
		trackNumAvailableBindingSites();


	}



	public void trackPolarMyosinBinding(){
		dataPW.println(MyosinMiniFilament.theMiniFilaments[0].getNumAHeadsBound() + "\t");
		dataPW2.println(MyosinMiniFilament.theMiniFilaments[0].getNumBHeadsBound() + "\t");
		/*for (int i=0; i<Actin.actinCt; i++){
			int myoHeadsBound=0;
			for (Monomer m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					myoHeadsBound+=1;
				}
			}
			if (Actin.theActins[i].initPolarity==0){
				dataPW.print(myoHeadsBound+"\t");
			}
			else if (Actin.theActins[i].initPolarity==Math.PI){
				dataPW2.print(myoHeadsBound+"\t");
			}*/

		//				if (i==0 && ActinVelocityEvaluator.startCounting){
		//					avgActinVelocityWhenBound=+ActinVelocityEvaluator.velocity[i];
		//					actinVelocityTrackCount++;
		//				}

	}
	public void trackNumAvailableBindingSites(){
		int numAPossibles=0; 
		int numBPossibles=0;
		for (int i=0; i<Myosin.myosinCt; i++){
			if (Myosin.theMyosins[i].iAmAnAMyo){
				numAPossibles=numAPossibles+Myosin.theMyosins[i].numPossibles;
			}
			else{
				numBPossibles=numBPossibles+Myosin.theMyosins[i].numPossibles;
			}


		}
		dataPW3.println(numAPossibles);
		dataPW4.println(numBPossibles);
		
	}


	/*public void boundMyoTangStretch (boolean startPowerStroke){//put this in the evaluator, feed in boundMyoPowerStroke as argument when called at every evaluation
		if (startPowerStroke){
			prefix=1;
			boundMyoPowerStroke=false;
		}
		tangdist=boundMyo.getTangDist();
	}
	 */



	/** Override this method to write out a datafile containing Evaluator subtype-specific
	 * data in whatever format you choose.  The method should return true if the dataFie was successfuly written,
	 * false otherwise.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}


}

