package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class ActinMyosinBindingEvaluator extends Evaluator {
	int [] timeUnboundIntervalCt;
	double startTime=0;
	double startTimeInterval;
	int currentStartTimeInterval;
	int numStartTimeIntervals;


	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		numIntervalsUnbound=new int [Actin.actinCt];
		tempMyoHeadsBound=new int [Actin.actinCt];
		numIntervalsMyosinUnbound=new int [Myosin.myosinCt];
		numMyoBindingEvents=new int [Myosin.myosinCt];
		numMyoUnbindingEvents=new int [Myosin.myosinCt];
		myoHeadsBound=new int [Actin.actinCt];
		timeUnboundTracker=new ValueTracker((int) (super.stopTime/interval));
		myoHeadsBoundTracker=new ValueTracker((int) (super.stopTime/interval));
		myoHeadsBoundIndividualActinsTracker= new ValueTracker((int) (super.stopTime/interval));
		numPossibleBindingSitesTracker=new ValueTracker((int) (super.stopTime/interval));
		maxHeadsOnSingleFilament = 0;
		avgActinVelocityWhenBound =0;
		actinVelocityTrackCount=0;
		aveDutyRatio=0;
		count=0;
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "actinMyosinBindingData.dat");

	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();

		} catch (IOException ioe) { System.out.println("ActinMyosinBindingEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		trackMyosinBinding();
		trackNumAvailableBindingSites();
		averageDutyRatio();

	}
	public String getHeaderString() {
		return new String(name + "\t" + "avg myo heads bound" + "\t" + "avg time unbound" + "\t" + "periods unbound/fil/s" + "\t" + " avg myos bound to single actin" + "\t" + "max myos bound to single actin" +"\t"+ "ave duty ratio"+"\t");
	}

	public String getDataString() {
		numUnboundPeriodsPerFilamentPerSecond /= Actin.actinCt;
		numUnboundPeriodsPerFilamentPerSecond /= Sim2D.simulationTime;
		aveDutyRatio /= count;
		return new String(name + "\t" + myoHeadsBoundTracker.runningAverageVal() + "\t" + timeUnboundTracker.runningAverageVal() + "\t" + numUnboundPeriodsPerFilamentPerSecond + "\t" + myoHeadsBoundIndividualActinsTracker.averageVal()+"\t" + maxHeadsOnSingleFilament + "\t" + aveDutyRatio + "\t");
	}
	/**
	 * track intervals of time actins are not bound by myosin
	 */


	int [] numIntervalsUnbound;
	ValueTracker timeUnboundTracker;
	ValueTracker myoHeadsBoundTracker;
	ValueTracker myoHeadsBoundIndividualActinsTracker;
	ValueTracker numPossibleBindingSitesTracker;
	static public int [] myoHeadsBound;
	int totalMyoHeadsBound;
	int [] tempMyoHeadsBound;
	double numUnboundPeriodsPerFilamentPerSecond;
	double maxHeadsOnSingleFilament;
	static public double avgActinVelocityWhenBound;
	static public int actinVelocityTrackCount;
	
	
	public void trackMyosinBinding(){
		myoHeadsBound=new int [Actin.actinCt];
		totalMyoHeadsBound=0;
		for (int i=0; i<Actin.actinCt; i++){
			for (Monomer m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					myoHeadsBound[i]+=1;
				}
			}
			totalMyoHeadsBound+=myoHeadsBound[i];
			dataPW.print(myoHeadsBound[i]+"\t");
			myoHeadsBoundIndividualActinsTracker.registerValue(myoHeadsBound[i]);
			if (myoHeadsBound[i]>maxHeadsOnSingleFilament){
				maxHeadsOnSingleFilament = myoHeadsBound[i];
			}
			if (myoHeadsBound[i]==0){
				numIntervalsUnbound[i]+=1;
				if (tempMyoHeadsBound[i] != 0){
					numUnboundPeriodsPerFilamentPerSecond +=1;
				}
			}
			else if (myoHeadsBound[i]!=0){
				if (numIntervalsUnbound[i]!=0){
					timeUnboundTracker.registerValue(numIntervalsUnbound[i]*interval);
					numIntervalsUnbound[i]=0;
				}
//				if (i==0 && ActinVelocityEvaluator.startCounting){
//					avgActinVelocityWhenBound=+ActinVelocityEvaluator.velocity[i];
//					actinVelocityTrackCount++;
//				}
			}
			tempMyoHeadsBound[i] = myoHeadsBound[i];

		}
		dataPW.println();
		myoHeadsBoundTracker.registerValue(totalMyoHeadsBound);
	}
	public void trackNumAvailableBindingSites(){
		for (int i=0; i<Myosin.myosinCt; i++){
			numPossibleBindingSitesTracker.registerValue(Myosin.theMyosins[i].numPossibles);
			//System.out.print(Myosin.theMyosins[i].numPossibles+"\t");
		}
		//System.out.println();
	}
	double aveDutyRatio=0;
	int count =0;
	public void averageDutyRatio(){
		aveDutyRatio+=Myosin.getAveDutyRatio();
		count+=1;

	}
	int [] numIntervalsMyosinUnbound;
	int [] numIntervalsMyosinBound;
	int [] numMyoBindingEvents;
	int [] numMyoUnbindingEvents;
	public void trackMyosinBindingRate(){
		for (int i=0; i<Myosin.myosinCt; i++){
			if (Myosin.theMyosins[i].boundMon == null){
				numIntervalsMyosinUnbound[i]+=1;
				if (numIntervalsMyosinBound[i]!=0){
					numMyoUnbindingEvents[i]+=1;
					numIntervalsMyosinBound[i]=0;
				}
			}
			else{
				if (numIntervalsMyosinUnbound[i]!=0){
					numMyoBindingEvents[i]+=1;
					numIntervalsMyosinUnbound[i]=0;
				}
				numIntervalsMyosinBound[i]+=1;
			}
		}
	}
	/*public void boundMyoTangStretch (boolean startPowerStroke){//put this in the evaluator, feed in boundMyoPowerStroke as argument when called at every evaluation
		if (startPowerStroke){
			prefix=1;
			boundMyoPowerStroke=false;
		}
		tangdist=boundMyo.getTangDist();
	}
	 */
	public void loadParameter(String tag, AMInputStream in)  throws Exception {

		if(tag.equalsIgnoreCase("NumStartTimeIntervals"))  {
			numStartTimeIntervals = in.nextInt();
		}

		else{
			super.loadParameter(tag, in);
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


}

