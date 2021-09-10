package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class NumMyoBindingSitesEvaluator extends Evaluator {
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
		numIntervalsMyosinUnbound=new int [Myosin.myosinCt];
		numMyoBindingEvents=new int [Myosin.myosinCt];
		numMyoUnbindingEvents=new int [Myosin.myosinCt];
		timeUnboundTracker=new ValueTracker(4000);
		myoHeadsBoundTracker=new ValueTracker(4000);
		numPossibleBindingSitesTracker=new ValueTracker(4000);
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "numBindingSites.dat");

	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("NumMyoBindingSitesEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		trackMyosinBinding();

	}
	public String getHeaderString() {
		return new String(name + "\t" + "avg myo heads bound" + "\t" + "avg time unbound");
	}

	public String getDataString() {
		return new String(name + "\t" + myoHeadsBoundTracker.runningAverageVal() + "\t" + timeUnboundTracker.runningAverageVal());
	}
	/**
	 * track intervals of time actins are not bound by myosin
	 */


	int [] numIntervalsUnbound;
	ValueTracker timeUnboundTracker;
	ValueTracker myoHeadsBoundTracker;
	ValueTracker numPossibleBindingSitesTracker;
	int totalMyoHeadsBound;
	public void trackMyosinBinding(){
		int [] myoBindingSites=new int [Actin.actinCt];
		totalMyoHeadsBound=0;
		for (int i=0; i<Actin.actinCt; i++){
			for (Monomer m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					myoBindingSites[i]+=1;
				}
			}
			//myoBindingSites[i]+=Actin.theActins[i].numBindingSites;
			//dataPW.print(myoBindingSites[i]+ "\t");
		}
		dataPW.println();
		
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

