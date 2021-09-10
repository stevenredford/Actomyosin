package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class ActinMyosinBindingSingleFilamentEvaluator extends Evaluator {
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
		timeUnboundTracker=new ValueTracker(4000);
	
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("ActinMyosinBindingSingleFilamentEvaluator.mkDataFile(): error creating File"); }
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "actinMyosinBindingData.dat");

	}


	public void evaluate(double tm)  {
		//trackMyosinBinding();


	}

	/**
	 * track intervals of time actins are not bound by myosin
	 */


	int [] numIntervalsUnbound;
	ValueTracker timeUnboundTracker;
	boolean bound;
	/*public void trackMyosinBinding(){
		
		for (int i=0; i<Actin.actinCt; i++){
			bound=false;
			double myoTangDist=0;
			for (Monomer m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (m.isFreeMyosin()){
					numIntervalsMyoUnbound+=1;
					if (numIntervalsMyoBound!=0){
						numMyoUnbindingEvents+=1;
						numIntervalsMyoBound=0;
					}
					
				}
				else{
					numIntervalsMyoBound+=1;
					if(numIntervalsMyoUnbound!=0){
						numMyoBindingEvents+=1;
						numIntervalsMyoUnbound=0;
					}
				}
			}
			dataPW.print(myoTangDist +"\t" + Actin.theActins[i].cm.x+ "\t");
			if (!bound){
				numIntervalsUnbound[i]+=1;
			}
			else if (bound && numIntervalsUnbound[i]!=0){
				timeUnboundTracker.registerValue(numIntervalsUnbound[i]*interval);
				numIntervalsUnbound[i]=0;
				//also track amount actin filament moves during unbound time
			}
			
		}
		dataPW.println();
	}*/


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

