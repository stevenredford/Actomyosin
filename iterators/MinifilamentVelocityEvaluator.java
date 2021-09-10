package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MinifilamentVelocityEvaluator extends Evaluator {

	double [] myoXVal;
	double [] oldMyoXVal;
	double velocity;
	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}
		

	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		myoXVal=new double [MyosinMiniFilament.miniFilamentCt];
		oldMyoXVal=new double [MyosinMiniFilament.miniFilamentCt];
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "miniFilamentVelocityData.dat");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("MinifilamentVelocityEvaluator.mkDataFile(): error creating File"); }
	}

	public void evaluate(double tm)  {
		trackVelocities();
	}
	

	/**
	 * track velocities of all minifilaments
	 */
	public void trackVelocities(){
		for (int i=0; i<MyosinMiniFilament.miniFilamentCt; i++){
			oldMyoXVal[i]=myoXVal[i];
			myoXVal[i]=MyosinMiniFilament.theMiniFilaments[i].cm.x;
			velocity=(myoXVal[i]-oldMyoXVal[i])/interval;
			dataPW.print(velocity+"\t");
		}
		dataPW.println();
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

