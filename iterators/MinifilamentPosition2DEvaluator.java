package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MinifilamentPosition2DEvaluator extends Evaluator {


	int timeIntervalNumber=0;
	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}
		

	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		timeIntervalNumber=0;
	
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "miniFilamentPositionData.dat");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
		
			
		} catch (IOException ioe) { System.out.println("MinifilamentPositionEvaluator.mkDataFile(): error creating File"); }
	}

	public void evaluate(double tm)  {
		trackVelocities();
	}
	

	/**
	 * track positions of all minifilaments
	 */
	
	public void trackVelocities(){
		timeIntervalNumber+=1;
		
			
		dataPW.println(1 +"\t"+timeIntervalNumber+"\t"+MyosinMiniFilament.theMiniFilaments[0].cm.x+"\t"+MyosinMiniFilament.theMiniFilaments[0].cm.y);
			
		
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

