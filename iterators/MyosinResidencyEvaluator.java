package iterators;


import java.util.*;
import analysis.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyosinResidencyEvaluator extends Evaluator {

	static int evaluationCt;
	
	static double [] aveMyos;
	static int maxVTracks = 2;
	static int longAveCt = 400;
	static int reallyLongAveCt = 2000;
	ValueTracker [] longAveMyos = new ValueTracker[maxVTracks];
	ValueTracker [] reallyLongAveMyos = new ValueTracker[maxVTracks];
	static double timeToStartDiffTrack = 10.0; //seconds
	static double maxLongAveDiff = 0;
	static double maxReallyLongAveDiff = 0;

	public void init(String path, String name) throws Exception {
		super.init(path,name);
		for (int i=0;i<maxVTracks;i++) {
			longAveMyos[i] = new ValueTracker(longAveCt);
			reallyLongAveMyos[i] = new ValueTracker(reallyLongAveCt);
		}
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
		dataFileName = new String(dataPath + File.separator + "myoResidencyData.dat");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("MyosinResidencyEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		getAverageMyosinsAttached();
	}



	/**
	 * get average myosins attached to filament over last report interval
	 */

	public void getAverageMyosinsAttached() {
		aveMyos = new double[Actin.actinCt];
		for (int i=0; i<Actin.actinCt;i++) {
			aveMyos[i] = Actin.theActins[i].aveMyosAttached();
			longAveMyos[i].registerValue(aveMyos[i]);
			reallyLongAveMyos[i].registerValue(aveMyos[i]);
		}
		dataPW.println(getDataString());
	}

	/** Override this method to write out a datafile containing Evaluator subtype-specific
	 * data in whatever format you choose.  The method should return true if the dataFie was successfully written,
	 * false otherwise.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equals("startAveragingForces"))  {
			
		}

		else super.loadParameter(tag, in);
	}


	public void writeDataFileHeader () {
		dataPW.println(getHeaderString());
	}
	
	public String getHeaderString() {
		String headStr = "time" + "\t";
		for (int i=0;i<Actin.actinCt;i++) {
			headStr += "actin" + i + "\t";
		}
		for (int i=0;i<Actin.actinCt;i++) {
			headStr += "actin" + i + "LongAve" + "\t";
		}
		for (int i=0;i<Actin.actinCt;i++) {
			headStr += "actin" + i + "ReallyLongAve" + "\t";
		}
		headStr += "maxLongAveDiff" + "\t";
		headStr += "maxReallyLongAveDiff" + "\t";

		return headStr;
	}

	public String getDataString() {
		String datStr = Sim2D.simulationTime + "\t";
		for (int i=0;i<Actin.actinCt;i++) {
			datStr += aveMyos[i] + "\t";
		}
		for (int i=0;i<Actin.actinCt;i++) {
			datStr += longAveMyos[i].averageVal() + "\t";
		}
		for (int i=0;i<Actin.actinCt;i++) {
			datStr += reallyLongAveMyos[i].averageVal() + "\t";
		}
		
		if (Sim2D.simulationTime > timeToStartDiffTrack) {
			double longAveDiff = Math.abs(longAveMyos[0].averageVal()-longAveMyos[1].averageVal());
			double reallyLongAveDiff = Math.abs(reallyLongAveMyos[0].averageVal()-reallyLongAveMyos[1].averageVal());
			if (longAveDiff > maxLongAveDiff) { maxLongAveDiff = longAveDiff; }
			if (reallyLongAveDiff > maxReallyLongAveDiff) { maxReallyLongAveDiff = reallyLongAveDiff; }
		}
		
		datStr += maxLongAveDiff + "\t";
		datStr += maxReallyLongAveDiff + "\t";
		
		return datStr;

	}

}

