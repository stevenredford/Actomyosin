package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoSpringLengthEvaluator extends Evaluator {
	

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

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("MyoSpringLengthEvaluator.mkDataFile(): error creating File"); }
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "myoSpringLengthData.dat");

	}


	public void evaluate(double tm)  {
		trackMyoSpringLengths();


	}

	/**
	 * track intervals of time actins are not bound by myosin
	 */



	
	public void trackMyoSpringLengths(){
		
		for (int i=0; i<Myosin.myosinCt; i++){
			double myoTangDist =0;
			
			if (Myosin.theMyosins[i].boundMon != null){
				myoTangDist=Myosin.theMyosins[i].getTangDist();
			}
			dataPW.print(myoTangDist+"\t");
		}
		dataPW.println();
	}


	public void loadParameter(String tag, AMInputStream in)  throws Exception {

			super.loadParameter(tag, in);
		
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

