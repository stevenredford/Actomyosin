package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class SpringForceActinPairEvaluator extends Evaluator {


	static public double [] maxForce;

	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		maxForce=new double[2];
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "pinForceData.dat");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("SpringForceActinPairEvaluator.mkDataFile(): error creating File"); }
	}

	public void evaluate(double tm)  {
		trackForceActinPair();
	}


	/**
	 * calculate sum of pin forces at each end of bundle
	 */
	public void trackForceActinPair(){

		for (int i=0; i<Actin.actinCt; i++){
			FilamentAttachment pEnd=Actin.theActins[i].getPEndAttachment();
			double force=Math.abs(pEnd.getPinForce().x);
			if(force>maxForce[i]){
				maxForce[i]=force;
			}
			dataPW.print(force + "\t");
		}
		dataPW.println();
	}

	/*public String getHeaderString() {
		return new String(name + "\t" + "max force fil 1" +"\t" + "max force fil 2");
	}
	public String getDataString() {
		return new String(name + "\t" + maxForce[1] + "\t" + maxForce[2]);
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

