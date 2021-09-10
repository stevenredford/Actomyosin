package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class ABMyoForceEvaluator extends Evaluator {


	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		avgAMyoForce=0;
		avgBMyoForce=0;
		numAForceValues=0;
		numBForceValues=0;
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "sumMyoForces.dat");

	}
	public String getHeaderString() {
		return new String(name + "\t" + "avg force on left myos" + "\t" + "avg force on right myos" + "\t");
	}

	public String getDataString() {
		avgAMyoForce/=numAForceValues;
		avgBMyoForce/=numBForceValues;
		return new String(name + "\t" + avgAMyoForce + "\t"+ avgBMyoForce+"\t");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();

		} catch (IOException ioe) { System.out.println("MyoForceIndividualActinsEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		sumMyosinSpringLengths();


	}

	double avgAMyoForce=0;
	double avgBMyoForce=0;
	int	numAForceValues=0;
	int	numBForceValues=0;


	public void sumMyosinSpringLengths(){

		double aMyoForces=0;
		double bMyoForces=0;
		for (int i=0;i<Myosin.myosinCt;i++){
			if (Myosin.theMyosins[i].boundMon != null){
				if (Myosin.theMyosins[i].iAmAnAMyo){
					aMyoForces+=Myosin.theMyosins[i].getforceMagTangSigned();
					avgAMyoForce+=Myosin.theMyosins[i].getforceMagTangSigned();
					numAForceValues+=1;
				}
				else{
					bMyoForces+=Myosin.theMyosins[i].getforceMagTangSigned();
					avgBMyoForce+=Myosin.theMyosins[i].getforceMagTangSigned();
					numBForceValues+=1;
				}
			}
		}


		dataPW.println(aMyoForces+"\t"+bMyoForces+"\t");


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

