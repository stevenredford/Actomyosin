package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoForceIndividualActinsEvaluator extends Evaluator {


	public void init(String path, String name) throws Exception {
		super.init(path,name);
		reset();
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		myoForceMagnitudePerHeadTracker = new ValueTracker((int) (super.stopTime/interval));
		myoForcePerHeadTracker= new ValueTracker((int) (super.stopTime/interval));
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "sumMyoForces.dat");

	}
	public String getHeaderString() {
		return new String(name + "\t" + "avg force magnitude per myo head" + "\t" + "avg force per myo head" + "\t");
	}

	public String getDataString() {

		return new String(name + "\t" + myoForceMagnitudePerHeadTracker.averageVal() + "\t"+ myoForcePerHeadTracker.averageVal()+"\t");
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

	/**
	 * track intervals of time actins are not bound by myosin
	 */



	ValueTracker myoForceMagnitudePerHeadTracker;
	ValueTracker myoForcePerHeadTracker;

	public void sumMyosinSpringLengths(){
		double [] myoSpringLengths=new double [Actin.actinCt];
		double [] myoForces=new double [Actin.actinCt];
		for (int i=0; i<Actin.actinCt; i++){
			double myosinForceMagnitudeOnActinI=0;
			double myosinForceMagnitudePerBoundHead=0;
			double myosinForcePerBoundHead=0;
			for (Monomer m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					myoSpringLengths[i]+=m.boundMyo.getTangDist();
					myoForces[i]+=m.boundMyo.getforceMagTangSigned();
				}
			}
			myosinForceMagnitudeOnActinI=Myosin.myoStrokeSpringConstant*myoSpringLengths[i];
			if (ActinMyosinBindingEvaluator.myoHeadsBound[i] > 0){
				myosinForceMagnitudePerBoundHead=myosinForceMagnitudeOnActinI/ActinMyosinBindingEvaluator.myoHeadsBound[i];
				myosinForcePerBoundHead=myoForces[i]/ActinMyosinBindingEvaluator.myoHeadsBound[i];
			}
			myoForceMagnitudePerHeadTracker.registerValue(myosinForceMagnitudePerBoundHead);
			myoForcePerHeadTracker.registerValue(myosinForcePerBoundHead);
			dataPW.print(myosinForceMagnitudeOnActinI+"\t");
		

		}
		dataPW.println();

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

