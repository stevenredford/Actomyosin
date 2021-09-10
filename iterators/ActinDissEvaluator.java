package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class ActinDissEvaluator extends Evaluator {


	double firstAttchTime;
	
	double firstDetachTime;
	
	double startEvalTime = 0;
	
	int numDetachEvents;

	boolean filamentDetached = true;
		
	Vector attachLifetimes, detachLifetimes;
	

	public void init(String path, String name) throws Exception {
		super.init(path,name);
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		attachLifetimes = new Vector();
		detachLifetimes = new Vector();
		numDetachEvents = 0;
		firstAttchTime = -1;
		firstDetachTime = -1;
		filamentDetached = true;
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoLifetimes.dat");

	}
	public String getHeaderString() {
		return new String("numDetachmentEvents\t" + "mean attached lifetime\t" + "mean detached lifetime\t");
	}

	public String getDataString() {
		double meanL = 0;
		int cnt = 0;
		Double l;
		attachLifetimes.add(new Double(Sim2D.simulationTime-firstAttchTime));
		for(Enumeration e = attachLifetimes.elements(); e.hasMoreElements(); ) {
			l = ((Double)e.nextElement()).doubleValue();
			dataPW.print(l + "\t");
			meanL += l;
			cnt++;
		}
		dataPW.println();
		meanL/=cnt;

		String s = new String(numDetachEvents + "\t");
		if(numDetachEvents>0)
			s = s + meanL + "\t";
		else s = s + Sim2D.simulationTime + "\t";

		meanL = 0; cnt = 0;
		detachLifetimes.add(new Double(Sim2D.simulationTime-firstAttchTime));
		for(Enumeration e = detachLifetimes.elements(); e.hasMoreElements(); ) {
			l = ((Double)e.nextElement()).doubleValue();
			dataPW.print(l + "\t");
			meanL += l;
			cnt++;
		}
		dataPW.println();
		meanL/= cnt;
		dataPW.close();
		if(numDetachEvents>0)
			s = s + meanL + "\t";
		else s = s + Sim2D.simulationTime + "\t";
		return s;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
		} catch (IOException ioe) { System.out.println("ActinDissEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		if(tm < startEvalTime) return;
		boolean myoHeadsBound = false;
		
		for (Monomer m=Actin.theActins[0].pEndMonomer; m!=null; m=m.next){
			if (!m.isFreeMyosin()){
				myoHeadsBound=true;
			}
		}
		
		if(myoHeadsBound) {
			if(!filamentDetached) {
				firstDetachTime = tm;
				attachLifetimes.add(new Double(tm-firstAttchTime));
				numDetachEvents++;
				setFilamentDetached(true);
			}
		}
		else {  // bound by at least one head.
			if(filamentDetached) {
				firstAttchTime = tm;
				if(firstDetachTime > 0) detachLifetimes.add(new Double(tm-firstDetachTime));
				setFilamentDetached(false);
			}
		}
	}

	public void setFilamentDetached(boolean d) {
		if(filamentDetached == d) return;
		filamentDetached = d;
	}


	/** Override this method to write out a datafile containing Evaluator subtype-specific
	 * data in whatever format you choose.  The method should return true if the dataFie was successfuly written,
	 * false otherwise.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equalsIgnoreCase("startEvalTime"))  {
			startEvalTime = in.nextDouble();
		}
		else super.loadParameter(tag, in);
	}
}

