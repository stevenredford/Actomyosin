package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoDissManyActinsEvaluator extends Evaluator {

	
	double [] firstAttachTime;
	
	double [] firstDetachTime;
	
	boolean fastFilamentRelaxation = false;
	
	Vector [] attachLifetimes, detachLifetimes;

	double startEvalTime = 0;
	
	int [] numDetachEvents;

	boolean [] miniFilamentDetached;
	
	int numActins;
	
	public void init(String path, String name) throws Exception {
		super.init(path,name);
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		numActins = Actin.actinCt;
		numDetachEvents = new int[numActins];
		firstAttachTime = new double[numActins];
		firstDetachTime = new double[numActins];
		miniFilamentDetached = new boolean[numActins];
		attachLifetimes = new Vector[numActins];
		detachLifetimes = new Vector[numActins];
		
		for(int i = 0; i < numActins; i++) {
			numDetachEvents[i] = 0;
			firstAttachTime[i] = firstDetachTime[i] = 0;
			attachLifetimes[i] = new Vector();
			detachLifetimes[i] = new Vector();
			miniFilamentDetached[i] = true;
		}
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoLifetimes.dat");

	}
	public String getHeaderString() {
		String s = new String(name + "\t");
		for(int i = 0; i < Actin.actinCt; i++) {
			s = s + "Actin" + i + "\t" + "numDetachmentEvents\t" + "mean attached lifetime\t" + "lastFirstAttachedTime\t" + "mean detached lifetime\t";
		}
		return s;
	}

	public String getDataString() {
		double meanL;
		int cnt;
		Double l;
		String s = new String(name + "\t");
		for(int i = 0; i < numActins; i++) {
			meanL=cnt=0;
			if(!miniFilamentDetached[i]) {
				attachLifetimes[i].add(new Double(Sim2D.simulationTime-firstAttachTime[i]));
			}
			for(Enumeration e = attachLifetimes[i].elements(); e.hasMoreElements(); ) {
				l = ((Double)e.nextElement()).doubleValue();
				dataPW.print(l + "\t");
				meanL += l;
				cnt++;
			}
			dataPW.println();
			if(cnt>0) meanL/=cnt;
			s = s + numDetachEvents + "\t";
			if(numDetachEvents[i]>0)
				s = s + meanL + "\t";
			else s = s + Sim2D.simulationTime + "\t";
			
			s = s + firstAttachTime[i] + "\t";
	
			meanL = cnt = 0;
			for(Enumeration e = detachLifetimes[i].elements(); e.hasMoreElements(); ) {
				l = ((Double)e.nextElement()).doubleValue();
				dataPW.print(l + "\t");
				meanL += l;
				cnt++;
			}
			dataPW.println();
			meanL/= cnt;
			dataPW.close();
			
			if(numDetachEvents[i]>0)
				s = s + meanL + "\t";
			else s = s + "0" + "\t";
		}
		return s;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
		} catch (IOException ioe) { System.out.println("MyoDissEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		if(tm < startEvalTime) return;
		
		double myoHeadsBound;
		Monomer m;
		
		for(int i = 0; i < numActins; i++) {
			myoHeadsBound = 0;
			for (m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					myoHeadsBound+=1;
				}
			}
			if(myoHeadsBound == 0) {
			if(!miniFilamentDetached[i]) {
				firstDetachTime[i] = tm;
				double attachTm = tm-firstAttachTime[i];
				attachLifetimes[i].add(new Double(attachTm));
				numDetachEvents[i]++;
				setMinifilamentDetached(i,true);
			}
		}
		else {  // bound by at least one head.
			if(miniFilamentDetached[i]) {
				firstAttachTime[i] = tm;
				if(firstDetachTime[i] > 0) detachLifetimes[i].add(new Double(tm-firstDetachTime[i]));
				setMinifilamentDetached(i,false);
			}
		}
		}
	}

	/** If minifilament is detached, we relax BOTH filaments.*/
	public void setMinifilamentDetached(int i, boolean d) {
		if(miniFilamentDetached[i] == d) return;
		miniFilamentDetached[i] = d;
		
		if(d) {
			if(fastFilamentRelaxation) {
				Actin.theActins[i].restoreAttachmentPositions();
			}
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
	public void loadParameter(String tag, AMInputStream in)  throws Exception {

		if(tag.equalsIgnoreCase("startEvalTime"))  {
			startEvalTime = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("fastFilamentRelaxation")) {
			fastFilamentRelaxation = in.nextBoolean();
		}

		else super.loadParameter(tag, in);

	}


}

