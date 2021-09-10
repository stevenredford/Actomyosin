package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoDissParallelActinsEvaluator extends Evaluator {

	/** A histogram that keeps running track of distribution of minifilament bound lifetimes.. */
	HistogramPlus lifetimeHistogram;
	
	double [] firstAttchTime;
	
	double startEvalTime = 0;
	
	int [] numDetachEvents;

	boolean [] miniFilamentDetached;
	
	Vector [] lifetimes;
	
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
		lifetimes = new Vector[numActins];
		numDetachEvents = new int[numActins];
		firstAttchTime = new double[numActins];
		miniFilamentDetached = new boolean[numActins];
		
		for(int i = 0; i < numActins; i++) {
			lifetimes[i] = new Vector();
			numDetachEvents[i] = 0;
			firstAttchTime[i] = 0;
			miniFilamentDetached[i] = false;
			setMinifilamentDetached(i,true);
		}
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoLifetimes.dat");

	}
	public String getHeaderString() {
		String s = new String(name + "\t");
		for(int i = 0; i < Actin.actinCt; i++) {
			s = s + "Actin" + i + "\t" + "numDetachmentEvents\t" + "mean lifetime\t";
		}
		return s;
	}

	public String getDataString() {
		double maxL;
		Double l;
		String s = new String(name + "\t");
		for(int i = 0; i < numActins; i++) {
			maxL = 0;
			for(Enumeration e = lifetimes[i].elements(); e.hasMoreElements(); ) {
				l = ((Double)e.nextElement()).doubleValue();
				if(l>maxL) maxL = l;
			}
			lifetimeHistogram = new HistogramPlus(40,0,maxL,dataPath,"lifetimes",null,false,true,false,false);
			for(Enumeration e = lifetimes[i].elements(); e.hasMoreElements(); ) {
				l = ((Double)e.nextElement()).doubleValue();
				lifetimeHistogram.addValue(l);
			}
			lifetimeHistogram.writeBinsToFile(dataPW);
			s = s + "Actin" + i + "\t" + numDetachEvents[i] + "\t";
			if(numDetachEvents[i]>0)
				s = s + lifetimeHistogram.getMean() + "\t";
			else s = s + Sim2D.simulationTime + "\t";
		}
		dataPW.close();
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
					lifetimes[i].add(new Double(tm-firstAttchTime[i]));
					numDetachEvents[i]++;
					setMinifilamentDetached(i,true);
				}
			}
			else {  // bound by at least one head.
				if(miniFilamentDetached[i]) {
					firstAttchTime[i] = tm;
					setMinifilamentDetached(i,false);
				}
			}
		}
	}

	/** If minifilament is detached from an individual actin , we relax only that filament's attachments.*/
	public void setMinifilamentDetached(int i, boolean d) {
		if(miniFilamentDetached[i] == d) return;
		
		miniFilamentDetached[i] = d;
		
		Actin.theActins[i].restoreAttachmentPositions();
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

