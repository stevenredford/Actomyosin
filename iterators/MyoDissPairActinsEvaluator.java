package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoDissPairActinsEvaluator extends Evaluator {

	/** A histogram that keeps running track of distribution of minifilament bound lifetimes.. */
	HistogramPlus lifetimeHistogram;
	
	double firstAttchTime;
	
	double startEvalTime = 0;
	
	int numDetachEvents;

	boolean miniFilamentDetached;
	
	Vector lifetimes;
	
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
		lifetimes = new Vector();
		numDetachEvents = 0;
		firstAttchTime = 0;
		miniFilamentDetached = false;
		setMinifilamentDetached(true);
		
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoActinPairAttachLifetimes.dat");

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
		maxL = 0;
		for(Enumeration e = lifetimes.elements(); e.hasMoreElements(); ) {
			l = ((Double)e.nextElement()).doubleValue();
			if(l>maxL) maxL = l;
		}
		lifetimeHistogram = new HistogramPlus(40,0,maxL,dataPath,"lifetimes",null,false,true,false,false);
		for(Enumeration e = lifetimes.elements(); e.hasMoreElements(); ) {
			l = ((Double)e.nextElement()).doubleValue();
			lifetimeHistogram.addValue(l);
		}
		lifetimeHistogram.writeBinsToFile(dataPW);
		s = s + numDetachEvents + "\t";
		if(numDetachEvents>0)
			s = s + lifetimeHistogram.getMean() + "\t";
		else s = s + Sim2D.simulationTime + "\t";
		dataPW.close();
		return s;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
		} catch (IOException ioe) { System.out.println("MyoDissPairActinsEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		if(tm < startEvalTime) return;
		
		boolean miniFilamentIsBound = true;
		boolean miniFilamentFullyDetached = true;
		double myoHeadsBound;
		Monomer m;
		boolean [] bound = new boolean[numActins];
		
		for(int i = 0; i < numActins; i++) {
			myoHeadsBound = 0;
			for (m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					myoHeadsBound+=1;
					bound[i] = true;
					miniFilamentFullyDetached = false;
				}
			}
			if(myoHeadsBound == 0) {
				miniFilamentIsBound = false;
				if(!miniFilamentDetached) {
					lifetimes.add(new Double(tm-firstAttchTime));
					numDetachEvents++;
					setMinifilamentDetached(true);
				}
			}
		}
		if(miniFilamentIsBound) {  // then it mut be bound to both filaments...
			if(miniFilamentDetached) {
				firstAttchTime = tm;
				setMinifilamentDetached(false);
			}
		}
		else if(miniFilamentFullyDetached) {
			for(int i = 0; i < Actin.actinCt; i++) {
				Actin.theActins[i].restoreAttachmentPositions();
			}
			MyosinMiniFilament.theMiniFilaments[0].restoreInitialPosition();
		}
		else { // one end detached)
			for(int i = 0; i < Actin.actinCt; i++) {
				if(bound[i]) {
					Actin.theActins[i].relaxAllAttachments();
				}
				else {
					Actin.theActins[i].restoreAttachmentPositions();
				}
					
			}
		}
			
	}
	
	public void setMinifilamentDetached(boolean d) {
		if(miniFilamentDetached == d) return;
		
		miniFilamentDetached = d;
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

