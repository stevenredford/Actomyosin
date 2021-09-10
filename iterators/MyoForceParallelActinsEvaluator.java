package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoForceParallelActinsEvaluator extends Evaluator {

	PrintWriter dataPW2;

	String dataFileName2;

	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus [] forceHistogram;

	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus [] cumForceHistogram;
	
	double  nextSampleTime;
	
	int intervalCt;
	
	double [] cumMyoForces;
	
	double [] cumPosition;
	
	double [] cumMyoHeadsBound;

	double [] lastPosition;
	
	double [] actinPositions;
	
	double [] finalForcePerHead;
					
	double [] finalMyoForces;
					
	double [] finalMyoHeadsBound;
	
	double [] avgForceSquared;
					
	boolean firstTime;
	
	double startAverageTime = 0;
	
	double lastAverageTime = -1;
	
	double sampleAverageInterval = 1e20;
	
	int numActins;
		
	int averageCt;

	public void init(String path, String name) throws Exception {
		super.init(path,name);
	}


	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		
		numActins = Actin.actinCt;
		intervalCt = averageCt = 0;
		firstTime = true;
		lastAverageTime = -1;
		nextSampleTime = sampleAverageInterval;
		
		forceHistogram = new HistogramPlus[numActins];
		cumForceHistogram = new HistogramPlus[numActins];
		cumMyoForces = new double[numActins];
		cumMyoHeadsBound = new double[numActins];
		cumPosition = new double[numActins];
		lastPosition = new double[numActins];
		actinPositions = new double[numActins];
		finalForcePerHead = new double[numActins];
		finalMyoForces = new double[numActins];
		finalMyoHeadsBound = new double[numActins];
		avgForceSquared = new double[numActins];
		
		for(int i = 0; i < numActins; i++) {
			forceHistogram[i] = new HistogramPlus(40,-10,10,dataPath,"myoStrains",null,false,true,false,false);
			cumForceHistogram[i] = new HistogramPlus(40,-10,10,dataPath,"cumMyoStrains",null,false,true,false,false);
		}
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoForces.dat");
		dataFileName2 = new String(dataPath + File.separator + "springForceHistogram.dat");
	}
	
	public String getHeaderString() {
		String s = new String(name + "\t");
		for(int i = 0; i < Actin.actinCt; i++) {
			s = s +  "Actin" + i + "\tcumMyoHeadsBound" + "\tcumMyoForces"+ "\tforcePerHead" +  "\tRMSD magnitude sum\t";
		}
		return s;
	}

	public String getDataString() {
		String s = new String(name + "\t");
		double [] rMSD = new double[numActins];
		for(int i = 0; i < numActins; i++) {
			if (averageCt > 0) {
				finalMyoHeadsBound[i]/=averageCt;
				finalMyoForces[i]/=averageCt;
				finalForcePerHead[i]/=averageCt;
				avgForceSquared[i]/=averageCt;
				rMSD[i] = Math.sqrt(avgForceSquared[i] - Math.pow(finalMyoForces[i], 2));
			}
			else {
				finalMyoHeadsBound[i]=finalMyoForces[i]=finalForcePerHead[i]=rMSD[i] =0;
			}
		}

		if(lastAverageTime<0) lastAverageTime=Sim2D.simulationTime;
		for(int i = 0; i< numActins; i++) {
			s = s +  "Actin" + i;
			s = s  + "\t" + finalMyoHeadsBound[i];
			s = s + "\t" + finalMyoForces[i];
			s = s + "\t" + finalForcePerHead[i];
			s = s + "\t" + rMSD[i];
			s = s + "\t";
		}
		return s;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			dataPW2 = new PrintWriter(new FileWriter(new File (dataFileName2)),true);
			writeDataFileHeader();

		} catch (IOException ioe) { System.out.println("MyoForceIndividualActinsEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		binMyosinForces(tm);
	}

	public void doFirstEval() {
		for(int i = 0; i < numActins; i++) {
			actinPositions[i] = Actin.theActins[1].cm.x;
		}
	}

	public void binMyosinForces(double tm){
		
		double myoForces;
		double myoHeadsBound;
		double forcePerHead;
		Monomer m;
		Actin a;
		FilamentAttachment bEnd, pEnd;
		
		intervalCt++;
		for(int i = 0; i < numActins; i++) {
			a = Actin.theActins[i];
			myoHeadsBound = myoForces = 0;
			for (m=Actin.theActins[i].pEndMonomer; m!=null; m=m.next){
				if (!m.isFreeMyosin()){
					forceHistogram[i].addValue(m.boundMyo.getforceMagTangSigned());
					cumForceHistogram[i].addValue(m.boundMyo.getforceMagTangSigned());
					myoHeadsBound+=1;
				}
			}
			bEnd= a.getBEndAttachment();
			pEnd = a.getPEndAttachment();
			myoForces = bEnd.getPinForce().x + pEnd.getPinForce().x;
			
			cumPosition[i] += a.cm.x;
			cumMyoForces[i] += myoForces;
			cumMyoHeadsBound[i] += myoHeadsBound;

			if(tm > nextSampleTime) {
				forceHistogram[i].writeBinsToFile(dataPW2);
				forceHistogram[i].clearBins();
				forcePerHead = cumMyoForces[i]/cumMyoHeadsBound[i];
				cumMyoForces[i]/=intervalCt;
				cumMyoHeadsBound[i]/=intervalCt;
				cumPosition[i]/=intervalCt;
				
				avgForceSquared[i] += Math.pow(myoForces, 2);
				finalForcePerHead[i]+=forcePerHead;
				finalMyoForces[i]+=cumMyoForces[i];
				finalMyoHeadsBound[i]+=cumMyoHeadsBound[i];
				dataPW.print("Actin " + i + "\t" + cumMyoHeadsBound[i] + "\t" + cumMyoForces[i] +"\t"+ forcePerHead +"\t"+ cumPosition[i] + "\t");
				cumMyoForces[i] = 0;
				cumMyoHeadsBound[i] = 0;
				lastPosition[i] = cumPosition[i];
				cumPosition[i] = 0;
			}
		}
		
		if(tm > nextSampleTime) {
			averageCt++;
			dataPW.print("\n");
			dataPW.flush();
			intervalCt = 0;
			nextSampleTime+=sampleAverageInterval;
		}
	}

	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {

		if(tag.equalsIgnoreCase("startAverageTime"))  {
			startAverageTime = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("sampleAverageInterval"))  {
			sampleAverageInterval = in.nextDouble();
		}
		else super.loadParameter(tag, in);

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

