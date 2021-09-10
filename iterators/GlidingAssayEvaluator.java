package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;

public class GlidingAssayEvaluator extends Evaluator {
	
	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus forceHistogram;

	PrintWriter dataPW2;

	String dataFileName2;

	double sampleAverageInterval;
	
	double  nextSampleTime;
	
	int intervalCt;
	
	double position;
	
	double lastPosition;
	
	double 	aveDutyRatio;

	double finalVelocity;

	double finalMyoHeadsBound;
					
	double finalDutyRatio;
	
	boolean firstTime;
	
	double startAverageTime = 0;
	
	double averageCt;

	double cumMyoHeadsBound;
	
	double cumVelocity;
	
	double firstPosition;
	
	boolean firstAverage;
	
	public void init(String path, String name) throws Exception {
		super.init(path,name);
	}

	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		forceHistogram.clearBins();
		intervalCt = 0;
		averageCt=0;
		cumMyoHeadsBound = 0;
		cumVelocity = 0;
		aveDutyRatio = 0;
		finalVelocity = 0;
		finalMyoHeadsBound = 0;
		finalDutyRatio = 0;
		position = 0;
		lastPosition = 0;
		firstTime = true;
		firstAverage = true;
		nextSampleTime = sampleAverageInterval;
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "gliding.dat");
		dataFileName2 = new String(dataPath + File.separator + "springForceHistogram.dat");
		forceHistogram = new HistogramPlus(40,-10,10,dataPath,"myoStrains",null,false,true,false,false);
	}

	
	/** Override this method to specify whether this Evaluator writes data files.  These
	 * files MUST be written to the path defined by the variable dataPata to be compiled
	 * by the iterator.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}


	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			dataPW2 = new PrintWriter(new FileWriter(new File (dataFileName2)),true);
			writeDataFileHeader();
		} catch (IOException ioe) { System.out.println("GlidingAssayEvaluator.mkDataFile(): error creating File"); }
	}

	/** Subclasses override to actually write file header */
	public void writeDataFileHeader() {
		dataPW.println("aveDutyRatio" + "\t" + "cumMyoHeadsBound" + "\t" + "cumVelocity");
	}

	public String getHeaderString() {
		String s = new String();
		s = s + "myoDensity\t" + "myoLineDensity\t" + "ViscosityMultiplier\t" + "GardelNumber\t" + "ResistiveLoad\t" + "velocity\t";
		return s;
	}

	public String getDataString() {
		String s = new String();
		s+= MyosinSurface.theMyoSurfaces[0].myoDensity + "\t";
		s+= MyosinSurface.theMyoSurfaces[0].myoLineDensity + "\t";
		s+= Constants.xWater + "\t";
		s+= Actin.getNoNameRatio() + "\t";
		s+= Actin.resistiveLoad +"\t";
		s+= finalVelocity + "\t";
		return s;
	}

	public void evaluate(double tm)  {
		evaluateGlide(tm);
	}

	public void evaluateGlide (double tm) {
		double myoHeadsBound = 0;
		
		intervalCt++;

		for (Monomer m=Actin.theActins[0].pEndMonomer; m!=null; m=m.next){
			if (!m.isFreeMyosin()){
				myoHeadsBound+=1;
				forceHistogram.addValue(m.boundMyo.getforceMagTangSigned());
			}
		}

		cumMyoHeadsBound += myoHeadsBound;
		aveDutyRatio+=Myosin.getAveDutyRatio();
		
		position = Actin.theActins[0].cm.x;
		cumVelocity += position-lastPosition;
		lastPosition = position;

		if(tm > nextSampleTime) {
			forceHistogram.writeBinsToFile(dataPW2);
			forceHistogram.clearBins();

			aveDutyRatio/=intervalCt;
			cumMyoHeadsBound/=intervalCt;
			cumVelocity/=sampleAverageInterval;
			
			if(tm > startAverageTime) {
				if(firstAverage) {
					firstPosition = position;
					firstAverage = false;
				}
				finalMyoHeadsBound+=cumMyoHeadsBound;
				finalDutyRatio+=aveDutyRatio;
				finalVelocity+=cumVelocity;
				averageCt++;
			}
			
			dataPW.println(aveDutyRatio + "\t" + cumMyoHeadsBound + "\t" + cumVelocity);
			intervalCt = 0;
			cumMyoHeadsBound = 0;
			cumVelocity = 0;
			aveDutyRatio = 0;
			nextSampleTime += sampleAverageInterval;
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

	public void cleanUp()  {}
				

}
