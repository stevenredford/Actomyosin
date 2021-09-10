package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoForceManyActinsEvaluator extends Evaluator {

	static public String
		FILAMENT = "Filament",
		MINIFILAMENT = "Minifilament";
	
	double sampleAverageInterval;
	
	double  nextSampleTime;
	
	int intervalCt;
	
	double 	aveDutyRatio;
	
	double [] cumSpringForces;
	
	double [] cumMyoHeadsBound;

	double cumVelocity;

	double position;
	
	double lastPosition;
	
	double finalDutyRatio;
	
	double finalBoundVelocity;
	
	double finalVelocity;

	double [] finalForcePerHead;
					
	double [] finalMyoForces;
					
	double [] finalMyoHeadsBound;
	
	boolean firstTime;
	
	double startAverageTime = 0;
	
	int averageCt;

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
		
		intervalCt = averageCt = 0;
		firstTime = true;
		nextSampleTime = sampleAverageInterval;
		aveDutyRatio = finalDutyRatio = 0;
		cumSpringForces = new double[numActins];
		cumMyoHeadsBound = new double[numActins];

		cumVelocity = 0;
		position = 0;
		lastPosition = 0;
		
		finalForcePerHead = new double[numActins];
		finalMyoForces = new double[numActins];
		finalMyoHeadsBound = new double[numActins];
		
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoForces.dat");
	}
	
	public void writeDataFileHeader() {
		dataPW.print("dutyRatio\tcumVelocity\t");
		for(int i = 0; i < numActins; i++) {
			dataPW.print("Actin " + i + "\t" + "cumMyoHeadsBound\t" + "cumMyoForces\t"+ "forcePerHead\t");
		}
		dataPW.println();
	}
	
	public String getHeaderString() {
		String s = new String(name + "\t" + "average duty ratio\t" + "velocity\t");
		for(int i = 0; i < Actin.actinCt; i++) {
			s = s +  "Actin" + i + "\tcumMyoHeadsBound" + "\tcumMyoForces"+ "\tforcePerHead\t";
		}
		return s;
	}

	public String getDataString() {
		String s = new String(name + "\t");
		if(averageCt > 0) {
			finalDutyRatio/=averageCt;
			finalVelocity/=averageCt;
			for(int i = 0; i < numActins; i++) {
				if (averageCt > 0) {
					finalMyoHeadsBound[i]/=averageCt;
					finalMyoForces[i]/=averageCt;
					finalForcePerHead[i]/=averageCt;
				}
			}
		}else {
			finalDutyRatio = finalVelocity = 0;
			for(int i = 0; i < numActins; i++) {
				finalMyoHeadsBound[i]=finalMyoForces[i]=finalForcePerHead[i]=0;
			}
		}

		s = s + finalDutyRatio + "\t" + finalVelocity + "\t";
		for(int i = 0; i< numActins; i++) {
			s = s +  "Actin" + i;
			s = s  + "\t" + finalMyoHeadsBound[i];
			s = s + "\t" + finalMyoForces[i];
			s = s + "\t" + finalForcePerHead[i];
			s = s + "\t";
		}
		return s;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();

		} catch (IOException ioe) { System.out.println("MyoForceIndividualActinsEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		binMyosinForces(tm);
	}

	public void binMyosinForces(double tm){
		
		double springForces;
		double myoHeadsBound;
		double forcePerHead;
		Monomer mon;
		Actin a;
		MyosinMiniFilament m = MyosinMiniFilament.theMiniFilaments[0];
		FilamentAttachment bEnd, pEnd;
		
		intervalCt++;
		
		aveDutyRatio+=Myosin.getAveDutyRatio();
		position = m.cm.x;
		cumVelocity += Point2D.getDiffX(position,lastPosition);
		lastPosition = position;
		
		if(tm > nextSampleTime) {
			aveDutyRatio/=intervalCt;
			cumVelocity/=intervalCt;
			
			if(tm > startAverageTime) {
				finalDutyRatio+=aveDutyRatio;
				finalBoundVelocity+=cumVelocity;
			}
			dataPW.print(aveDutyRatio + "\t" + cumVelocity + "\t");
		}
	
		for(int i = 0; i < numActins; i++) {
			a = Actin.theActins[i];
			bEnd= a.getBEndAttachment();
			pEnd = a.getPEndAttachment();
			springForces = bEnd.getPinForce().x + pEnd.getPinForce().x;
			
			myoHeadsBound = 0;
			for (mon=a.pEndMonomer; mon!=null; mon=mon.next){
				if (!mon.isFreeMyosin()){
					myoHeadsBound+=1;
				}
			}
			cumMyoHeadsBound[i] += myoHeadsBound;

			double d = Point2D.getDiffX(a.cm.x,m.cm.x);
			double dL = 0.5*(a.length - m.totalLength);
			if(a.uVect.x>0) {
				if(d + dL < 100) {
					a.treadmill(30);
				}
			}
			else {
				if(dL-d < 100) {
					a.treadmill(30);
				}
			}

			cumSpringForces[i] += springForces;
				
			if(tm > nextSampleTime) {
				forcePerHead = cumSpringForces[i]/cumMyoHeadsBound[i];
				cumSpringForces[i]/=intervalCt;
				cumMyoHeadsBound[i]/=intervalCt;
				
				if(tm > startAverageTime) {
					finalForcePerHead[i]+=forcePerHead;
					finalMyoForces[i]+=cumSpringForces[i];
					finalMyoHeadsBound[i]+=cumMyoHeadsBound[i];
				}
				dataPW.print("Actin " + i + "\t" + cumMyoHeadsBound[i] + "\t" + cumSpringForces[i] +"\t"+ forcePerHead +"\t");
				cumSpringForces[i] = 0;
				cumMyoHeadsBound[i] = 0;
			}
		}
		
		if(tm > nextSampleTime) {
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

