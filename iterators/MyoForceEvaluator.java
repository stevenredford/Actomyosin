package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoForceEvaluator extends Evaluator {
	
	static public String
		FILAMENT = "Filament",
		MINIFILAMENT = "Minifilament";
	
	PrintWriter dataPW2;
	PrintWriter dataPW3;
	PrintWriter	dataPW4;
	PrintWriter	dataPW5;
	PrintWriter	dataPW6;

	
	String dataFileName2;
	String dataFileName3;
	String dataFileName4;
	String dataFileName5;
	String dataFileName6;

	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus forceHistogram;

	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus unnormalizedForceHistogram;
	
	/*histograms that separately keep track of forces on motors with different adp release rates*/
	HistogramPlus motorType1ForceHistogram;
	HistogramPlus motorType2ForceHistogram;
	
	/**a force histogram that does not average over sampleAverageInterval*/
	HistogramPlus instantaneousForceHistogram;
	
	double sampleAverageInterval;
	
	double  nextSampleTime;
	
	int intervalCt;
	
	double position;
	
	double lastPosition;
	
	double aveDutyRatio;
		
	double finalBoundVelocity;
	
	double finalVelocity;
					
	double finalForcePerHead;
					
	double finalMyoForces;
					
	double finalMyoHeadsBound;
					
	double finalDutyRatio;
	
	double startAverageTime = 0;
	
	int averageCt;

	boolean miniFilamentBound = false;
	
	boolean miniFilamentDetached = false;
	
	boolean recordActinPositions=false;
	
	double cumMyoHeadsBound;
	
	double cumSpringForces;
	
	double cumCrossbridgeForceSum;
	
	double cumBoundVelocity;
	
	double cumVelocity;
	
	double cumForcePerHead;
	
	double nBoundSteps;
	
	double totalBoundSteps;
	
	double meanNegativeForce=0;
	
	double meanPositiveForce=0;
	
	double meanMotorType1NegativeForce=0;
	
	double meanMotorType2NegativeForce=0;
	
	
	double stDevBoundVelocity=0;
	
	int stDevCt=0;
	

	
	String velocityTrackingTarget = MINIFILAMENT;
	
	public void init(String path, String name) throws Exception {
		super.init(path,name);
	}

	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		super.reset();
		forceHistogram.clearBins();
		unnormalizedForceHistogram.clearBins();
		motorType1ForceHistogram.clearBins();
		motorType2ForceHistogram.clearBins();
		instantaneousForceHistogram.clearBins();
		intervalCt = averageCt = 0;
		aveDutyRatio = finalDutyRatio = 0;
		cumSpringForces = 0;
		cumCrossbridgeForceSum=0;
		cumMyoHeadsBound = 0;
		cumForcePerHead = 0;
		cumBoundVelocity = 0;
		cumVelocity = 0;
		nBoundSteps = 0;
		totalBoundSteps = 0;
		aveDutyRatio = 0;
		finalVelocity = 0;
		finalBoundVelocity = 0;
		finalForcePerHead = 0;
		finalMyoForces = 0;
		finalMyoHeadsBound = 0;
		position = 0;
		lastPosition = 0;
		stDevBoundVelocity=0;
		stDevCt=0;
		miniFilamentDetached = false;
		miniFilamentBound = false;
		nextSampleTime = sampleAverageInterval;
		MyosinMiniFilament.turnResistiveLoadOff();
	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoForces.dat");
		dataFileName2 = new String(dataPath + File.separator + "springForceHistogram.dat");
		dataFileName3 = new String(dataPath + File.separator + "unnormalizedSpringForceHistogram.dat");
		dataFileName4 = new String(dataPath + File.separator + "motorType1SpringForceHistogram.dat");
		dataFileName5 = new String(dataPath + File.separator + "motorType2SpringForceHistogram.dat");
		dataFileName6 = new String(dataPath + File.separator + "instantaneousSpringForceHistogram.dat");
		forceHistogram = new HistogramPlus(40,-10,10,dataPath,"myoStrains",null,false,true,false,false);
		motorType1ForceHistogram = new HistogramPlus(40,-10,10,dataPath,"myoType1Strains",null,false,true,false,false);
		motorType2ForceHistogram = new HistogramPlus(40,-10,10,dataPath,"myoType2Strains",null,false,true,false,false);
		unnormalizedForceHistogram = new HistogramPlus(40,-10,10,dataPath,"cumMyoStrains",null,true,false,false,false);
		instantaneousForceHistogram = new HistogramPlus(40,-10,10,dataPath,"instantaneousMyoStrains",null,true,false,false,false);
	}
	public String getHeaderString() {
		String s = new String();
		s = s + "NMyosinHeads\t" + "resistiveLoad\t" + "aveDutyRatio\t" + "cumMyoHeadsBound\t" + "cumMyoForces\t"+ "forcePerHead\t" + "velocity\t" + "boundVelocity\t" + "stDevVelocity\t"+"meanNegativeForceMag\t"+"meanPositiveForceMag\t"+"meanNegativeForceMagMotor1\t"+"meanNegativeForceMagMotor2\t";
		return s;
	}

	public String getDataString() {
		String s = new String();
		if (averageCt > 0) {
			finalDutyRatio/=averageCt;
			finalMyoHeadsBound/=averageCt;
			finalMyoForces/=averageCt;
			finalForcePerHead/=averageCt;
			finalVelocity/=averageCt;
			finalBoundVelocity/=averageCt;
			meanNegativeForce/=averageCt;
			meanPositiveForce/=averageCt;
			meanMotorType1NegativeForce/=averageCt;
			meanMotorType2NegativeForce/=averageCt;
			stDevBoundVelocity/=stDevCt;
			stDevBoundVelocity=stDevBoundVelocity-Math.pow(finalVelocity, 2);
			stDevBoundVelocity=Math.sqrt(stDevBoundVelocity);
		}
		else 	finalDutyRatio=finalMyoHeadsBound=finalMyoForces=finalForcePerHead=finalVelocity=finalBoundVelocity=0;
		s = s + MyosinMiniFilament.nMyosinHeads + "\t" + MyosinMiniFilament.theMiniFilaments[0].sizeOfResistiveLoad() + "\t" + finalDutyRatio + "\t" + finalMyoHeadsBound + "\t" + finalMyoForces + "\t"+
			finalForcePerHead + "\t" + finalVelocity + "\t" + finalBoundVelocity + "\t"+stDevBoundVelocity+"\t"+meanNegativeForce+"\t"+meanPositiveForce+"\t"+meanMotorType1NegativeForce+"\t"+meanMotorType2NegativeForce;
		return s;
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			dataPW2 = new PrintWriter(new FileWriter(new File (dataFileName2)),true);
			dataPW3 = new PrintWriter(new FileWriter(new File (dataFileName3)),true);
			dataPW4 = new PrintWriter(new FileWriter(new File (dataFileName4)),true);
			dataPW5 = new PrintWriter(new FileWriter(new File (dataFileName5)),true);
			dataPW6 = new PrintWriter(new FileWriter(new File (dataFileName6)),true);

			writeDataFileHeader();

		} catch (IOException ioe) { System.out.println("MyoForceIndividualActinsEvaluator.mkDataFile(): error creating File"); }
	}

	/** Subclasses override to actually write file header */
	public void writeDataFileHeader() {
		dataPW.println("aveDutyRatio" + "\t" + "cumMyoHeadsBound" + "\t" + "cumSpringForces" +"\t"+ "instantaneousSpringForce"+"\t"+"cumForcePerHead" +"\t"+ "cumBoundVelocity" + "\t" + "cumVelocity"+"\t"+"actinPosition"+"\t"+"cumCrossbridgeForceSum"+"\t"+"negativeCrossbridgeForceMagnitude"+"\t"+"positiveCrossbridgeForceMagnitude"+"\t"+"negativeCrossbridgeForceMagnitudeMotor1"+"\t"+"negativeCrossbridgeForceMagnitudeMotor2"+"\t"+"instantaneousVelocity");
	}

	public void evaluate(double tm)  {
		binMyosinForces(tm);
	}

	public void doFirstEval() {	}

	public void binMyosinForces(double tm){
		
		double instantaneousVelocity=0;
		double springForces = 0;
		double crossbridgeForceSum=0;
		double myoHeadsBound = 0;
		double vOffset = 0;
		double v;
		
		intervalCt++;
		double negativeForceMagnitude=0;
		double positiveForceMagnitude=0;
		double negativeForceMagnitudeMotorType1=0;
		double negativeForceMagnitudeMotorType2=0;

		
		for (Monomer m=Actin.theActins[0].pEndMonomer; m!=null; m=m.next){
			if (!m.isFreeMyosin()){
				forceHistogram.addValue(m.boundMyo.getforceMagTangSigned());
				unnormalizedForceHistogram.addValue(m.boundMyo.getforceMagTangSigned());
				if (tm>nextSampleTime){
					instantaneousForceHistogram.addValue(m.boundMyo.getforceMagTangSigned());
				}
				myoHeadsBound+=1;
				crossbridgeForceSum+=m.boundMyo.getforceMagTangSigned();
				if (m.boundMyo.usedReleaseRateMult==Myosin.adpReleaseRateMult) {
					motorType1ForceHistogram.addValue(m.boundMyo.getforceMagTangSigned());
					if (m.boundMyo.getforceMagTangSigned()<0){
						negativeForceMagnitudeMotorType1=negativeForceMagnitudeMotorType1+Math.abs(m.boundMyo.getforceMagTangSigned());
					}

				}
				else if (m.boundMyo.usedReleaseRateMult==Myosin.secondADPReleaseRateMult){
					motorType2ForceHistogram.addValue(m.boundMyo.getforceMagTangSigned());
					if (m.boundMyo.getforceMagTangSigned()<0){
						negativeForceMagnitudeMotorType2=negativeForceMagnitudeMotorType2+Math.abs(m.boundMyo.getforceMagTangSigned());
					}

				}
				if (m.boundMyo.getforceMagTangSigned()<0){
					negativeForceMagnitude=negativeForceMagnitude+Math.abs(m.boundMyo.getforceMagTangSigned());
				}
				else{
					positiveForceMagnitude=positiveForceMagnitude+Math.abs(m.boundMyo.getforceMagTangSigned());
				}
			}
		}
		
		
		springForces = Actin.theActins[0].getBEndAttachment().getPinForce().x + Actin.theActins[0].getPEndAttachment().getPinForce().x;
		cumMyoHeadsBound += myoHeadsBound;
		aveDutyRatio+=Myosin.getAveDutyRatio();
		
		if(myoHeadsBound > 0) {
			
			double d = Point2D.getDiffX(Actin.theActins[0].cm.x,MyosinMiniFilament.theMiniFilaments[0].cm.x);
			double dL = 0.5*(Actin.theActins[0].length - MyosinMiniFilament.theMiniFilaments[0].totalLength);
			if(d + dL < 100 && !recordActinPositions) {
				Actin.theActins[0].treadmill(30);
				vOffset = 30*Actin.monLength;
			}
			
			if(velocityTrackingTarget.equals(MINIFILAMENT)) {
				position = MyosinMiniFilament.theMiniFilaments[0].cm.x;
			}
			else {
				position = Actin.theActins[0].cm.x;
			}

			if(miniFilamentDetached) {
				miniFilamentDetached = false;
				MyosinMiniFilament.turnResistiveLoadOn();
			}
			else {
				v = Point2D.getDiffX(position,lastPosition);
				if(velocityTrackingTarget.equals(FILAMENT)) {
					 v -= vOffset;
				}
				instantaneousVelocity=v/interval;
				stDevBoundVelocity+=Math.pow(instantaneousVelocity,2);
				stDevCt++;
				cumBoundVelocity += v;
				cumForcePerHead += springForces/myoHeadsBound;
				cumSpringForces += springForces;
				cumCrossbridgeForceSum+=crossbridgeForceSum;
				nBoundSteps+=1;
			}
		}
		else {
			if(!miniFilamentDetached) {
				cumSpringForces += springForces;
				cumCrossbridgeForceSum+=crossbridgeForceSum;
				miniFilamentDetached = true;
				MyosinMiniFilament.turnResistiveLoadOff();
			}
		}
		lastPosition = position;

		if(tm > nextSampleTime) {
			forceHistogram.writeBinsToFile(dataPW2);
			forceHistogram.clearBins();
			unnormalizedForceHistogram.writeBinsToFile(dataPW3);
			unnormalizedForceHistogram.clearBins();
			motorType1ForceHistogram.writeBinsToFile(dataPW4);
			motorType1ForceHistogram.clearBins();
			motorType2ForceHistogram.writeBinsToFile(dataPW5);
			motorType2ForceHistogram.clearBins();
			instantaneousForceHistogram.writeBinsToFile(dataPW6);
			instantaneousForceHistogram.clearBins();
			aveDutyRatio/=intervalCt;
			cumMyoHeadsBound/=intervalCt;
			cumCrossbridgeForceSum/=intervalCt;
			if(nBoundSteps > 0) {
				cumForcePerHead/=nBoundSteps;
				cumSpringForces/=nBoundSteps;
				cumBoundVelocity/=(nBoundSteps*interval);
			}
			cumVelocity = nBoundSteps*cumBoundVelocity/intervalCt;
			
			if(tm > startAverageTime) {
				finalForcePerHead+=cumForcePerHead;
				finalMyoForces+=cumSpringForces;
				finalMyoHeadsBound+=cumMyoHeadsBound;
				finalDutyRatio+=aveDutyRatio;
				finalVelocity+=cumVelocity;
				finalBoundVelocity+=cumBoundVelocity;
				meanNegativeForce+=negativeForceMagnitude;
				meanPositiveForce+=positiveForceMagnitude;
				meanMotorType1NegativeForce+=negativeForceMagnitudeMotorType1;
				meanMotorType2NegativeForce+=negativeForceMagnitudeMotorType2;
				averageCt++;
			}
			double actinPosition=0;
			if(recordActinPositions){
				actinPosition=Actin.theActins[0].cm.x;
			}
			
			
			dataPW.print(aveDutyRatio + "\t" + cumMyoHeadsBound + "\t" + cumSpringForces +"\t"+ springForces+"\t"+cumForcePerHead +"\t"+ cumBoundVelocity + "\t" + cumVelocity +"\t"+ actinPosition+"\t"+cumCrossbridgeForceSum+"\t"+negativeForceMagnitude+"\t"+positiveForceMagnitude+"\t"+negativeForceMagnitudeMotorType1+"\t"+negativeForceMagnitudeMotorType2+"\t");
			if (myoHeadsBound>0){
				dataPW.print(instantaneousVelocity);
			}
			else{
				dataPW.print("");
			}
			dataPW.println();
			intervalCt = 0;
			totalBoundSteps += nBoundSteps;
			nBoundSteps = 0;
			miniFilamentDetached = false;
			cumSpringForces = 0;
			cumCrossbridgeForceSum=0;
			cumMyoHeadsBound = 0;
			aveDutyRatio = 0;
			cumBoundVelocity = 0;
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
		else if(tag.equalsIgnoreCase("recordActinPositions"))  {
			recordActinPositions = in.nextBoolean();
		}
		else if(tag.equalsIgnoreCase("velocityTrackingTarget") ) {
			String s = in.nextString();
			if(s.equalsIgnoreCase(MINIFILAMENT)) {
				velocityTrackingTarget = MINIFILAMENT;
			}
			else if(s.equalsIgnoreCase(FILAMENT))  {
				velocityTrackingTarget = FILAMENT;
			}
			else throw new Exception("MyoForceEvaluator.loadParameter(): got bad input for velocityTrackingTarget:  " +  velocityTrackingTarget);
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

