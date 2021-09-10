package iterators;


import java.util.*;
import java.io.*;

import analysis.ValueTracker;

import main.*;
import io.*;
import util.*;

public class MyoDissEvaluator extends Evaluator {

	/** A histogram that keeps running track of distribution of minifilament bound lifetimes.. */
	HistogramPlus lifetimeHistogram;

	boolean fastFilamentRelaxation = false;

	double firstAttchTime;

	double firstDetachTime;

	double startEvalTime = 0;

	int numDetachEvents;

	boolean miniFilamentDetached = true;

	Vector attachLifetimes, detachLifetimes;
	
	PrintWriter	dataPW2;
	
	String	dataFileName2;
	


	double timeToFirstLongAttachment = Sim2D.simulationTime;

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
		miniFilamentDetached = true;
		totalMotor1ReleaseEvents = 0;
		totalMotor2ReleaseEvents = 0;


	}

	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = new String(dataPath + File.separator + "MyoLifetimes.dat");
		dataFileName2 = new String(dataPath + File.separator + "ATPaseRates.dat");

	}
	public String getHeaderString() {
		return new String("numDetachmentEvents\t" + "mean attached lifetime\t" + "first last attachment time\t" + "time to first long attachment\t" + "mean detached lifetime\t");
	}



	public String getDataString() {
		double meanL = 0;
		int cnt = 0;
		Double l;
		if(!miniFilamentDetached) {
			attachLifetimes.add(new Double(Sim2D.simulationTime-firstAttchTime));
		}
		for(Enumeration e = attachLifetimes.elements(); e.hasMoreElements(); ) {
			l = ((Double)e.nextElement()).doubleValue();
			dataPW.print(l + "\t");
			meanL += l;
			cnt++;
		}
		dataPW.println();
		if(cnt>0) meanL/=cnt;
		String s = new String(numDetachEvents + "\t");
		if(numDetachEvents>0)
			s = s + meanL + "\t";
		else s = s + Sim2D.simulationTime + "\t";

		s = s + firstAttchTime + "\t" + timeToFirstLongAttachment + "\t";

		meanL = 0; cnt = 0;
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
		else s = s + "0" + "\t";
		return s;
	}
	public void writeDataFileHeader() {
		dataPW.println("ATPaseRateMotor1"+"\t"+"ATPaseRateMotor2");
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			dataPW2 = new PrintWriter(new FileWriter(new File (dataFileName2)),true);
		} catch (IOException ioe) { System.out.println("MyoDissEvaluator.mkDataFile(): error creating File"); }
	}


	public void evaluate(double tm)  {
		if(tm < startEvalTime) return;
		if(MyosinMiniFilament.theMiniFilaments[0].getNumBoundMyosins() == 0) {
			if(!miniFilamentDetached) {
				firstDetachTime = tm;
				double attachTm = tm-firstAttchTime;
				if(attachTm > 30) timeToFirstLongAttachment = firstAttchTime;
				attachLifetimes.add(new Double(attachTm));
				numDetachEvents++;
				setMinifilamentDetached(true);
			}
		}
		else {  // bound by at least one head.
			if(miniFilamentDetached) {
				firstAttchTime = tm;
				if(firstDetachTime > 0) detachLifetimes.add(new Double(tm-firstDetachTime));
				setMinifilamentDetached(false);
			}
		}
		getNumMotorReleaseEvents();
		double motor1ATPaseRate=totalMotor1ReleaseEvents/interval;
		totalMotor1ReleaseEvents=0;
		double motor2ATPaseRate=totalMotor2ReleaseEvents/interval;
		totalMotor2ReleaseEvents=0;
		dataPW2.println(motor1ATPaseRate+"\t"+motor2ATPaseRate);
	}
	double totalMotor1ReleaseEvents = 0;
	double totalMotor2ReleaseEvents = 0;

	public void getNumMotorReleaseEvents(){
		
		for (int i=0;i<Myosin.myosinCt; i++){
			Myosin m = Myosin.theMyosins[i];
			if (!m.iAmAnAMyo){
				if (m.usedReleaseRateMult==Myosin.adpReleaseRateMult) {
					totalMotor1ReleaseEvents+=m.getNumADPReleaseEvents();
				}
				else if (m.usedReleaseRateMult==Myosin.secondADPReleaseRateMult){
					totalMotor2ReleaseEvents+=m.getNumADPReleaseEvents();
				}
				m.zeroNumADPReleaseEvents();
			}
		}

	}

	public void setMinifilamentDetached(boolean d) {
		if(miniFilamentDetached == d) return;
		miniFilamentDetached = d;
		if(d) {
			if(MyosinMiniFilament.theMiniFilaments[0].pinned && fastFilamentRelaxation) {
				Actin.theActins[0].restoreAttachmentPositions();
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

