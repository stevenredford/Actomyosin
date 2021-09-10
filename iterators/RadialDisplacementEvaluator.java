package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import heatMap.*;

public class RadialDisplacementEvaluator extends Evaluator {
	
	Point2D uVect = new Point2D();
	Point2D oldCM = new Point2D();
	Point2D newCM = new Point2D();
	Point2D dispVec = new Point2D();
	
	PrintWriter totalDispPW;
	PrintWriter instDispPW;
	PrintWriter cumDispPW;

	
	public int radialArraySize = 25;
	public double radialVsDistInc;
	Point2D cenOfDisc;

	public double [] radialMonomerVs = new double [radialArraySize];
	public double [] cumRadialMonomerVs = new double [radialArraySize];
	public int [] radialMonomerVsCt = new int [radialArraySize];
	public double [] radialMyoVs = new double [radialArraySize];
	public double [] cumRadialMyoVs = new double [radialArraySize];
	public int [] radialMyoVsCt = new int [radialArraySize];
			
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equalsIgnoreCase("RadialArraySize")) {
			radialArraySize = in.nextInt();
		}
		else {
			super.loadParameter(tag, in);
		}
	}
	
	public void reset()  {
		super.reset();
		radialVsDistInc = 0.5*Sim2D.xDimension/(radialArraySize*Math.sqrt(2.));
		cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
		for (int i=0;i<radialArraySize;i++) {
			cumRadialMonomerVs[i] = 0;
			cumRadialMyoVs[i] = 0;
		}
	}
	
	private void resetAllMarks()  {
		for (int i=0;i<Monomer.monomerCt;i++) {
			resetMonomerMarks(Monomer.theMonomers[i]);
		}
		for (int i=0;i<Myosin.myosinCt;i++) {
			resetMyoMarks(Myosin.theMyosins[i]);
		}
	}
	
	private void resetMonomerMarks(Monomer m) {
		m.clearMarks();
		m.addMarker("radialVelocity");
		m.setMark("radialVelocity");
		m.addMarker("totalRadialDisplacement");
		m.setMark("totalRadialDisplacement");
	}
	
	private void resetMyoMarks(Myosin m) {
		m.clearMarks();
		m.addMarker("radialVelocity");
		m.setMark("radialVelocity");
		m.addMarker("totalRadialDisplacement");
		m.setMark("totalRadialDisplacement");
	}
			
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "RadialDisplacements.dat";
	}
			
	public void mkDataFile () {
		try {
			totalDispPW = new PrintWriter(new FileWriter(new File (dataPath + File.separator + "totalRadDisp.dat")),true);
			instDispPW = new PrintWriter(new FileWriter(new File (dataPath + File.separator + "InstRadDisp.dat")),true);
			cumDispPW = new PrintWriter(new FileWriter(new File (dataPath + File.separator + "cumRadDisp.dat")),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("RadialDisplacementEvaluator.mkDataFile(): error creating File"); }
	}

	public void mapOfRadialVs(String mark, boolean resetMark) {
		clearRadialVs();
		Monomer mon;
		Myosin myo;
		for (int i=0;i<Monomer.monomerCt;i++) {
			mon = Monomer.theMonomers[i];
			oldCM = mon.getMark(mark);
			mapToRadialV(mon,oldCM);
			if(resetMark) mon.setMark(mark);
		}
		for (int i=0;i<Myosin.myosinCt;i++) {
			myo = Myosin.theMyosins[i];
			oldCM = myo.getMark(mark);
			mapToRadialV(myo,oldCM);
			if(resetMark) myo.setMark(mark);
		}
		averageRadialVs();
	}
	
	public void clearRadialVs () {
		for (int i=0;i<radialArraySize;i++) {
			radialMonomerVs[i] = 0;
			radialMonomerVsCt[i] = 0;
			radialMyoVs[i] = 0;
			radialMyoVsCt[i] = 0;
		}
	}
	
	public void averageRadialVs () {
		for (int i=0;i<radialArraySize;i++) {
			if (radialMonomerVsCt[i]==0) {
				radialMonomerVs[i] = 0;
			} else {
				radialMonomerVs[i] /= radialMonomerVsCt[i];
			}
			if (radialMyoVsCt[i]==0) {
				radialMyoVs[i] = 0;
			} else {
				radialMyoVs[i] /= radialMyoVsCt[i];
			}
		}
	}
		
	public void mapToRadialV (Monomer m, Point2D oldCM) {
		newCM = m.getLocation();
		double dist = Point2D.getDistance(oldCM,cenOfDisc);
		uVect.getVector(oldCM,cenOfDisc);
		uVect.uVect();
		dispVec.getVector(oldCM, newCM);
		double dotV = Point2D.dot(dispVec, uVect);
		int rLoc = (int)(dist/radialVsDistInc);
		if(rLoc >= radialArraySize) rLoc = radialArraySize-1;
		radialMonomerVs[rLoc]+=dotV;
		radialMonomerVsCt[rLoc]++;
	}

	public void mapToRadialV (Myosin m, Point2D oldCM) {
		newCM = m.getLocation();
		double dist = Point2D.getDistance(oldCM,cenOfDisc);
		uVect.getVector(oldCM,cenOfDisc);
		uVect.uVect();
		dispVec.getVector(oldCM, newCM);
		double dotV = Point2D.dot(dispVec, uVect);
		int rLoc = (int)(dist/radialVsDistInc);
		if(rLoc >= radialArraySize) rLoc = radialArraySize-1;
		radialMyoVs[rLoc]+=dotV;
		radialMyoVsCt[rLoc]++;
	}
		
	
	/** This Evaluator writes data files. */
	public boolean hasData()  {
		return true;
	}

	/** This Evaluator writes movie frames. */
	public boolean hasMovie()  {
		return false;
	}
	
	public void doFirstEval() {
		resetAllMarks();
	}
		
	public void evaluate(double tm)  {
		writeRadialVelocityData();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}
	
	public void writeRadialVelocityData () {
		double ratio;
		mapOfRadialVs("radialVelocity",true);
		instDispPW.print("radialMonomerVs" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			instDispPW.print(radialMonomerVs[i] + "\t");
		}
		instDispPW.print("\n");
		instDispPW.print("radialMyoVs" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			instDispPW.print(radialMyoVs[i] + "\t");
		}
		instDispPW.print("\n");
		instDispPW.print("MonomerV/MyosinV" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			if(Math.abs(radialMyoVs[i]) > 0) {
				ratio = radialMonomerVs[i]/radialMyoVs[i];
			} else ratio = 0;
			instDispPW.print(ratio + "\t");
		}
		instDispPW.print("\n");
		instDispPW.flush();
		
		cumDispPW.print("radialMonomerVs" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			cumRadialMonomerVs[i] += radialMonomerVs[i];
			cumDispPW.print(cumRadialMonomerVs[i] + "\t");
		}
		cumDispPW.print("\n");
		cumDispPW.print("radialMyoVs" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			cumRadialMyoVs[i] += radialMyoVs[i];
			cumDispPW.print(cumRadialMyoVs[i] + "\t");
		}
		cumDispPW.print("\n");
		cumDispPW.print("MonomerV/MyosinV" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			if(Math.abs(cumRadialMyoVs[i]) > 0) {
				ratio = cumRadialMonomerVs[i]/cumRadialMyoVs[i];
			} else ratio = 0;
			cumDispPW.print(ratio + "\t");
		}
		cumDispPW.print("\n");
		cumDispPW.flush();
		

		mapOfRadialVs("totalRadialDisplacement",false);
		totalDispPW.print("radialMonomerVs" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			totalDispPW.print(radialMonomerVs[i] + "\t");
		}
		totalDispPW.print("\n");
		totalDispPW.print("radialMyoVs" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			totalDispPW.print(radialMyoVs[i] + "\t");
		}
		totalDispPW.print("\n");
		totalDispPW.print("MonomerV/MyosinV" + "\t");
		for(int i = 0; i < radialArraySize; i++) {
			if(Math.abs(radialMyoVs[i]) > 0) {
				ratio = radialMonomerVs[i]/radialMyoVs[i];
			} else ratio = 0;
			totalDispPW.print(ratio + "\t");
		}
		totalDispPW.print("\n");
		totalDispPW.flush();
	}
	
	public String getHeaderString() {
		return new String("\t" + name + "\t" + "totalMonomerDisplacement" + "\t" + "totalMyosinDisplacement" + "\t" + "Total ratio" + "\t" + "cumulativeMonomerDisplacement" + "\t" + "cumulativeMyosinDisplacement" + "\t" + "Cumulative ratio");
	}
	
	public String getDataString() {
		String s = new String("");
		mapOfRadialVs("totalRadialDisplacement",false);
		double avgMyoD = 0, avgMonD = 0, avgRatio = 0;
		for(int i = 0; i < radialArraySize; i++) {
			avgMonD += radialMonomerVs[i];
			avgMyoD += radialMyoVs[i];
		}
		avgMonD /= radialArraySize;
		avgMyoD /= radialArraySize;
		if(avgMyoD > 0) {
			avgRatio = avgMonD/avgMyoD;
		}
		s = s + new String("\t" + name + "\t" + avgMonD + "\t" + avgMyoD + "\t" + avgRatio);

		avgMyoD = 0; avgMonD = 0; avgRatio = 0;
		for(int i = 0; i < radialArraySize; i++) {
			avgMonD += cumRadialMonomerVs[i];
			avgMyoD += cumRadialMyoVs[i];
		}
		avgMonD /= radialArraySize;
		avgMyoD /= radialArraySize;
		if(avgMyoD > 0) {
			avgRatio = avgMonD/avgMyoD;
		}
		s = s + new String("\t" + avgMonD + "\t" + avgMyoD + "\t" + avgRatio);
		return s;
	}
	
	public void cleanUp()  {}
				

}
