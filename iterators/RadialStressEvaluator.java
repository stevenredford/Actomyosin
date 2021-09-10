package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import heatMap.*;
import parameters.*;

public class RadialStressEvaluator extends Evaluator {
	

	Point2D rVec = new Point2D();
	Point2D cm = new Point2D();
	double dist;
	double stress;
	double rDot;
	double cDot;
	
	public int radialArraySize = 25;
	public double radialInc = 0.5*Sim2D.xDimension/(radialArraySize*Math.sqrt(2.));
	Point2D cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
	
	public double [] radialActinStress = new double [radialArraySize];
	public double [] circumferentialActinStress = new double [radialArraySize];
	public double [] radialMyoStress = new double [radialArraySize];
	public double [] circumferentialMyoStress = new double [radialArraySize];
	public double [] radialLinkerStress = new double [radialArraySize];
	public double [] circumferentialLinkerStress = new double [radialArraySize];
	double totalRadActinStress = 0, totalRadMyosinStress = 0, totalRadLinkerStress = 0, totalRadialStress = 0;
	double totalCircActinStress = 0, totalCircMyosinStress = 0, totalCircLinkerStress = 0, totalCircumferentialStress = 0;
			
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
	}
	
	public void init(String path, String name) throws Exception  {
		super.init(path,name);
		Parameters.setParameter("monitorStress",true);
	}
	
	public void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "RadialStress.dat";
	}
	
	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("RadialStressEvaluator.mkDataFile(): error creating File"); }
	}

	public void clearRadialStress () {
		for (int i=0;i<radialArraySize;i++) {
			radialActinStress[i] = circumferentialActinStress[i] = 0;
			radialMyoStress[i] = circumferentialMyoStress[i] = 0;
			radialLinkerStress[i] = circumferentialLinkerStress[i] = 0;
		}
	}
		
	public void mapToRadialStress (Monomer m) {
		cm = m.getLocation();
		dist = Point2D.getDistance(cm,cenOfDisc);
		int rLoc = (int)(dist/radialInc);
		if(rLoc >= radialArraySize) rLoc = radialArraySize-1;
		rVec.getVector(cm,cenOfDisc);
		rVec.uVect();
		rDot = Math.abs(Point2D.dot(m.getFilament().uVect, rVec));
		cDot = Math.sqrt(1-rDot*rDot);
		radialActinStress[rLoc]+=m.avgMonStress*rDot*rDot*Actin.monLength;
		circumferentialActinStress[rLoc]+=m.avgMonStress*cDot*rDot*Actin.monLength;
	}

	/** Maps radial and circumferential stress on each stress segment of this minifilament. */
	public void mapToRadialStress (MyosinMiniFilament m) {
		Point2D [] centers = m.getStressCenters();
		int centerIndex = MyosinMiniFilament.nMyosinHeads/2-1;
		Point2D mVect = m.uVect;
		double len;
		for(int i = 0; i < centers.length; i++) {
			dist = Point2D.getDistance(centers[i],cenOfDisc);
			int rLoc = (int)(dist/radialInc);
			if(rLoc >= radialArraySize) rLoc = radialArraySize-1;
			rVec.getVector(centers[i],cenOfDisc);
			rVec.uVect();
			rDot = Math.abs(Point2D.dot(mVect, rVec));
			cDot = Math.sqrt(1-rDot*rDot);
			if(i == centerIndex) {
				len = m.length;
			} else len = MyosinMiniFilament.myoSpacing;
			radialMyoStress[rLoc]+=m.getInternalStress(i)*rDot*rDot*len;
			circumferentialMyoStress[rLoc]+=m.getInternalStress(i)*cDot*rDot*len;
		}
	}

	/** Maps radial and circumferential stress of this minifilament. */
	public void mapToRadialStress (Crosslinker  l) {
		if(l.isFree()) return;
		cm = l.getCentroid();
		dist = Point2D.getDistance(cm,cenOfDisc);
		int rLoc = (int)(dist/radialInc);
		if(rLoc >= radialArraySize) rLoc = radialArraySize-1;
		rVec.getVector(cm,cenOfDisc);
		rVec.uVect();
		rDot = Math.abs(Point2D.dot(l.uVect, rVec));
		cDot = Math.sqrt(1-rDot*rDot);
		radialLinkerStress[rLoc]+=l.forceAv*rDot*rDot*l.getLength();
		circumferentialLinkerStress[rLoc]+=l.forceAv*cDot*rDot*l.getLength();
	}

	public void mapOfRadialStress() {
		clearRadialStress();
		for (int i=0;i<Monomer.monomerCt;i++) {
			mapToRadialStress(Monomer.theMonomers[i]);
		}
		for (int i=0;i<MyosinMiniFilament.miniFilamentCt;i++) {
			mapToRadialStress(MyosinMiniFilament.theMiniFilaments[i]);
		}
		for (int i=0;i<Crosslinker.crosslinkerCt;i++) {
			mapToRadialStress(Crosslinker.theCrosslinkers[i]);
		}
		double r = radialInc/2;
		double circ;
		double den;
		for(int i = 0; i < radialArraySize; i++) {
			circ = 2*Math.PI*r;
			den = circ*radialInc;
			totalRadActinStress += radialActinStress[i];
			totalRadMyosinStress += radialMyoStress[i];
			totalRadLinkerStress += radialLinkerStress[i];
			totalCircActinStress += circumferentialActinStress[i];
			totalCircMyosinStress += circumferentialMyoStress[i];
			totalCircLinkerStress += circumferentialLinkerStress[i];
			radialActinStress[i]/=den;
			radialMyoStress[i]/=den;
			radialLinkerStress[i]/=den;
			circumferentialActinStress[i]/=den;
			circumferentialMyoStress[i]/=den;
			circumferentialLinkerStress[i]/=den;
			r += radialInc;
		}
	}
	
	
	/** This Evaluator writes data files. */
	public boolean hasData()  {
		return true;
	}

	public String getHeaderString() {
		return new String("\t" + name + "\t" + "totalRadActinStress" + "\t" + "totalRadMyosinStress" + "\t" + "totalRadLinkerStress" + "\t" + "totalRadialStress" +
							  "\t" + "totalCircActinStress" + "\t" + "totalCircMyosinStress" + "\t" + "totalCircLinkerStress" + "\t" + "totalCircumferentialStress");
	}
	
	public String getDataString() {
		double den = Math.PI*radialArraySize*radialArraySize*evaluationCt;
		totalRadActinStress/=den;
		totalRadMyosinStress/=den;
		totalRadLinkerStress/=den;
		totalCircActinStress/=den;
		totalCircMyosinStress/=den;
		totalCircLinkerStress/=den;
		totalRadialStress += totalRadActinStress + totalRadMyosinStress + totalRadLinkerStress;
		totalCircumferentialStress += totalCircActinStress + totalCircMyosinStress + totalCircLinkerStress;
		return new String("\t" + name + "\t" + totalRadActinStress + "\t" + totalRadMyosinStress + "\t" + totalRadLinkerStress + "\t" + totalRadialStress +
							  "\t" + totalCircActinStress + "\t" + totalCircMyosinStress + "\t" + totalCircLinkerStress + "\t" + totalCircumferentialStress);
	}
	/** This Evaluator writes movie frames? */
	public boolean hasMovie()  {
		return false;
	}
	
	public void evaluate(double tm)  {
		writeRadialStressData();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}
	
	public void writeRadialStressData () {
		mapOfRadialStress();
		dataPW.print("radialActinStress\t");
		for(int i = 0; i < radialActinStress.length; i++) {
			dataPW.print(radialActinStress[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.print("radialMyoStress\t");
		for(int i = 0; i < radialMyoStress.length; i++) {
			dataPW.print(radialMyoStress[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.print("radialLinkerStress\t");
		for(int i = 0; i < radialLinkerStress.length; i++) {
			dataPW.print(radialLinkerStress[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.print("circumferentialActinStress\t");
		for(int i = 0; i < circumferentialActinStress.length; i++) {
			dataPW.print(circumferentialActinStress[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.print("circumferentialMyoStress\t");
		for(int i = 0; i < circumferentialMyoStress.length; i++) {
			dataPW.print(circumferentialMyoStress[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.print("circumferentialLinkerStress\t");
		for(int i = 0; i < circumferentialLinkerStress.length; i++) {
			dataPW.print(circumferentialLinkerStress[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.flush();
	}
	
	public void cleanUp()  {}
				

}
