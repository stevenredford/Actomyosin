/**
 * InternalStressTracker.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;



import java.text.DecimalFormat;
import java.awt.*;
import java.awt.geom.AffineTransform;
import util.*;


/** A utitlity class that tracks local tension and compression along a filament over time for implementing
 * e.g. buckling calculations etc.
 */

public class InternalStressTracker
{
	
	static double buckleFactor;
	
	/** Max number of force points along filament. */
	int maxAxialForces = 500;
	
	/** Current number of force points along filament. */
	int axialForceCt;
	
	double maxBuckleProgress;
	
	Monomer maxBuckleMon = null;
	
	/** stores force for each force applied to an actin filament */
	public double [] axialForces;
	
	/** stores point of application of axial forces. */
	public Point2D [] axialForcePts;
	
	/** store the segment length between each force pt */
	public double [] forceSegLength;
	
	
	public double [] buckleRates;
	
	/** positive number is tension, negative compression */
	public double [] curInternalForces;

	/** positive number is tension, negative compression */
	public double [] aveInternalForces;
	
	/** Number of monomers in each stress segment. */
	public int [] numMonsInSeg;
	
	/** Fraction of critical buckling force on each segment. */
	public double [] buckleFraction;
	
	static DecimalFormat buckleFracFormat = new DecimalFormat ("#0.0#; #0.0#");

	
	public void initForceTracking()  {
		axialForceCt = 0;
		axialForces = new double [maxAxialForces];
		axialForcePts = new Point2D [maxAxialForces];
		forceSegLength = new double [maxAxialForces];
		buckleRates = new double [maxAxialForces];
		curInternalForces = new double [maxAxialForces];
		aveInternalForces = new double [maxAxialForces];
		numMonsInSeg = new int [maxAxialForces];
		buckleFraction = new double [maxAxialForces];
	}

	public void trackStresses(Actin a) {
		findAxialForces(a);
		findInternalForces(a);
		getForceSegLengths();
		findBucklingFraction(a);
		mapStressToMonomers(a);
		//aveInternalForces(a);
	}

		// Load axial force components and points of application
	public void findAxialForces (Actin a) {
		// load the pin force at pointed-end... always index 0, even if no force
		if (a.isPPinned()) {
			axialForces[0] = Point2D.dot(a.pAttachment.getPinForce(), a.uVect);
		} else {
			axialForces[0] = 0;
		}
		axialForcePts[0] = a.pEnd;
		axialForceCt = 1;
		
		// load all the forces
		Monomer curM = a.pEndMonomer;
		boolean atEnd = false;
		while(!atEnd && curM!=null){
			curM.monStressSeg = axialForceCt-1;
			if (!curM.isFree()) {
				axialForces[axialForceCt] = Point2D.dot(curM.myForce,a.uVect);
				axialForcePts[axialForceCt] = curM.getLocation();
				axialForceCt++;
			}
			
			if (curM == a.bEndMonomer) {
				atEnd = true;
			} else {
				curM = curM.next;
			}
		}
		
		// load the pin force at barbed-end... always a force added here, even if zero
		if (a.isBPinned()) {
			axialForces[axialForceCt] = Point2D.dot(a.bAttachment.getPinForce(), a.uVect);
		} else {
			axialForces[axialForceCt] = 0;
		}
		axialForcePts[axialForceCt] = a.bEnd;
		axialForceCt++;
	}

	public void findInternalForces (Actin a) {
		double bPointForceSum = 0;	// force on barbed-end point in internal segment
		double pPointForceSum = 0;	// force on pointed-end point in internal segment

		for (int i=0;i<axialForceCt;i++) { bPointForceSum += axialForces[i]; }  // sum all axial forces
		for (int i=0;i<axialForceCt;i++) {

			bPointForceSum -= axialForces[i];
			pPointForceSum += axialForces[i];
			
			if (pPointForceSum > 0 && bPointForceSum < 0) {  // this is the only condition for internal compression
				if (pPointForceSum > Math.abs(bPointForceSum)) {
					curInternalForces[i] = bPointForceSum;
				} else {
					curInternalForces[i] = -pPointForceSum;
				}
			} else if (pPointForceSum < 0 && bPointForceSum > 0) { // condition for internal tension
				if (Math.abs(pPointForceSum) > bPointForceSum) {
					curInternalForces[i] = bPointForceSum;
				} else {
					curInternalForces[i] = -pPointForceSum;
				}
			} else {
				curInternalForces[i] = 0;  // shouldn't ever set value here... just in case though
			}
		}


		double buckleForce = buckleFactor/Math.pow(a.physicalLength, 2);
		double buckleFor300nm = buckleFactor/Math.pow(300, 2);

	}
	
	public void findBucklingFraction (Actin a) {
		double buckleForSeg;
		for (int i=0;i<axialForceCt-1;i++) {
			if (forceSegLength[i] == 0) {
				buckleFraction[i] = 0;
				buckleRates[i] = 0;
			} else {
				buckleForSeg = buckleFactor/Math.pow(forceSegLength[i],2);
				buckleFraction[i] = -curInternalForces[i]/buckleForSeg;
				buckleRates[i] = getBuckleRates(a,forceSegLength[i]);
			}
		}
	}
	
	public void severFromBuckling (Actin a) {
		if(maxBuckleProgress >= 1) {
			int greatestStressSeg = maxBuckleMon.monStressSeg;
			double segLength = Point2D.getDistance(axialForcePts[greatestStressSeg],axialForcePts[greatestStressSeg+1]);
			double arcToSever = Point2D.getDistance(a.pEnd, axialForcePts[greatestStressSeg]) + 0.5*segLength;
			Monomer severMon = a.getMonomerAt(arcToSever);
			a.addSeverCandidate(severMon);
		}
	}

	public double getBuckleRates (Actin a, double segLength) {
		double heightAbove = 10;  // nm  assumed height of filament above surface (see Howard pg 107)
		double cPerp = 4*Math.PI*Constants.viscosity/(Math.log(2*heightAbove/a.radius));
		double timeC = (cPerp/a.actinEI)*Math.pow(segLength/(1.5*Math.PI), 4);
		return Sim2D.getDeltaT()/timeC;
	}

	public void mapStressToMonomers (Actin a) {
		Monomer curM = a.pEndMonomer;
		boolean atEnd = false;
		maxBuckleProgress = 0;
		double tmp = 0;
		while(!atEnd && curM!=null){
			tmp = curM.registerStress(curInternalForces[curM.monStressSeg],buckleFraction[curM.monStressSeg],buckleRates[curM.monStressSeg]);
			if(tmp > maxBuckleProgress) {
				maxBuckleProgress = tmp;
				maxBuckleMon = curM;
			}
			if (curM == a.bEndMonomer) {
				atEnd = true;
			} else {
				curM = curM.next;
			}
		}
	}
	
	public void aveInternalForces (Actin a) {
		Monomer curM = a.pEndMonomer;
		boolean atEnd = false;
		
		// zero stresses and counts
		for (int i=0; i<axialForceCt-1; i++) {
			aveInternalForces[i] = 0;
			numMonsInSeg[i] = 0;
		}
		
		while(!atEnd && curM!=null){
			aveInternalForces[curM.monStressSeg] += curM.avgMonStress;
			numMonsInSeg[curM.monStressSeg]++;
			
			if (curM == a.bEndMonomer) {
				atEnd = true;
			} else {
				curM = curM.next;
			}
		}
		
		// average values
		for (int i=0; i<axialForceCt-1; i++) {
			if (numMonsInSeg[i] == 0) {
				aveInternalForces[i] = 0;
			} else {
				aveInternalForces[i] = aveInternalForces[i]/numMonsInSeg[i];
			}
		}
	}
	
	public void getForceSegLengths() {
		for (int i=0;i<axialForceCt-1;i++) {
			forceSegLength[i] = Point2D.getDistance(axialForcePts[i],axialForcePts[i+1]);
		}
	}
	public void showInternalStress (Actin a,Graphics g, double scale) {
		int greatestStressSeg = 0;
		for (int i=0;i<axialForceCt-1;i++) {
			if (buckleFraction[i] > buckleFraction[greatestStressSeg]) { greatestStressSeg = i; }
		}
		if (buckleFraction[greatestStressSeg] > 1) {
			g.setColor(Color.yellow);
			g.fillOval(axialForcePts[greatestStressSeg].getPixX(),axialForcePts[greatestStressSeg].getPixY(),10,10);
			g.fillOval(axialForcePts[greatestStressSeg+1].getPixX(),axialForcePts[greatestStressSeg+1].getPixY(),10,10);
			g.setColor(Color.white);
			Point2D [][] sgs = Point2D.segmentLine(axialForcePts[greatestStressSeg],axialForcePts[greatestStressSeg+1],a.uVect);
			for(int i = 0; i < sgs.length; i++) {
				g.drawLine(sgs[i][0].getPixX(), sgs[i][0].getPixY(), sgs[i][1].getPixX(), sgs[i][1].getPixY());
			}
			printBuckleFraction (axialForcePts[greatestStressSeg],axialForcePts[greatestStressSeg+1],g,greatestStressSeg);
		}
	}
	public void printBuckleFraction (Point2D p1, Point2D p2, Graphics g, int ind) {
		//System.out.println ("BuckleFraction for " + ind + " = " + buckleFraction[ind]);
		if (buckleFraction[ind] > 0.1) {
			int xOff = -10;
			int yOff = -3;
			int xMid = (p1.getPixX() + p2.getPixX())/2 + xOff;
			int yMid = (p1.getPixY() + p2.getPixY())/2 + yOff;
			g.setColor(Color.white);
			g.setFont(new Font ("TimesRoman", Font.PLAIN, 14));
			g.drawString(buckleFracFormat.format(buckleFraction[ind]), xMid, yMid);
			g.setFont(new Font ("TimesRoman", Font.PLAIN, 10));
		}
	}
	
	public void printBuckleFractionVertical (Point2D p1, Point2D p2, Graphics g, int ind) {
		Graphics2D g2d = (Graphics2D)g;
		int yOff = -10;
		int xMid = (p1.getPixX() + p2.getPixX())/2;
		int yMid = (p1.getPixY() + p2.getPixY())/2 + yOff;
	    // clockwise 90 degrees
	    AffineTransform at = new AffineTransform();
	    // thanks to M.C. Henle for the bug fix!
	    at.setToRotation(-Math.PI/2.0);
	    g2d.setTransform(at);
	    g2d.drawString("Vertical text", xMid, yMid);
		    
	}

}

