package collision;

/**
 * MyoCollisionDetector.java
 *
 * @author Created by Omnicore CodeGuide
 */

import util.*;
import main.*;

public class MyoCollisionDetector extends CollisionDetector
{
	
//	double nextMyoBindCheck;
//
//	double maxMyoBindProb;
//
//	double avMaxMyoBindProb;
//
//	double avMaxMove;
//
//	double myoBindCheckInterval;

	/** Diagnostic: Number of times a myosin checked an actin filament for collision */
	public int myoCks = 0;
	
	/** Diagnostic: Number of times a myosin bound an actin filament */
	public int myoBounds = 0;

	/** Used in binder/actin collision checking to ensure each filament checked only once.*/
	boolean[] actinChecked = new boolean[5000];
	
	/** Keeps track of which actins were checked in a certain bin so we can reset actinChecked efficiently. */
	int[] checkedActins = new int[5000];
	
	int checkedActinCt;
	
	BindingInfo info = new BindingInfo();
	
	/** A list of possible filament binding partners computed for each Myosin head. */
	private Actin [] possibles = new Actin[1000];
	
	/** A list of target binding positions on candidate filaments expressed as monomer position. */
	private int []  targetPositions = new int[1000];

	/** A list of cumulative binding probabilities across list of all possible binding partners. */
	private double [] cumCandidateProb = new double[1000];

	/** Number of possible binding partners for the current Myo. */
	private int numPossibles;
	
	/** This defines the number of monomers to either side of tyhe preferred binding position
	 * a myosin head can bind.  This is determined by the span of the probability distribution
	 * for myosin head binding which is statically computed and stored by the Myosin class.
	 *
	 **/
	private int headReach;
	
	private double [] cumP;
	
	public void checkCollisions(double simTime) {
		checkMyosinMesh();
	}

	/** run through myosin mesh, checking each grid element for registered myosins.  For each, consider possible binding partners. */
	public void checkMyosinMesh(){
		cumP = Myosin.cumBindingProbabilities;
		headReach = cumP.length/2;
		Myosin myosin;
		for(int i = 0; i < Myosin.myosinCt; i++) {
			myosin = Myosin.theMyosins[i];
			numPossibles = 0;
			if(myosin.canBind()){
				sampleActinMeshTargetsForMyosin(myosin);
			}
			checkFilamentCandidates(myosin);
			myosin.numPossibles=numPossibles;
		}
	}
	
	private void addPossible(BindingInfo info) {
		double pTot = numPossibles != 0 ? cumCandidateProb[numPossibles-1] : 0;
		int offset = (int)Math.round(info.arcL/Actin.monLength);
		double p;
					
		if(offset >= -headReach && offset <= info.actin.nMonomers-1+headReach) {
			if(offset < headReach) { // preferred binding point near pEnd
				p = cumP[offset+headReach];
			}
			else {
				p = 1.0;
			}
			if(offset > info.actin.nMonomers-1-headReach) { // preferred binding point near bEnd
				p -= cumP[offset - (info.actin.nMonomers-headReach)];
			}
			possibles[numPossibles] = info.actin;
			targetPositions[numPossibles] = offset;
			pTot += Sim2D.getDeltaT()*Myosin.myoUniformBindingProb*p;
			cumCandidateProb[numPossibles] = pTot;
			numPossibles++;
		}
		
	}
		
	
	private void sampleActinMeshTargetsForMyosin(Myosin myo) {

		Actin actin;
		checkedActinCt = 0;
		int x, y;
		Point2D intPt = new Point2D();

		Point2D center = myo.attachPt;

		int startBinX = Mesh.getBinX(center.x - Myosin.minBindingDistance);
		int stopBinX = Mesh.getBinX(center.x + Myosin.minBindingDistance);
		int startBinY = Mesh.getBinY(center.y - Myosin.minBindingDistance);
		int stopBinY = Mesh.getBinY(center.y + Myosin.minBindingDistance);

		// Sample all intersected bins in ACTIN_MESH for candidate filaments,
		// making sure every filament is checked only once...
		for(int xx=startBinX;xx<=stopBinX;xx++){
			for(int yy=startBinY;yy<=stopBinY;yy++){

				x = Mesh.wrapX(xx); y = Mesh.wrapY(yy);

				if(Mesh.ACTIN_MESH.timeStamps[x][y]==Sim2D.counter) {

					for(int i=0;i<Mesh.ACTIN_MESH.activeCts[x][y];i++){
						int actinId = (int)Mesh.ACTIN_MESH.meshpoints[x][y][i];
						if (!actinChecked[actinId]) {
							actinChecked[actinId]=true;
							checkedActins[checkedActinCt]=actinId;
							checkedActinCt++;
							actin=Actin.theActins[actinId];
							myoCks++;
//System.out.println("tm = " + Sim2D.simulationTime + " myoChks = " + myoCks  + " actinID = " + actinId);

							if (actinCollision(actin,intPt,myo.attachPt,Myosin.minBindingDistance, info,headReach*Actin.monLength)) {
								if(ckActinVsMyoAngle(actin,myo)) {
									getMyoBindingSite(actin,info);
									addPossible(info);
								}
							}
						}
					}
				}
			}
		}
		// reset the actins that were checked for this myo
		for(int j=0;j<checkedActinCt;j++) {
			actinChecked[checkedActins[j]]=false;
		}
	}
	
//	private void sampleActinMeshTargetsForMyosin(Myosin myo) {
//
//		Actin actin;
//		Point2D intPt = new Point2D();
//
//		for(int i = 0; i < Actin.actinCt; i++) {
//			actin = Actin.theActins[i];
//			if (actinCollision(actin,intPt,myo.attachPt,Myosin.minBindingDistance, info,headReach*Actin.monLength)) {
//				if(ckActinVsMyoAngle(actin,myo)) {
//					getMyoBindingSite(actin,info);
//					addPossible(info);
//				}
//			}
//		}
//		// reset the actins that were checked for this myo
//		for(int j=0;j<checkedActinCt;j++) {
//			actinChecked[checkedActins[j]]=false;
//		}
//	}

	private void checkFilamentCandidates(Myosin myo) {
		if(numPossibles>0) {
			myo.closeEnoughToBind = true;
			double prob = Math.random();
			double residualProb = prob;
			if (prob < cumCandidateProb[numPossibles-1]) { // myo head binds somewhere.
				int choice = 0; // randomly select which filament
				while(prob > cumCandidateProb[choice]) {
					residualProb = prob - cumCandidateProb[choice];
					choice++;
				}
				checkFilamentCandidate(myo,choice,residualProb); // now ask which monomer along this filament
			}
		} else {
			myo.closeEnoughToBind = false;
		}
	}
	
	private void checkFilamentCandidate(Myosin myo, int fil, double residualP) {

		double arcL;
		
		int choice = 0;
		int offset = 0;
		
		if(targetPositions[fil] < headReach) {
			offset = headReach-targetPositions[fil];
			residualP += cumP[offset-1];
		}
		
		while(residualP > cumP[choice+offset]) {
			choice++;
		}
			
		arcL = (targetPositions[fil] - headReach + choice + offset)*Actin.monLength;
		Monomer m = possibles[fil].getMonomerAt(arcL);
		if(m.isFree()) {
			myo.bindMon(m);
			myoBounds++;
		}
	}
			
				
	/** Given the position of closest point on actin filament in arc length and the distance to that closest point
	 * from this myosin head (stored in info), find the preferred binding position along */
	private void getMyoBindingSite(Actin actin, BindingInfo info) {
		info.arcL += Math.sqrt(Myosin.minBindingDistance*Myosin.minBindingDistance - info.minD*info.minD);
		info.actin = actin;
	}

	/** Check whether the myosin  has the right angle with respect to the actin to allow binding. */
	public boolean ckActinVsMyoAngle(Actin actin, Myosin m) {
		if (m.myHolder instanceof MyosinSurface) return true; // no angle constraint for myosins attached to a surface
		double cosAngle = Point2D.dot(m.getAim(),actin.uVect);
		if(cosAngle > Myosin.myoAngleBias) {
			return true;
		}
		else return false;
	}

}


