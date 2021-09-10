package collision;

/**
 * XLCollisionDetector.java
 *
 * @author Created by Omnicore CodeGuide
 */

import main.*;
import util.*;


public class XLCollisionDetector extends CollisionDetector
{
	
	double nextXLBindCheck;
	
	double maxXLBindProb;
	
	double avMaxXLBindProb;
	
	double avMaxMove;
	
	double xlBindCheckInterval;
	
	/** A list of possible Monomer binding partners computed for each XL. */
	private Monomer [] xlpossibles = new Monomer[1000];
	
	/** arclength position for each monomer candidate. */
	private double [] arcLengths = new double[1000];
	
	/** min distances for each monomer candidate. */
	private double [] distances = new double[1000];
	
	/** A list of total binding probability for each possible binding partner. */
	private double [] possibleFilProb = new double[1000];
	
	/** Number of possible binding partners for the current XL. */
	private int numPossibles;
	
	BindingInfo info = new BindingInfo();
	
	/** Used in binder/actin collision checking to ensure each filament checked only once.*/
	boolean[] actinChecked = new boolean[5000];
	
	/** Keeps track of which actins were checked in a certain bin so we can reset actinChecked efficiently. */
	int[] checkedActins = new int[5000];
	
	int checkedActinCt;
	
	
	public void init() {
		nextXLBindCheck = xlBindCheckInterval = Sim2D.getDeltaT();
		avMaxXLBindProb = avMaxMove = 0;
	}
		

	public void checkCollisions(double simTime) {
		if(simTime > nextXLBindCheck) {
			double maxMove = Actin.getMaxFilamentMoveSinceLastCheck();
			avMaxMove = 0.9*avMaxMove + 0.1*maxMove;
			maxXLBindProb = checkCrosslinkerMesh(xlBindCheckInterval);
			avMaxXLBindProb = 0.9*avMaxXLBindProb + 0.1*maxXLBindProb;
			if(avMaxXLBindProb < Crosslinker.xlTargetBindingProbPerTimestep &&
			  avMaxMove < Crosslinker.targetMonomerMoveFraction*Crosslinker.reach && xlBindCheckInterval < Crosslinker.xlBindCheckMaxInterval) {
				xlBindCheckInterval*=1.1;
				avMaxXLBindProb*=1.1;
				avMaxMove*=1.1;
			}
			if(avMaxXLBindProb > Crosslinker.xlMaxBindingProbPerTimestep ||
			  avMaxMove > Crosslinker.maxMonomerMoveFraction*Crosslinker.reach) {
				xlBindCheckInterval*=0.9;
				avMaxXLBindProb*=0.9;
				avMaxMove*=0.9;
			}
			nextXLBindCheck += xlBindCheckInterval;
		}
	}
	
	private void addPossible(BindingInfo i) {
		xlpossibles[numPossibles] = i.monomer;
		distances[numPossibles] = i.minD;
		arcLengths[numPossibles] = i.arcL;
		numPossibles++;
	}
		
	
	/** run through crosslinker mesh, checking each grid element for registered linkers.  For each, consider possible binding partners. */
	public double checkCrosslinkerMesh(double interval) {
		Crosslinker linker;
		XLVariant linkerProps; 
		double [][] pTable;//=new double[10] [10];  //changed here
		double bindingProb, maxBindingProb = 0;
		for(int i = 0; i < Crosslinker.crosslinkerCt; i++) {
			linker=Crosslinker.theCrosslinkers[i];
			numPossibles = 0;
			if(linker.isFree()) {
				linkerProps = linker.getXLVariant();
				linkerProps.getSpringConstant();
				pTable = linkerProps.getDistanceProbabilityTable();
				linker.updateBindingSites();
				bindingProb = sampleActinMeshTargetsForLinker(linker, interval,linkerProps,pTable);
				if(bindingProb > maxBindingProb) {
					maxBindingProb = bindingProb;
				}
				checkFilamentCandidates(linker, linkerProps,pTable);
			}
		}
		return maxBindingProb;
	}
	
	private double sampleActinMeshTargetsForLinker(Crosslinker linker, double interval, XLVariant linkerProps,double [][] pTable) {

		Actin actin;
		checkedActinCt = 0;
		
		int x, y;

		Point2D center = linker.getSinglyBoundSite();
		
		int startBinX = Mesh.getBinX(center.x - linker.reach);
		int stopBinX = Mesh.getBinX(center.x + linker.reach);
		int startBinY = Mesh.getBinY(center.y - linker.reach);
		int stopBinY = Mesh.getBinY(center.y + linker.reach);
		
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
							ckActinVsXL(actin,linker);
						}
					}
				}
			}
		}
		
		// reset the actins that were checked for this bin
		for(int j=0;j<checkedActinCt;j++) {
			actinChecked[checkedActins[j]]=false;
		}
		
		// make a cumulative probability array for the filaments that were found to be close enough
		double pTot=0;
	
		for (int i=0; i<numPossibles; i++){
			pTot+=pTable[(int)(distances[i]/linkerProps.getDistInterval())][2]*xlBindCheckInterval;
			possibleFilProb[i]=pTot;
		}
		return numPossibles > 0 ? possibleFilProb[numPossibles-1] : 0.0;
	}
		
	/** Check whether the actin is in range of the crosslinker and if so add the closest monomer on that filament to the possibles list. */
	public void ckActinVsXL (Actin actin, Crosslinker linker) {
		if(linker.isBoundTo(actin)){ return; }
		Point2D intPt = new Point2D();
		BindingInfo info = new BindingInfo();
		if(linker.site1Bound()) {
			if(ckActinVsXLAngle(actin,linker,linker.bound1.getFilament())) {
				if(!actinCollision(actin,intPt,linker.site1, linker.reach,info,0)) {
					return;
				}
			}
		}
		else if(linker.site2Bound()) {
			if(ckActinVsXLAngle(actin,linker,linker.bound2.getFilament())) {
				if(!actinCollision(actin,intPt,linker.site2, linker.reach,info,0)) {
					return;
				}
			}
		}
			
		if(info.arcL < 0 || info.arcL > actin.length) {
			return;
		}
		info.monomer=actin.getMonomerAt(info.arcL);
		if (info.monomer == null) { return ; }
		addPossible(info);
		return;
	}

	
	
	
	// Test random number against total Prob to see if binding occurs.
	// If so, randomly select candidate filament and then which monomer on that filament.
	private void checkFilamentCandidates(Crosslinker linker, XLVariant linkerProps, double [][] pTable){
		if(numPossibles>0) {
			double Prob = Math.random();
			if (Prob < possibleFilProb[numPossibles-1]){
				int jj = 0;
				while(Prob > possibleFilProb[jj]) {
					jj++;
				}
				checkFilamentCandidate(linker,jj,linkerProps,pTable);
			}
		}
	}
	
	//draw random monomer on the filament and bind it
	private void checkFilamentCandidate(Crosslinker linker,int fil, XLVariant linkerProps, double [][] pTable) {
		Monomer monCandidate = null;
		boolean foundMon = false;
		int counter = 0;
		double A, Prob;
		int binChoice = (int)(Math.round(distances[fil]/linkerProps.getDistInterval()));
		while((!foundMon)&&(counter<10)){
			Prob = Math.random()*pTable[binChoice][2];
			int i = 0;
			while(Prob>pTable[binChoice][i+3]){
				i++;
			}
			double flip=Math.random();
			if(flip>0.5){
				monCandidate = xlpossibles[fil].myFilament.getMonomerAt(arcLengths[fil]+i*Actin.monLength);
			} else {
				monCandidate = xlpossibles[fil].myFilament.getMonomerAt(arcLengths[fil]-i*Actin.monLength);
			}
			//  NOTE: Must test monCandidate == null because we don't avoid checking for
			//  monomers beyond the end of a filament.   Should do this more cleanly...
			if (monCandidate != null && monCandidate.isFreeLinker()){
				foundMon = true;
				if(linker.site1Bound()) {
					linker.setBound2(monCandidate);
				}  else {
					linker.setBound1(monCandidate);
				}
				return;
			}
			counter++;
		}
	}

	
	/** Check whether the crosslinker has the right angle with respect to the actin to allow binding. */
	private static boolean ckActinVsXLAngle(Actin actin, Crosslinker linker,Actin otherActin) {
		double cosAngle = Point2D.dot(actin.uVect,otherActin.uVect);
		if(cosAngle > linker.getMinBindingAngle() || -cosAngle > linker.getMinBindingAngle()) {
			return true;
		}
		else return false;
	}

}

