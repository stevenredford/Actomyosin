package collision;


/**
 * MyoVCollisionDetector.java
 *
 * @author Created by Omnicore CodeGuide
 */

import main.*;
import util.*;

public class MyoVCollisionDetector extends CollisionDetector
{


	/** Diagnostic: Number of times a myosinV checked an actin filament for collision */
	public int myoVCks = 0;

	/** Diagnostic: Number of times a myosinV got close enough to possibly bind an actin filament */
	public int myoVHits = 0;
	
	/** Diagnostic: Number of times a myosinV bound an actin filament */
	public int myoVBounds = 0;

	/** Used in binder/actin collision checking to ensure each filament checked only once.*/
	boolean[] actinChecked = new boolean[5000];
	
	/** Keeps track of which actins were checked in a certain bin so we can reset actinChecked efficiently. */
	int[] checkedActins = new int[500];
	
	int checkedActinCt;
	
	/** Number of possible binding partners for the current thing. */
	private int numPossibles;
	
	private BindingInfo [] myoVPossibles = new BindingInfo[1000];    // for crosslinkers

	public void checkCollisions(double simTime) {
		checkMyosinVMesh();
	}

	/** run through myosin mesh, checking each grid element for registered myosins.  For each, consider possible binding partners. */
	public void checkMyosinVMesh(){
		MyosinV myosinV;
		
		for(int i = 0; i < MyosinV.myosinVCt; i++) {
			numPossibles = 0;
			myosinV = MyosinV.theMyosinVs[i];
			if(myosinV.canBind()){
				sampleActinMeshTargetsForMyosinV(myosinV);
			}
			checkFilamentCandidates(myosinV);
		}
	}

	private void sampleActinMeshTargetsForMyosinV(MyosinV myoV) {
	
		Actin actin;
		checkedActinCt = 0;
		int x, y;
		Point2D intPt = new Point2D();
		BindingInfo info;

		Point2D center = myoV.mfAttachmentPt;
		
		int startBinX = Mesh.getBinX(center.x - MyosinV.minBindingDistance);
		int stopBinX = Mesh.getBinX(center.x + MyosinV.minBindingDistance);
		int startBinY = Mesh.getBinY(center.y - MyosinV.minBindingDistance);
		int stopBinY = Mesh.getBinY(center.y + MyosinV.minBindingDistance);
		
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
							myoVCks++;
							info = new BindingInfo();
							
							if (actinCollision(actin,intPt,myoV.mfAttachmentPt,MyosinV.minBindingDistance, info,0)) {
								if(ckActinVsMyoVAngle(actin,myoV)) {
									getMyoVBindingSite(actin,info);
									myoVHits++;
									if (!info.monomer.isFreeMyosinV()){
										myoVPossibles[numPossibles++] = info;
									}
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
	
	
		
	private void checkFilamentCandidates(MyosinV myoV) {
		Actin actin;
		if(numPossibles>0) {
			int choice = (int)(Math.random()*numPossibles);
			if(choice == numPossibles) choice--;
			actin = myoVPossibles[choice].monomer.myFilament;
			if (Math.random() < (myoV.getMyosinVBindingProb(actin) * Sim2D.getDeltaT())) {
				myoV.bindMon(myoVPossibles[choice].monomer);
				myoVBounds++;
			}
		}
	}
	
	/** Given the position of closest Point2D on actin filament in arc length and the distance to that closest Point2D
	 * from this myosinV head (stored in info), find the preferred binding position alo */
	private void getMyoVBindingSite(Actin actin, BindingInfo info) {
		info.arcL += Math.sqrt(MyosinV.minBindingDistance*MyosinV.minBindingDistance - info.minD*info.minD);
		if(info.arcL > actin.length) {
			info.arcL = actin.length;
		}
		info.monomer = actin.getMonomerAt(info.arcL);
	}
	
	
	/** Check whether the myosin V has the right angle with respect to the actin to allow binding. */
	public boolean ckActinVsMyoVAngle(Actin actin, MyosinV m) {
		double cosAngle = Point2D.dot(m.getAim(),actin.uVect);
		if(cosAngle > MyosinV.myoVAngleBias) {
			return true;
		}
		else return false;
	}

}


