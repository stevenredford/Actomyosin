package collision;

/**
 * MiniFilamentCollisionDetector.java
 *
 * @author Created by Omnicore CodeGuide
 */

import main.*;
import util.*;


public class MiniFilamentCollisionDetector extends CollisionDetector
{
	
	
	/** Used in MINIFILAMENT/MINIFILAMENT collision checking to ensure each minifilament checked only once.*/
	boolean[] MFChecked = new boolean[5000];
	
	/** Keeps track of which Minifilaments were checked for a certain Minifilament
	 * so we can reset MFChecked efficiently. */
	int[] checkedMFs = new int[5000];
	
	int checkedMFCt;
	

	public void checkCollisions(double simTime) {
		checkMiniFilamentCollisions();
	}

	/** run through myosin minifilament mesh, checking each grid element for registered minifilaments.  For each, consider possible collision partners. */
	public void checkMiniFilamentCollisions(){
		for(int i = 0; i < MyosinMiniFilament.miniFilamentCt; i++) {
			checkedMFCt = 0;
			checkMFMesh(MyosinMiniFilament.theMiniFilaments[i]);
		}
	}
		
		
	private void checkMFMesh(MyosinMiniFilament mFil) {
	
		int x, y;
		int thatMFilID;
		Point2D center = mFil.cm;
		
		int startBinX = Mesh.getBinX(center.x - Myosin.myoCollisionRadius);
		int stopBinX = Mesh.getBinX(center.x + Myosin.myoCollisionRadius);
		int startBinY = Mesh.getBinY(center.y - Myosin.myoCollisionRadius);
		int stopBinY = Mesh.getBinY(center.y + Myosin.myoCollisionRadius);
		
		// Sample all intersected bins in MFIL_MESH for candidate minifilaments,
		// making sure every minifilament is checked only once...
		for(int xx=startBinX;xx<=stopBinX;xx++){
			for(int yy=startBinY;yy<=stopBinY;yy++){
				
				x = Mesh.wrapX(xx); y = Mesh.wrapY(yy);
				
				if(Mesh.MFIL_MESH.timeStamps[x][y]==Sim2D.counter) {
					
					for(int i=0;i<Mesh.MFIL_MESH.activeCts[x][y];i++){
						thatMFilID = (int)Mesh.MFIL_MESH.meshpoints[x][y][i];
						if(mFil.myMiniFilamentNumber > thatMFilID) {
							if (!MFChecked[thatMFilID]) {
								MFChecked[thatMFilID]=true;
								checkedMFs[checkedMFCt]=thatMFilID;
								checkedMFCt++;
								miniFilamentCollision(mFil,MyosinMiniFilament.theMiniFilaments[thatMFilID]);
							}
						}
					}
				}
			}
		}
		// reset the minifilaments that were checked for this minifilament
		for(int j=0;j<checkedMFCt;j++) {
			MFChecked[checkedMFs[j]]=false;
		}
	}
	
	public void miniFilamentCollision(MyosinMiniFilament m1, MyosinMiniFilament m2) {
		if (m1.offCortex() | m2.offCortex()) { return; }
		if (Point2D.getDistanceSqrd(m1.cm, m2.cm) > MyosinMiniFilament.collisionRadSqrd) { return; }
		double impD = Myosin.myoCollisionRadius - Point2D.getDistance(m1.cm, m2.cm);
		double mag = (0.01*impD/Sim2D.getDeltaT())/(1/m1.myosinfilTransGamma+1/m2.myosinfilTransGamma);
		m1.tempP.getVector(m2.cm,m1.cm);
		m1.tempP.scale(mag);
		m1.addForce(m1.tempP);
		m1.tempP.scale(-1);
		m2.addForce(m1.tempP);
	}
	
}

