package collision;

/**
 * MyosinVTailCollisionDetector.java
 *
 * @author Created by Omnicore CodeGuide
 */

import main.*;
import util.*;

public class MyosinVTailCollisionDetector extends CollisionDetector
{
	
	/** KLUDGE */
	public static double bigNumber = 1e20;

	/** Diagnostic: Number of times a myosinVtail checked a barrier for collision */
	public static int mVTailBarrierCks = 0;

	/** Diagnostic: Number of times a myosinVtail got close enough to possibly reflect off a barrier */
	public static int mVTailBarrierHits = 0;
	
	/** Diagnostic: Number of times a myosinVtail reflected off a barrier */
	public static int myoVReflects = 0;
	
	/** Used in myosinVtail/myosinVtail collision checking to ensure each myosinVtail checked only once.*/
	boolean[] MVTailChecked = new boolean[5000];
	
	/** Keeps track of which myosinVtails were checked for a certain myosinVtails
	 * so we can reset MVTailChecked efficiently. */
	int[] checkedMVTails = new int[500];
	
	int checkedMVTailCt;

	public void checkCollisions(double simTime) {
		checkMVTailCollisions();
	}
	
	/** run through MVTail mesh, checking each grid element for registered MVTails.
	 * For each, consider possible collision partners. */
	public void checkMVTailCollisions(){
		for(int i = 0; i < MyosinVTail.myosinVTailCt; i++) {
			checkMVTailMesh(MyosinVTail.theMyosinVTails[i]);
		}
	}


	/** run through myosin V tail mesh, checking each grid element for registered myosinsV.  For each, consider possible collision partners. */
	public void checkMVTailMesh(MyosinVTail tail){
	
		int x, y;
		int thatTailID;
		checkedMVTailCt = 0;
		
		Point2D center = tail.cm;
		
		int startBinX = Mesh.getBinX(center.x - MyosinV.myoVCollisionRadius);
		int stopBinX = Mesh.getBinX(center.x + MyosinV.myoVCollisionRadius);
		int startBinY = Mesh.getBinY(center.y - MyosinV.myoVCollisionRadius);
		int stopBinY = Mesh.getBinY(center.y + MyosinV.myoVCollisionRadius);
		
		// Sample all intersected bins in MVTAIL_MESH for candidate tails,
		// making sure every tail is checked only once...
		for(int xx=startBinX;xx<=stopBinX;xx++){
			for(int yy=startBinY;yy<=stopBinY;yy++){
				
				x = Mesh.wrapX(xx); y = Mesh.wrapY(yy);

				if(Mesh.MVTAIL_MESH.timeStamps[x][y]==Sim2D.counter) {
					
					for(int i=0;i<Mesh.MVTAIL_MESH.activeCts[x][y];i++){
						thatTailID = (int)Mesh.MVTAIL_MESH.meshpoints[x][y][i];
						if(tail.myMyosinVTailNumber > thatTailID) {
							if (!MVTailChecked[thatTailID]) {
								MVTailChecked[thatTailID]=true;
								checkedMVTails[checkedMVTailCt]=thatTailID;
								checkedMVTailCt++;
								myosinVTailCollision(tail,MyosinVTail.theMyosinVTails[thatTailID]);
							}
						}
					}
				}
			}
		}
		// reset the actins that were checked for this MyosinVtail
		for(int j=0;j<checkedMVTailCt;j++) {
			MVTailChecked[checkedMVTails[j]]=false;
		}
	}

	
	// Check MVTail mesh against Barrier mesh
	public void checkBarrierCrossings(){
//		System.out.println("Collision.checkBarrierMesh: method started...");
		Point2D intPt = new Point2D();
		MyosinVTail mvtail;
		for(int x=0;x<Mesh.nXBins;x++){
			for(int y=0;y<Mesh.nYBins;y++) {
				if(Mesh.MVTAIL_MESH.timeStamps[x][y]-1==Sim2D.counter){
					for(int i=0;i<Mesh.MVTAIL_MESH.activeCts[x][y];i++){
						int mvtailID=(int)Mesh.MVTAIL_MESH.meshpoints[x][y][i];
						mvtail=MyosinVTail.theMyosinVTails[mvtailID];
//						System.out.println("Collision.checkForBarrierCrossings: checking mvtail "+mvtailID+" at mesh pt (x,y) = ("+x+","+y+") at time = "+Sim2D.simulationTime+" s");
						checkBarrierMesh(x,y,mvtail);
					}
				}
			}
		}
//		System.out.println("Collision.checkBarrierMesh: method complete...");
	}
	
	public void checkBarrierMesh(int x,int y,MyosinVTail mvtail){
		BarrierElement barrier;
		if(((Mesh.BARRIER_MESH.timeStamps[x][y])-1)!=Sim2D.counter){ return; }
		
		for(int i=0;i<Mesh.BARRIER_MESH.activeCts[x][y];i++){

			int barrierId = (int)Mesh.BARRIER_MESH.meshpoints[x][y][i];
			barrier=BarrierElement.theBarrierElements[barrierId];
			mVTailBarrierCks++;
								
			checkThisMVTailVsThisBarrier(x, y, barrier, mvtail);
		}
	}

	public void checkThisMVTailVsThisBarrier(int x, int y, BarrierElement barrier, MyosinVTail mvtail) {
		ReflectInfo info = new ReflectInfo();
		Point2D intersectPt = new Point2D();
		getMVTailPositions(mvtail,info);

		info.barrierX1 = barrier.aEnd.x;
		info.barrierY1 = barrier.aEnd.y;
		info.barrierX2 = barrier.bEnd.x;
		info.barrierY2 = barrier.bEnd.y;
		if (barrierCrossing(barrier,intersectPt,info)) {

			mVTailBarrierHits++;
			findReflectedPosition(barrier,info);
			correctMVTailPosition(x,y,mvtail,info);

		}
	}
	
	/** Check whether a given step to pt crosses a barrier **/
	public static boolean barrierCrossing(BarrierElement barrier, Point2D ptOfIntersect, ReflectInfo info) {
		Point2D lastPt = new Point2D(info.lastX,info.lastY);
		Point2D nextPt = new Point2D(info.nextX,info.nextY);
		Point2D barrierPt1 = new Point2D(info.barrierX1,info.barrierY1);
		Point2D barrierPt2 = new Point2D(info.barrierX2,info.barrierY2);
		
		checkStepCrossLine(barrierPt1, barrierPt2, lastPt, nextPt, ptOfIntersect);
		
		if (ptOfIntersect.x != bigNumber && ptOfIntersect.x != -bigNumber) {
			return true;
		} else {
		return false;
		}
	}
		
	/** check if a step from stepPt1 to stepPt2, crosses a line defined by linePt1 and linePt2 */
	public static double checkStepCrossLine(Point2D p1, Point2D p2, Point2D p3, Point2D p4, Point2D intPt) {

		  double xD1,yD1,xD2,yD2,xD3,yD3;
		  double dot,deg,len1,len2;
		  double segmentLen1,segmentLen2;
		  double ua,ub,div;

		  // calculate differences
		  xD1=p2.x-p1.x;
		  xD2=p4.x-p3.x;
		  yD1=p2.y-p1.y;
		  yD2=p4.y-p3.y;
		  xD3=p1.x-p3.x;
		  yD3=p1.y-p3.y;

		  // calculate the lengths of the two lines
		  len1=Math.sqrt(xD1*xD1+yD1*yD1);
		  len2=Math.sqrt(xD2*xD2+yD2*yD2);

		  // calculate angle between the two lines.
		  dot=(xD1*xD2+yD1*yD2); // dot product
		  deg=dot/(len1*len2);

		  // if abs(angle)==1 then the lines are parallel,
		  // so no intersection is possible
		  if(Math.abs(deg)==1) {
			  intPt.set(-bigNumber,-bigNumber);
			  return 0;
		  }

		  // find intersection Pt between two lines
		  div=yD2*xD1-xD2*yD1;
		  ua=(xD2*yD3-yD2*xD3)/div;
		  ub=(xD1*yD3-yD1*xD3)/div;
		  intPt.x=p1.x+ua*xD1;
		  intPt.y=p1.y+ua*yD1;

		  // calculate the combined length of the two segments
		  // between Pt-p1 and Pt-p2
		  xD1=intPt.x-p1.x;
		  xD2=intPt.x-p2.x;
		  yD1=intPt.y-p1.y;
		  yD2=intPt.y-p2.y;
		  segmentLen1=Math.sqrt(xD1*xD1+yD1*yD1)+Math.sqrt(xD2*xD2+yD2*yD2);

		  // calculate the combined length of the two segments
		  // between Pt-p3 and Pt-p4
		  xD1=intPt.x-p3.x;
		  xD2=intPt.x-p4.x;
		  yD1=intPt.y-p3.y;
		  yD2=intPt.y-p4.y;
		  segmentLen2=Math.sqrt(xD1*xD1+yD1*yD1)+Math.sqrt(xD2*xD2+yD2*yD2);

		  // if the lengths of both sets of segments are the same as
		  // the lengths of the two lines the point is actually
		  // on the line segment.

		  // if the point isn’t on the line, return null
		  if (Math.abs(len1-segmentLen1)>0.01 || Math.abs(len2-segmentLen2)>0.01) {
		    intPt.set(-bigNumber,-bigNumber);
		  	return 0;
		  }
		  return ua;
	}
	
	/** Check whether the mvtail lies close enough to the barrier segment delimited by segStart and segStop to possibly reflect.
	 * If so, returns loads info packet with distance to intersection point and arclength of that point
	 * along the filament segment and returns true.  Otherwise returns false.
	 **/
	public static boolean barrierSegCross(Point2D segStart, Point2D segStop, ReflectInfo info){
		Point2D lastPt = new Point2D(info.lastX,info.lastY);
		Point2D nextPt = new Point2D(info.nextX,info.nextY);
		Point2D ptOfIntersect = new Point2D();
		System.out.println("Collision.barrierSegCross: method running...");
		checkStepCrossLine(segStop, segStop, lastPt, nextPt, ptOfIntersect);

		if (ptOfIntersect.x == 99999 || ptOfIntersect.x == -99999) {
//			System.out.println("Collision.barrierCrossing: barrier crossed!");
			return false;
		}
//		System.out.println("Collision.barrierCrossing: barrier NOT crossed!");
		return true;
	}

	public static void getMVTailPositions(MyosinVTail mvtail,ReflectInfo info) {
		info.lastX = mvtail.lastXpos;
		info.lastY = mvtail.lastYpos;
		info.nextX = mvtail.nextXpos;
		info.nextY = mvtail.nextYpos;
//		System.out.println("Collision.getLastMVTailPosition(): mvtail "+mvtail.myMyosinVTailNumber+" got positions at time "+Sim2D.simulationTime+" s");
	}
	
	public static void findReflectedPosition(BarrierElement barrier, ReflectInfo info) {
		double theta;
		double tanratio;
		if ((info.barrierX2 - info.barrierX1) == 0){
			theta = Math.PI/2;
		} else {
			tanratio = (info.barrierY2 - info.barrierY1)/(info.barrierX2 - info.barrierX1);
			theta = Math.atan(tanratio);
		}
		double xoffset = info.barrierX1;
		double yoffset = info.barrierY1;
		double xin = info.nextX - xoffset;
		double yin = info.nextY - yoffset;
		double xout = xin*Math.cos(2*theta) + yin*Math.sin(2*theta);
		double yout = xin*Math.sin(2*theta) - yin*Math.cos(2*theta);
		info.reflectedX = xout + xoffset;
		info.reflectedY = yout + yoffset;
//		System.out.println("Collision.findReflectedPosition: theta = "+(2*3.14*theta)+", xoffset = "+xoffset+", xin = "+xin+", xout = "+xout);
	}
	
	public static void correctMVTailPosition(int x,int y,MyosinVTail mvtail, ReflectInfo info) {
		checkMVTailCorrection(x,y,mvtail,info);
		info.xCorrection = info.reflectedX - info.nextX;
		info.yCorrection = info.reflectedY - info.nextY;
		mvtail.reflectPositionOffBarrier(info.xCorrection,info.yCorrection);
//		System.out.println("Collision.correctMVTailPosition: xCorr = "+info.xCorrection+", yCorr = "+info.yCorrection);
	}
	
	
	public static void checkMVTailCorrection(int x, int y, MyosinVTail mvtail, ReflectInfo info) {
		BarrierElement barrier;
		if(((Mesh.BARRIER_MESH.timeStamps[x][y])-1)!=Sim2D.counter){ return; }
		
		for(int i=0;i<Mesh.BARRIER_MESH.activeCts[x][y];i++){

			int barrierId = (int)Mesh.BARRIER_MESH.meshpoints[x][y][i];
			barrier=BarrierElement.theBarrierElements[barrierId];
			mVTailBarrierCks++;
								
			Point2D intersectPt = new Point2D();
			
			info.barrierX1 = barrier.aEnd.x;
			info.barrierY1 = barrier.aEnd.y;
			info.barrierX2 = barrier.bEnd.x;
			info.barrierY2 = barrier.bEnd.y;
				
			Point2D lastPt = new Point2D(info.lastX,info.lastY);
			Point2D reflectPt = new Point2D(info.reflectedX,info.reflectedY);
			Point2D barrierPt1 = new Point2D(info.barrierX1,info.barrierY1);
			Point2D barrierPt2 = new Point2D(info.barrierX2,info.barrierY2);
			
			checkStepCrossLine(barrierPt1, barrierPt2, lastPt, reflectPt, intersectPt);
			
			if (intersectPt.x != bigNumber && intersectPt.x != -bigNumber) {
				return;
			} else {
				info.reflectedX = info.lastX;
				info.reflectedY = info.lastY;
				return;
			}
		}
		
	}

	public static void myosinVTailCollision(MyosinVTail m1,
		MyosinVTail m2) {
		if (m1.offCortex() | m2.offCortex()) {
			return;
		}
		if (Point2D.getDistanceSqrd(m1.cm, m2.cm) > MyosinVTail.collisionRadSqrd) {
			return;
		}
		double impD = MyosinV.myoVCollisionRadius
				- Point2D.getDistance(m1.cm, m2.cm);
		double mag = (0.01 * impD / Sim2D.getDeltaT())
				/ (1 / m1.myosinVtailTransGamma + 1 / m2.myosinVtailTransGamma);
		m1.tempP.getVector(m2.cm, m1.cm);
		m1.tempP.scale(mag);
		m1.addForce(m1.tempP);
		m1.tempP.scale(-1);
		m2.addForce(m1.tempP);
	}

}

class ReflectInfo {
	double lastX = -1;
	double lastY = -1;
	double nextX = -1;
	double nextY = -1;
//	double crossX = -1;
//	double crossY = -1;
	double barrierX1 = -1;
	double barrierY1 = -1;
	double barrierX2 = -1;
	double barrierY2 = -1;
	double reflectedX = -1;
	double reflectedY = -1;
	double xCorrection = -1;
	double yCorrection = -1;
}

