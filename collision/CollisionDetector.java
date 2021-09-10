package collision;


import util.*;
import main.*;

/**
 * CollisionDetector.java
 *
 * @author Created by Omnicore CodeGuide
 */


abstract public class CollisionDetector
{

	abstract public void checkCollisions(double simTime);
	
	public void checkBarrierCrossings() {}
	
	public void init() {}

	/**
	 * Method lineAndPoint
	 *
	 * Computes point on line segment (p1->p2) closest to point pP.  Assigns point coordinates to intPt and
	 * returns the position of intPt along (p1->p2), scaled so that points along (p1->p2) lie between 0 and 1..
	 *
	 * @param    p1                  a  Point
	 * @param    p2                  a  Point
	 * @param    pP                  a  Point
	 * @param    intPt               a  Point
	 *
	 * @return   a double
	 *
	 */
	public static double lineAndPoint (Point2D p1, Point2D p2, Point2D pP, Point2D intPt) {
		intPt.set(-1,-1);
		
		double xDiff = Point2D.getDiffX(p1,p2);
		double yDiff = Point2D.getDiffY(p1,p2);
		double magSqrd = Math.pow(xDiff,2) + Math.pow(yDiff,2);
		
		double u = ((pP.x-p1.x)*xDiff + (pP.y-p1.y)*yDiff)/magSqrd;
		
		intPt.set(p1.x + u*xDiff, p1.y + u*yDiff);
		intPt.wrapPoint();
		return u;
	}
	

	
	/** Determine whether an intersection exists between two line segments defined by (p11, p12) and (p21, p22)
	 * if so, store the result in intPt and returns the position of intPt along targetsegment (t1,t2),
	 * scaled so that positions on the segment lie in [0,1].  If not, return -1.
	 */
	public static double lineSegmentIntersectionTest(Point2D p11, Point2D p12,Point2D p21, Point2D p22, Point2D intPt) {
		Point2D v1 = Point2D.diff(p12,p11);
		Point2D v2 = Point2D.diff(p22,p21);
		double vCross = Point2D.Cross(v2,v1);
		if(vCross == 0) return -1;
		Point2D pDiff = Point2D.diff(p21,p11);
		double lam1 = Point2D.Cross(v2,pDiff)/vCross;
		double lam2 = Point2D.Cross(v1,pDiff)/vCross;
		if(lam1>=0 && lam1<=1 && lam2>= 0 && lam2<=1)  { // intersection occurred
			intPt.set(p11);
			intPt.inc(lam1,v1);
			return lam1;
		}
		return -1;
	}
		
	/** Check whether the given Point2D pt lies within distance a minD of the filament by checking each segment of the
	 * fliament in turn. If so, set the coordinates of intPt, store distance to intersection point and its arc length
	 * along the filament in info and return true. If not, returns false.
	 **/
	public static boolean actinCollision(Actin actin,Point2D intPt, Point2D pt, double minD, BindingInfo info, double overHang) {
		switch (actin.segCt) {
		case 1:
			if(actinSegCollision(actin.segs[0][0],actin.segs[0][1],intPt,pt, minD,info)) {
				info.arcL *= actin.segLength[0];
				return true;
			}
			return false;
		case 2:
			if(actinSegCollision(actin.segs[0][0],actin.segs[0][1],intPt, pt, minD,info)) {
				if(info.arcL <=1) {  // if > 1, its on the next segment
					info.arcL *= actin.segLength[0];
					if(info.arcL > -overHang) {
						return true;
					}
				}
			}
			if(actinSegCollision(actin.segs[1][0],actin.segs[1][1],intPt, pt, minD,info)) {
				info.arcL = actin.segLength[0] + info.arcL*actin.segLength[1];
				return true;
			}
			return false;
		case 3:
			if(actinSegCollision(actin.segs[0][0],actin.segs[0][1],intPt, pt, minD,info)) {
				if(info.arcL <=1) { // if > 1, its on the next segment
					info.arcL *= actin.segLength[0];
					if(info.arcL > -overHang) {
						return true;
					}
					return true;
				}
			}
			if(actinSegCollision(actin.segs[1][0],actin.segs[1][1],intPt, pt, minD,info)) {
				if(info.arcL <=1) { // if > 1, its on the next segment
					info.arcL = actin.segLength[0] + info.arcL*actin.segLength[1];
					return true;
				}
			}
			if(actinSegCollision(actin.segs[2][0],actin.segs[2][1],intPt, pt, minD,info)) {
				info.arcL = actin.segLength[0] + actin.segLength[1] + info.arcL*actin.segLength[2];
				return true;
			}
			return false;
		}
		return false;
	}
	

	/** Check whether the line segment specified by (pt1,pt2) intersects the (possibly multi-segmented actin filament
	 * The line segment must not cross the coundaries of the periodic domain. If an intersection is found,
	 * store the coordiantes of the intersestion point in intPt and return the distance along the filament
	 * from its pointed end to the intersection point.If not, return -1.
	 **/
	public static double actinSegCollision(Actin actin, Point2D pt1, Point2D pt2, Point2D intPt) {
		double arcLength = -1;
		double arcLNormed = 0;
		switch (actin.segCt) {
		case 1:
			arcLNormed = lineSegmentIntersectionTest(actin.segs[0][0],actin.segs[0][1],pt1,pt2,intPt);
			if (arcLNormed >=0) {
				arcLength = arcLNormed*actin.segLength[0];
			}
			return arcLength;
		case 2:
			arcLNormed = lineSegmentIntersectionTest(actin.segs[0][0],actin.segs[0][1],pt1,pt2,intPt);
			if (arcLNormed >=0) {
				arcLength = arcLNormed*actin.segLength[0];
				return arcLength;
			}
			arcLNormed = lineSegmentIntersectionTest(actin.segs[1][0],actin.segs[1][1],pt1,pt2,intPt);
			if (arcLNormed >=0) {
				arcLength = actin.segLength[0] + arcLNormed*actin.segLength[1];
			}
			return arcLength;
		case 3:
			arcLNormed = lineSegmentIntersectionTest(actin.segs[0][0],actin.segs[0][1],pt1,pt2,intPt);
			if (arcLNormed >=0) {
				arcLength = arcLNormed*actin.segLength[0];
				return arcLength;
			}
			arcLNormed = lineSegmentIntersectionTest(actin.segs[1][0],actin.segs[1][1],pt1,pt2,intPt);
			if (arcLNormed >=0) {
				arcLength = actin.segLength[0] + arcLNormed*actin.segLength[1];
				return arcLength;
			}
			arcLNormed = lineSegmentIntersectionTest(actin.segs[2][0],actin.segs[2][1],pt1,pt2,intPt);
			if (arcLNormed >=0) {
				arcLength = actin.segLength[0] + actin.segLength[1] + arcLNormed*actin.segLength[2];
			}
			return arcLength;
		}
		return arcLength;
	}
	
	/** Check whether the myosin lies close enough to the filament segment delimited by segStart and segStop to possibly bind.
	 * If so, returns loads info packet with distance to intersection point and arclength of that point
	 * along the filament segment and returns true.  Otherwise returns false.
	 **/
	public static boolean actinSegCollision(Point2D segStart, Point2D segStop,Point2D intPt, Point2D pt, double minD, BindingInfo info){
		info.arcL = lineAndPoint(segStart,segStop,pt,intPt);
		if ((info.minD = Point2D.getDistance(intPt, pt)) < minD) return true;
		return false;
	}
	

}

class BindingInfo {
	
	double arcL = -1;
	double minD = -1;
	Monomer monomer = null;
	Actin actin = null;
}


