/**
 * FixedForceSubstrate.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;

import util.Point2D;


/** A fixed attachment point that never moves and ignores forces..*/
public class FixedForceSubstrate implements ForceSubstrate {
	
	public FixedForceSubstrate(Point2D pt) {
		attachPt = new Point2D(pt);
	}
		
	public Point2D attachPt;
	
	public Point2D getGlobalCoords(double localCoords) {
		return attachPt;
	}
	
	public void addForce(Point2D x, Point2D F) {}

	public void linkDetached(Link l) {}

}

