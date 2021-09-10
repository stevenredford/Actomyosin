/**
 * ForceSubstrate.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;

import util.Point2D;

public interface ForceSubstrate
{
	
	public Point2D getGlobalCoords(double localCoords);
	
	public void addForce(Point2D x, Point2D F);
	
	public void linkDetached(Link l);
	
}

