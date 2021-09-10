/**
 * Link.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;

import util.Point2D;
import analysis.*;
import io.*;

abstract public class Link
{

	ForceSubstrate  substrate1 = null;
	
	ForceSubstrate  substrate2 = null;
	
	double localPt1;
	
	double localPt2;
	
	/** point holds the current force */
	Point2D F = new Point2D();
	
	public Link(ForceSubstrate s1, ForceSubstrate s2, double c1, double c2) {
		substrate1 = s1;
		substrate2 = s2;
		localPt1 = c1;
		localPt2 = c2;
	}
	
	public void breakLink() {
		if(substrate1 != null) substrate1.linkDetached(this);
		if(substrate2 != null) substrate2.linkDetached(this);
		substrate1 = substrate2 = null;
		localPt1 = localPt2 = -1;
	}

	public Point2D getForce() {
		return F;
	}
	
	abstract public void doForces();
}

