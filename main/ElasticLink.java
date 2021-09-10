/**
 * ElasticLink.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;

import parameters.*;
import util.Point2D;
import io.*;
import analysis.*;

public class ElasticLink extends Link
{

	static String className = new String("main.ElasticLink");
	
	double springConstant=0;
	
	public ElasticLink(ForceSubstrate s1, ForceSubstrate s2, double c1, double c2) {
		super(s1,s2,c1,c2);
	}
			
	public void doForces() {
		Point2D attachPt1 = substrate1.getGlobalCoords(localPt1);
		Point2D attachPt2 = substrate2.getGlobalCoords(localPt2);
		F.getVector(attachPt1, attachPt2);
		F.scale(springConstant); //now force vector
		substrate1.addForce(attachPt1,F);
		F.scale(-1);
		substrate2.addForce(attachPt2,F);
	}
}

