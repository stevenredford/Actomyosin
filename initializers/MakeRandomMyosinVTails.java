package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import io.*;
import main.*;

public class MakeRandomMyosinVTails  extends Initializer
{
	int numMyosinVTails;							// boundary defaults:
	boolean useBoundaries = false;					// 1) don't use boundaries.
	double xLowerBound = -Sim2D.xDimension/2;	// 2) if using boundaries, each
	double xUpperBound = Sim2D.xDimension/2;	//    boundary must be specified
	double yLowerBound = -Sim2D.yDimension/2;	//    or else just us the edge
	double yUpperBound = Sim2D.yDimension/2;	//    the world.
		
	public void init() {
		double initX;
		double initY;
		double randomAng;
		double xRange = xUpperBound - xLowerBound;
		double yRange = yUpperBound - yLowerBound;
		
		if (useBoundaries == true) {
			for(int i = 0; i < numMyosinVTails; i++) {
				initX = Math.random()*xRange + xLowerBound - Sim2D.xDimension/2;
				initY = Math.random()*yRange + yLowerBound - Sim2D.yDimension/2;
				randomAng = 2*Math.PI*Math.random();
				MyosinVTail.makeMyosinVTail(initX, initY, Math.cos(randomAng), Math.sin(randomAng));
			}
			
		} if (useBoundaries == false) {
			for(int i = 0; i < numMyosinVTails; i++) {
				initX = Math.random()*Sim2D.xDimension;
				initY = Math.random()*Sim2D.yDimension;
				randomAng = 2*Math.PI*Math.random();
				MyosinVTail.makeMyosinVTail(initX, initY, Math.cos(randomAng), Math.sin(randomAng));
			}
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("NumMyosinVTails"))  {
			numMyosinVTails = in.nextInt();
		} else if(tag.equals("XLowerBound"))  {
			xLowerBound = in.nextDouble();
		} else if(tag.equals("XUpperBound"))  {
			xUpperBound = in.nextDouble();
		} else if(tag.equals("YLowerBound"))  {
			yLowerBound = in.nextDouble();
		} else if(tag.equals("YUpperBound"))  {
			yUpperBound = in.nextDouble();
		} else if(tag.equals("UseBoundaries"))  {
			useBoundaries = in.nextBoolean();
		} else throw new Exception("MakeRandomMyosinVTails.loadParameter(): got bad tag = " + tag);
	}
	
}

