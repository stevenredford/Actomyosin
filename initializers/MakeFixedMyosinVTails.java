package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import io.*;
import main.*;


public class MakeFixedMyosinVTails  extends Initializer
{
	int numMyosinVTails;
	
	double initXPosition;
	
	double initYPosition;
	
	double xSpacing;
	
	double ySpacing;
	
	double angleIncrement;
	
	boolean pinned = false;
	
	public void init() {
		double initX = Sim2D.xDimension/2+initXPosition;
		double initY = Sim2D.yDimension/2+initYPosition;
		double fixedAng = Math.PI/2;
		double xuvect = Math.cos(fixedAng);
		double yuvect = Math.sin(fixedAng);
		
		for(int i = 0; i < numMyosinVTails; i++) {
			MyosinVTail.makeMyosinVTail(initX, initY, xuvect, yuvect);
			initX += xSpacing;
			initY += ySpacing;
			fixedAng += angleIncrement;
			xuvect = Math.cos(fixedAng);
			yuvect = Math.sin(fixedAng);
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("NumMyosinVTails"))  {
			numMyosinVTails = in.nextInt();
		}
		else if(tag.equals("InitXPosition"))  {
			initXPosition = in.nextDouble();
		}
		else if(tag.equals("InitYPosition"))  {
			initYPosition = in.nextDouble();
		}
		else if(tag.equals("XSpacing"))  {
			xSpacing = in.nextDouble();
		}
		else if(tag.equals("YSpacing"))  {
			ySpacing = in.nextDouble();
		}
		else if(tag.equals("AngleIncrement"))  {
			angleIncrement = in.nextDouble();
		}
		else if(tag.equals("Pinned"))  {
			pinned = in.nextBoolean();
		}
		else throw new Exception("MakeMyosinVTails.loadParameter(): got bad tag = " + tag);
	}
}

