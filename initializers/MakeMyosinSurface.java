package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.Random;

import io.*;
import main.*;
import util.*;
import parameters.*;

public class MakeMyosinSurface extends Initializer
{
	double x1;
	double y1;
	double x2;
	double y2;
	Point2D startPt,stopPt;
	double myoDensity = 0;
	double myoLineDensity = 0;
	boolean uniform = true;
	String myoDensityMode = null;
	
	
	public void init() {
		
		if(myoDensityMode != null) {
			if(myoDensityMode.equalsIgnoreCase(MyosinSurface.SURFACE_DENSITY)) {
				if (myoDensity != 0) {
					Parameters.setParameter("myoSurfaceDensity",myoDensity);
				}
				new MyosinSurface(x1+Sim2D.xDimension/2,y1+Sim2D.yDimension/2,x2,y2);
			}
			else { // LINE_DENSITY
				if (myoLineDensity != 0) {
					Parameters.setParameter("myoLineDensity",myoLineDensity);
				}
				startPt = new Point2D(x1+Sim2D.xDimension/2,y1+Sim2D.yDimension/2);
				stopPt = new Point2D(x2+Sim2D.xDimension/2,y2+Sim2D.yDimension/2);
				new MyosinSurface(startPt,stopPt,uniform);
			}
		}

		else if (myoDensity != 0) {
			Parameters.setParameter("myoSurfaceDensity",myoDensity);
			new MyosinSurface(x1+Sim2D.xDimension/2,y1+Sim2D.yDimension/2,x2,y2);
		}
		
		else if (myoLineDensity != 0) {
			Parameters.setParameter("myoLineDensity",myoLineDensity);
			startPt = new Point2D(x1+Sim2D.xDimension/2,y1+Sim2D.yDimension/2);
			stopPt = new Point2D(x2+Sim2D.xDimension/2,y2+Sim2D.yDimension/2);
			new MyosinSurface(startPt,stopPt,uniform);
		}

	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {

		if(tag.equals("CenterXPosition"))  {
			x1 = in.nextDouble();
		}
		else if(tag.equals("CenterYPosition"))  {
			y1 = in.nextDouble();
		}
		else if(tag.equals("XDimension"))  {
			x2 = in.nextDouble();
		}
		else if(tag.equals("YDimension"))  {
			y2 = in.nextDouble();
		}
		else if(tag.equals("StartX"))  {
			x1 = in.nextDouble();
		}
		else if(tag.equals("StartY"))  {
			y1 = in.nextDouble();
		}
		else if(tag.equals("StopX"))  {
			x2 = in.nextDouble();
		}
		else if(tag.equals("StopY"))  {
			y2 = in.nextDouble();
		}
		else if(tag.equals("MyosinDensity"))  {
			myoDensity = in.nextDouble();
		}
		else if(tag.equals("MyoLineDensity"))  {
			myoLineDensity = in.nextDouble();
		}
		else if(tag.equals("myoDensityMode")) {
			String s = in.nextString();
			if(s.equalsIgnoreCase(MyosinSurface.SURFACE_DENSITY) || s.equalsIgnoreCase(MyosinSurface.LINE_DENSITY)) {
				myoDensityMode = s;
			}
			else throw new Exception("MakeMyosinSurface.loadParameter:  got bad input fro myoDensityMode = " + s);
		}
		
		else throw new Exception("MakeMyosinHolder.loadParameter(): got bad tag = " + tag);
		
	}
}


