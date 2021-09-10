package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class MakeSpecifiedBarrierElements  extends Initializer
{
	
	Vector templatevecs = new Vector();
	Vector templatecoords = new Vector();
	
	
	public void init() {
		for(Enumeration e = templatevecs.elements(); e.hasMoreElements(); ) {
			((TemplateVec)e.nextElement()).initVec();
		}
		for(Enumeration f = templatecoords.elements(); f.hasMoreElements(); ) {
			((TemplateCoords)f.nextElement()).initCoords();
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {

		if(tag.equals("BarrierVector") || tag.equals("BarrierCoords")) {
			if(tag.equals("BarrierVector"))  {
				TemplateVec zz = new TemplateVec();
				zz.loadTemplateVec(in);
				templatevecs.add(zz);
			}
			if(tag.equals("BarrierCoords"))  {
				TemplateCoords z = new TemplateCoords();
				z.loadTemplateCoords(in);
				templatecoords.add(z);
			}
		}
		else throw new Exception("MakeSpecifiedBarrierElements.loadParameter(): got bad tag = " + tag);
	}
	
	class TemplateVec {
		
		double xPosition;
		double yPosition;
		double orientation;
		double length;
		boolean pinned;
		
		public void loadTemplateVec(AMInputStream in)  throws Exception {
			
			String tag = in.nextTag();
			
			while(!tag.equals("endBarrierVector")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}
				else if(tag.equals("InitYPosition"))  {
					yPosition = in.nextDouble();
				}
				else if(tag.equals("Length"))  {
					length = in.nextDouble();
				}
				else if(tag.equals("Orientation"))  {
					orientation = 2.0*Math.PI*in.nextDouble()/360.0;
				}
				else if(tag.equals("Pinned"))  {
					pinned = in.nextBoolean();
				}
				else throw new Exception("MakeSpecifiedBarrierElements.Template.loadParameter(): got bad tag (Vector) = " + tag);
				tag = in.nextTag();
			}
		}
		
		public void initVec() {
			double initX = Sim2D.xDimension/2+xPosition;
			double initY = Sim2D.yDimension/2+yPosition;
			new BarrierElement(initX, initY,orientation,length,pinned);
		}
	

	}

	class TemplateCoords {
		
		double xA;
		double yA;
		double xB;
		double yB;
		
		public void loadTemplateCoords(AMInputStream in)  throws Exception {
			
			String tag = in.nextTag();
			
			while(!tag.equals("endBarrierCoords")) {
				if(tag.equals("XA"))  {
					xA = in.nextDouble();
				}
				else if(tag.equals("YA"))  {
					yA = in.nextDouble();
				}
				else if(tag.equals("XB"))  {
					xB = in.nextDouble();
				}
				else if(tag.equals("YB"))  {
					yB = in.nextDouble();
				}
				else throw new Exception("MakeSpecifiedBarrierElements.Template.loadParameter(): got bad tag (Coords) = " + tag);
				tag = in.nextTag();
			}
		}
	
		public void initCoords() {
			new BarrierElement(xA,yA,xB,yB);
		}
	}
}
	

