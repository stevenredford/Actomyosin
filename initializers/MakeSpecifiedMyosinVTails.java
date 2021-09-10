package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class MakeSpecifiedMyosinVTails  extends Initializer
{
	Vector templates = new Vector();
	
	
	public void init() {
		for(Enumeration e = templates.elements(); e.hasMoreElements(); ) {
			((Template)e.nextElement()).init();
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("MyosinVTail"))  {
			Template tt = new Template();
			tt.loadTemplate(in);
			templates.add(tt);
		}
		else throw new Exception("MakeSpecifiedMyosinVTails.loadParameter(): got bad tag = " + tag);
	}
	
	class Template {
		
		double xPosition;
		double yPosition;
		double orientation;
		boolean pinned;
		
		public void loadTemplate(AMInputStream in)  throws Exception {
			
			String tag = in.nextTag();
			
			while(!tag.equals("endMyosinVTail")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}
				else if(tag.equals("InitYPosition"))  {
					yPosition = in.nextDouble();
				}
				else if(tag.equals("Orientation"))  {
					orientation = Math.PI*in.nextDouble()/360.0;
				}
				else if(tag.equals("Pinned"))  {
					pinned = in.nextBoolean();
				}
				else throw new Exception("MakeSpecifiedMyosinVTails.Template.loadParameter(): got bad tag = " + tag);
				tag = in.nextTag();
			}
		}
		
		public void init() {
			double initX = Sim2D.xDimension/2+xPosition;
			double initY = Sim2D.yDimension/2+yPosition;
			double xuvect = Math.cos(orientation);
			double yuvect = Math.sin(orientation);
			MyosinVTail.makeMyosinVTail(initX, initY, xuvect, yuvect);
		}

	}
}

