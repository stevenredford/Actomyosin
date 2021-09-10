package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class MakeSpecifiedMinifilaments  extends Initializer
{
	Vector templates = new Vector();
	
	
	public void init() {
		for(Enumeration e = templates.elements(); e.hasMoreElements(); ) {
			((Template)e.nextElement()).init();
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("MiniFilament"))  {
			Template t = new Template();
			t.loadTemplate(in);
			templates.add(t);
		}
		else throw new Exception("MakeSpecifiedMinifilaments.loadParameter(): got bad tag = " + tag);
	}
	
	class Template {
		
		double xPosition;
		double yPosition;
		double orientation;
		boolean pinned = false;
		double viscosityTweak=1;
		double myoRateMult=1.0;
		
		public void loadTemplate(AMInputStream in)  throws Exception {
			
			String tag = in.nextTag();
			
			while(!tag.equals("endMiniFilament")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}
				else if(tag.equals("InitYPosition"))  {
					yPosition = in.nextDouble();
				}
				else if(tag.equals("Orientation"))  {
					orientation = Math.PI*in.nextDouble()/180.0;
				}
				else if(tag.equals("Pinned"))  {
					pinned = in.nextBoolean();
				}
				else if(tag.equalsIgnoreCase("viscosityTweak"))  {
					viscosityTweak = in.nextDouble();
				}
				else throw new Exception("MakeSpecifiedMinifilaments.Template.loadParameter(): got bad tag = " + tag);
				tag = in.nextTag();
			}
		}
		
		public void init() {
			double initX = Sim2D.xDimension/2+xPosition;
			double initY = Sim2D.yDimension/2+yPosition;
			double xuvect = Math.cos(orientation);
			double yuvect = Math.sin(orientation);
			MyosinMiniFilament newbie = new MyosinMiniFilament(initX, initY, xuvect, yuvect,pinned);
			newbie.setTweakDragFactor(viscosityTweak);
		}

	}
}

