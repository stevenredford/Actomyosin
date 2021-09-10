package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;
import util.*;

public class MakeParallelFilamentArray  extends Initializer
{
	
	static public Point2D leftEndpoint = new Point2D();
	
	static public Point2D rightEndpoint = new Point2D();
	
	int numFilaments;
	
	int numAnchorsLeft;
	
	int numAnchorsRight;
	
	int counterLeft=1;
	
	int counterRight=1;
	
	double anchorSpacingLeft;
	
	double anchorSpacingRight;
	
	double averageLength;
	
	double stdvLength;
	
	FilamentAttachment minusAttachmentTemplate = null;
	
	FilamentAttachment plusAttachmentTemplate = null;
		
	double initY;

	
	Random lengthGen = new Random((long)Math.random());
	
	public void init() throws Exception {
		
			for(int i = 0; i < numFilaments; i++) {
			
				double randomAng = Math.random() < 0.5 ? 0:Math.PI;
				int monomers = 0;
					while (monomers <= 0) {
						double length=averageLength + stdvLength*lengthGen.nextGaussian();
						monomers = (int)(length/Actin.monLength);
					}
					double initX = Math.random()*Sim2D.xDimension;
				if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
				new Actin (initX, initY,randomAng, monomers,monomers,true,null,null);
			
		}
		
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("leftXPosition"))  {
			leftEndpoint.x = in.nextDouble();
		}
		else if(tag.equals("leftYPosition"))  {
			leftEndpoint.y = in.nextDouble();
		}
		else if(tag.equals("rightXPosition"))  {
			rightEndpoint.x = in.nextDouble();
		}
		else if(tag.equals("rightYPosition"))  {
			rightEndpoint.y = in.nextDouble();
		}
		else if(tag.equals("averageLength"))  {
			averageLength = in.nextDouble();
		}
		else if(tag.equals("stdvLength"))  {
			stdvLength = in.nextDouble();
		}
		else if(tag.equals("numFilaments"))  {
			numFilaments = in.nextInt();
		}
		else if(tag.equals("MinusEndAttachment"))  {
			String attach_name = "no_name";
			try {
				attach_name = in.nextString();
				Class c = Class.forName("main." + attach_name);
				minusAttachmentTemplate = (FilamentAttachment)c.newInstance();
				minusAttachmentTemplate.load(in);
			} catch(Exception e) {
				System.out.println("MakeParallelFilamentArray.loadParameter: Problem loading MinusAttachment " + attach_name + ":" + e.toString());
				throw(e);
			}
		}
		else if(tag.equals("PlusEndAttachment"))  {
			String attach_name = "no_name";
			try {
				attach_name = in.nextString();
				Class c = Class.forName("main." + attach_name);
				plusAttachmentTemplate = (FilamentAttachment)c.newInstance();
				plusAttachmentTemplate.load(in);
			} catch(Exception e) {
			System.out.println("MakeParallelFilamentArray.loadParameter: Problem loading PlusAttachment " + attach_name + ":" + e.toString());
				throw(e);
			}
		}
		else if(tag.equals("initY"))  {
			initY = in.nextDouble();
		}
		else if(tag.equals("numAnchorsLeft"))  {
			numAnchorsLeft = in.nextInt();
		}
		else if(tag.equals("numAnchorsRight"))  {
			numAnchorsRight = in.nextInt();
		}
		else if(tag.equals("anchorSpacingLeft"))  {
			anchorSpacingLeft = in.nextDouble();
		}
		else if(tag.equals("anchorSpacingRight"))  {
			anchorSpacingRight = in.nextDouble();
		}
		
	
		else throw new Exception("MakeParallelFilamentArray.loadParameter(): got bad tag = " + tag);
	}
	

}

