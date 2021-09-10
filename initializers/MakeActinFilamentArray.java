package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.Random;

import io.*;
import main.*;

public class MakeActinFilamentArray  extends Initializer
{
	int numFilaments;
	
	double initXPosition;
	
	double initYPosition;
	
	double initAngle;
	
	double xSpacing;
	
	double ySpacing;
	
	double angleIncrement;
	
	FilamentAttachment minusAttachmentTemplate = null;
	
	FilamentAttachment plusAttachmentTemplate = null;
	
	/** Random number generator to generate initial filament lengths. */
	Random lengthGen = new Random((long)Math.random());
	
	public void init() throws Exception  {
		double initX = Sim2D.xDimension/2+initXPosition;
		double initY = Sim2D.yDimension/2+initYPosition;
		double fixedAng = initAngle+Math.PI/2;
		double xuvect = Math.cos(fixedAng);
		double yuvect = Math.sin(fixedAng);
		
		FilamentAttachment mAtt, pAtt;
		
		int monomers;
		for(int i = 0; i < Actin.numInitialActins; i++) {
			if(Actin.useUniformInitialLengthDistribution) {
				monomers = (int)(Actin.avgActinLength/Actin.monLength);
			}
			else {
				monomers = (int)((Actin.avgActinLength*(1 + Actin.stdDevActinLength*Math.abs(lengthGen.nextGaussian())))/Actin.monLength);
			}
			
			
			if(minusAttachmentTemplate != null) {
				mAtt = minusAttachmentTemplate.makeNewInstance();
			}
			else mAtt = null;
			
			if(plusAttachmentTemplate != null) {
				pAtt = plusAttachmentTemplate.makeNewInstance();
			}
			else pAtt = null;
		
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			new Actin(initX, initY,fixedAng, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
			initX += xSpacing;
			initY += ySpacing;
			fixedAng += angleIncrement;
			xuvect = Math.cos(fixedAng);
			yuvect = Math.sin(fixedAng);
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("NumFilaments"))  {
			Actin.numInitialActins = in.nextInt();
		}
		else if(tag.equals("InitXPosition"))  {
			initXPosition = in.nextDouble();
		}
		else if(tag.equals("InitYPosition"))  {
			initYPosition = in.nextDouble();
		}
		else if(tag.equals("InitAngle"))  {
			initAngle = in.nextDouble();
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
		else if(tag.equals("avgActinLength"))  {
			Actin.avgActinLength = in.nextDouble();
		}
		else if(tag.equals("stdDevActinLength"))  {
			Actin.stdDevActinLength = in.nextDouble();
		}
		else if(tag.equals("useUniformInitialLengthDistribution"))  {
			Actin.useUniformInitialLengthDistribution = in.nextBoolean();
		}
		else if(tag.equals("Static"))  {
			if(in.nextBoolean()) {
				Actin.setActinDynamicsMode(Actin.STATIC_FILAMENTS);
			}
		}
		else if(tag.equals("MinusEndAttachment"))  {
			String attach_name = "no_name";
			try {
				attach_name = in.nextString();
				Class c = Class.forName("main." + attach_name);
				minusAttachmentTemplate = (FilamentAttachment)c.newInstance();
				minusAttachmentTemplate.load(in);
			} catch(Exception e) {
				System.out.println("MakeActinFilamentArray.loadParameter: Problem loading MinusAttachment " + attach_name + ":" + e.toString());
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
			System.out.println("MakeActinFilamentArray.loadParameter: Problem loading PlusAttachment " + attach_name + ":" + e.toString());
				throw(e);
			}
		}
		else throw new Exception("MakeActinFilamentArray.loadParameter(): got bad tag = " + tag);
	}
}


