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

public class MakeDiscOfActinFilaments  extends Initializer
{
	static String
		NUMBER = "NUMBER",
		DENSITY = "DENSITY";
	
	String targetMode;
	double discRadius = 1.0;
	boolean outward = false;
	boolean inward = false;
		
	/** Random number generator to generate initial filament lengths. */
	Random lengthGen = new Random((long)Math.random());
	
	FilamentAttachment minusAttachmentTemplate = null;
	
	FilamentAttachment plusAttachmentTemplate = null;

	public void init() throws Exception {
		
		if(targetMode.equals(NUMBER)) {
			makeFilamentsToNumber();
		}
		else if(targetMode.equals(DENSITY)) {
			makeFilamentsToDensity();
		}
	}
	
	public void makeFilamentsToDensity() throws Exception {
		double curActinDensity = 0;
		double totalMons = 0;
		double angle;
		
		FilamentAttachment mAtt=null, pAtt=null;
		
		while (curActinDensity < Actin.targetActinDensity) {
			double initX = (2*Math.random()-1)*discRadius;
			double initY = (2*Math.random()-1)*discRadius;
			while (Math.sqrt(initX*initX + initY*initY) > discRadius) {
				initY = (2*Math.random()-1)*discRadius;
			}
			initX += Sim2D.xDimension/2;
			initY += Sim2D.yDimension/2;
			if (outward) {
				angle = getAnglePointingOutward(new Point2D(initX,initY));
			} else if (inward) {
				angle = getAnglePointingInward(new Point2D(initX,initY));
			} else {
				angle = 2*Math.PI*Math.random();
			}
			int finalMons = (int)((Actin.avgActinLength*(1 + Actin.stdDevActinLength*lengthGen.nextGaussian()))/Actin.monLength);
			if (finalMons > Actin.maxMon) {
				finalMons = Actin.maxMon;
				//FileOps.reportln ("makeActinToADensity.init(): requested too long a filament, setting to max length for this arena of " + finalMons + " monomers");
			}
			if (finalMons < Actin.minMon) { finalMons = Actin.minMon; }
			
			if(minusAttachmentTemplate != null) {
				mAtt = minusAttachmentTemplate.makeNewInstance();
			} else mAtt = null;
			if(plusAttachmentTemplate != null) {
				pAtt = plusAttachmentTemplate.makeNewInstance();
			} else pAtt = null;
		
			new Actin (initX, initY,angle, finalMons, finalMons,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
			totalMons += finalMons;
			curActinDensity = totalMons/(discRadius*discRadius*Math.PI*Sim2D.zDimension);
			
		}
	}

	private void makeFilamentsToNumber() throws Exception {
		double angle;
		FilamentAttachment mAtt=null, pAtt=null;
		
		for(int i = 0; i < Actin.numInitialActins; i++) {
			double initX = (2*Math.random()-1)*discRadius;
			double initY = (2*Math.random()-1)*discRadius;
			while (Math.sqrt(initX*initX + initY*initY) > discRadius) {
				initY = (2*Math.random()-1)*discRadius;
			}
			initX += Sim2D.xDimension/2;
			initY += Sim2D.yDimension/2;
			if (outward) {
				angle = getAnglePointingOutward(new Point2D(initX,initY));
			} else if (inward) {
				angle = getAnglePointingInward(new Point2D(initX,initY));
			} else {
				angle = 2*Math.PI*Math.random();
			}
			int monomers = 0;
			
			if(Actin.useUniformInitialLengthDistribution) {
				monomers = (int)(Actin.avgActinLength/Actin.monLength);
			}
			else {
				while (monomers <= 0) {
					monomers = (int)((Actin.avgActinLength*(1 + Actin.stdDevActinLength*Math.abs(lengthGen.nextGaussian())))/Actin.monLength);
				}
			}
			
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			
			if(minusAttachmentTemplate != null) {
				mAtt = minusAttachmentTemplate.makeNewInstance();
			}
			else mAtt = null;
			
			if(plusAttachmentTemplate != null) {
				pAtt = plusAttachmentTemplate.makeNewInstance();
			}
			else pAtt = null;
			
			new Actin (initX, initY,angle, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
		}
	}
	
	private double getAnglePointingOutward (Point2D aPt) {
		Point2D cPt = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
		//double ptDist = Point2D.getDistance(cPt, aPt);
		double xDiffNorm = Point2D.getDiffX(cPt, aPt);
		double yDiffNorm = Point2D.getDiffY(cPt, aPt);
		double angOut = Math.atan2(yDiffNorm,xDiffNorm);
		return angOut;
	}
	
	private double getAnglePointingInward (Point2D aPt) {
		Point2D cPt = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
		//double ptDist = Point2D.getDistance(aPt, cPt);
		double xDiffNorm = Point2D.getDiffX(aPt, cPt);
		double yDiffNorm = Point2D.getDiffY(aPt, cPt);
		double angOut = Math.atan2(yDiffNorm,xDiffNorm);
		return angOut;
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("TargetDensity"))  {
			Actin.targetActinDensity = in.nextDouble();
		}
		else if(tag.equals("TargetNumber"))  {
			Actin.numInitialActins = in.nextInt();
		}
		else if(tag.equals("TargetMode"))  {
			targetMode = (in.nextString()).toUpperCase();
			if(!(targetMode.equals("NUMBER") || targetMode.equals("DENSITY")))
				throw new Exception("MakeRandomFilaments.loadParameter(): got bad targetMode = " + targetMode);
		}
		else if(tag.equals("discRadius"))  {
			discRadius = in.nextDouble();
			if (discRadius > Sim2D.xDimension) { System.out.println("Unworkable discRadius in MakeDiscOfActinFilaments!!"); System.exit(0); }
		}
		else if (tag.equals("outward")) {
			if(in.nextBoolean()) {
				outward = true;
			}
		}
		else if (tag.equals("inward")) {
			if(in.nextBoolean()) {
				inward = true;
			}
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
				System.out.println("MakeDiscOfActinFilaments.loadParameter: Problem loading MinusAttachment " + attach_name + ":" + e.toString());
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
			System.out.println("MakeDiscOfActinFilaments.loadParameter: Problem loading PlusAttachment " + attach_name + ":" + e.toString());
				throw(e);
			}
		}
		else throw new Exception("MakeDiscOfActinFilaments.loadParameter(): got bad tag = " + tag);
	}
}

