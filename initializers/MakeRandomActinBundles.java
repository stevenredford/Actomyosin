package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */
import java.util.Random;

import io.*;
import main.*;

public class MakeRandomActinBundles  extends Initializer
{
	static String
		NUMBER = "NUMBER",
		DENSITY = "DENSITY";
	
	String targetMode;
	
	FilamentAttachment attachmentTemplate = null;
		
	/** Random number generator to generate initial filament lengths. */
	Random lengthGen = new Random((long)Math.random());
	
	public void init()  throws Exception {
		
		if(targetMode.equals(NUMBER)) {
			makeBundlesToNumber();
		}
		else if(targetMode.equals(DENSITY)) {
			makeBundlesToDensity();
		}
	}
	
	public void makeBundlesToDensity() throws Exception {
		double curActinDensity = 0;
		double totalMons = 0;
		FilamentAttachment att;
		
		while (curActinDensity < Actin.targetActinDensity) {
			double initX = Math.random()*Sim2D.xDimension;
			double initY = Math.random()*Sim2D.yDimension;
			double randomAng = 2*Math.PI*Math.random();
			int finalMons = (int)((Actin.avgActinLength*(1 + Actin.stdDevActinLength*Math.abs(lengthGen.nextGaussian())))/Actin.monLength);
			if (finalMons > Actin.maxMon) {
				finalMons = Actin.maxMon;
				System.out.println ("makeActinToADensity.init(): requested too long a filament, setting to max length for this arena of " + finalMons + " monomers");
			}
			if (finalMons < Actin.minMon) { finalMons = Actin.minMon; }
			int initMons = (int)(Math.random()*(finalMons-5) + 5);
			
			if(attachmentTemplate != null) {
				att = attachmentTemplate.makeNewInstance();
			} else att = null;
			new Bundle (initX, initY,randomAng, initMons, finalMons,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),att,att);
			totalMons += finalMons;
			curActinDensity = totalMons/(Sim2D.xDimension*Sim2D.yDimension*Sim2D.zDimension);
		}
	}

	private void makeBundlesToNumber() throws Exception {
		FilamentAttachment att;
		for(int i = 0; i < Actin.numInitialActins; i++) {
			double initX = Math.random()*Sim2D.xDimension;
			double initY = Math.random()*Sim2D.yDimension;
			double randomAng = 2*Math.PI*Math.random();
			int monomers;
			
			if(Actin.useUniformInitialLengthDistribution) {
				monomers = (int)(Actin.avgActinLength/Actin.monLength);
			}
			else {
				monomers = (int)((Actin.avgActinLength*(1 + Actin.stdDevActinLength*Math.abs(lengthGen.nextGaussian())))/Actin.monLength);
			}
			
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			
			if(attachmentTemplate != null) {
				att = attachmentTemplate.makeNewInstance();
			}
			else att = null;
	
			new Bundle (initX, initY,randomAng, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),att,att);
		}
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
		else if(tag.equals("EndAttachment"))  {
			String attach_name = "no_name";
			try {
				attach_name = in.nextString();
				Class c = Class.forName("main." + attach_name);
				attachmentTemplate = (FilamentAttachment)c.newInstance();
				attachmentTemplate.load(in);
			} catch(Exception e) {
			System.out.println("MakeRandomActinBundles.loadParameter: Problem loading Attachment " + attach_name + ":" + e.toString());
				throw(e);
			}
		}
		else throw new Exception("MakeRandomActinBundles.loadParameter(): got bad tag = " + tag);
	}
}

