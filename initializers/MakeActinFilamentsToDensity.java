package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */
import java.util.Random;

import io.*;
import main.*;

public class MakeActinFilamentsToDensity  extends Initializer
{
	
	/** Random number generator to generate initial filament lengths. */
	Random lengthGen = new Random((long)Math.random());
	
	public void init() {
		double curActinDensity = 0;
		double totalMons = 0;
		while (curActinDensity < Actin.targetActinDensity) {
			double initX = Math.random()*Sim2D.xDimension;
			double initY = Math.random()*Sim2D.yDimension;
			double randomAng = 2*Math.PI*Math.random();
			int finalMons = (int)((Actin.avgActinLength*(1 + Actin.stdDevActinLength*lengthGen.nextGaussian()))/Actin.monLength);
			if (finalMons > Actin.maxMon) {
				finalMons = Actin.maxMon;
				FileOps.reportln ("makeActinToADensity.init(): requested too long a filament, setting to max length for this arena of " + finalMons + " monomers");
			}
			if (finalMons < Actin.minMon) { finalMons = Actin.minMon; }
			new Actin (initX, initY,randomAng, finalMons, finalMons,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),null,null);
			totalMons += finalMons;
			curActinDensity = totalMons/(Sim2D.xDimension*Sim2D.yDimension*Sim2D.zDimension);
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("avgActinLength"))  {
			Actin.avgActinLength = in.nextDouble();
		}
		else if(tag.equals("stdDevActinLength"))  {
			Actin.stdDevActinLength = in.nextDouble();
		}
		else if(tag.equals("TargetDensity"))  {
			Actin.targetActinDensity = in.nextDouble();
		}
		else throw new Exception("MakeActinFilamentsToDensity.loadParameter(): got bad tag = " + tag);
	}
}

