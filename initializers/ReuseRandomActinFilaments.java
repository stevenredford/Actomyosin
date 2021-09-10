package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */
import java.util.Random;

import io.*;
import main.*;

public class ReuseRandomActinFilaments  extends MakeRandomActinFilaments
{
		
	static double [][] cList;
	
	public void init() throws Exception {
		makeFilaments();
	}
	
	public void makeFilaments() throws Exception {
		FilamentAttachment mAtt, pAtt;
		int monomers;
		for(int i = 0; i < Actin.numInitialActins; i++) {
			if(Actin.useUniformInitialLengthDistribution) {
				monomers = (int)(Actin.avgActinLength/Actin.monLength);
			}
			else {
				monomers = (int)((Actin.avgActinLength + Actin.stdDevActinLength*lengthGen.nextGaussian())/Actin.monLength);
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

			new Actin (cList[i][0], cList[i][1],cList[i][2], monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
		}
	}

	static private void makeFilamentCoordinates() {
		cList = new double[Actin.numInitialActins][3];
		for(int i = 0; i < Actin.numInitialActins; i++) {
			cList[i][0] = Math.random()*Sim2D.xDimension;
			cList[i][1] = Math.random()*Sim2D.yDimension;
			cList[i][2] = 2*Math.PI*Math.random();
		}
	}
			
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("TargetDensity"))  {
			throw new Exception("ReuseRandomActinFilaments only makes a fixed number of filaments");
		}
		else if(tag.equals("TargetNumber"))  {
			Actin.numInitialActins = in.nextInt();
			makeFilamentCoordinates();
		}
		else super.loadParameter(tag,in);
	}
}

