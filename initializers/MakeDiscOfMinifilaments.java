package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import io.*;
import main.*;


public class MakeDiscOfMinifilaments  extends Initializer
{
	double discRadius = 1.0;
		
	public void init() {
		double initX;
		double initY;
		double randomAng;
		
		for(int i = 0; i < MyosinMiniFilament.numInitialMinifilaments; i++) {
			initX = (2*Math.random()-1)*discRadius;
			initY = (2*Math.random()-1)*discRadius;
			while (Math.sqrt(initX*initX + initY*initY) > discRadius) {
				initY = (2*Math.random()-1)*discRadius;
			}
			initX += Sim2D.xDimension/2;
			initY += Sim2D.yDimension/2;
			randomAng = 2*Math.PI*Math.random();
			MyosinMiniFilament.makeMiniFilament(initX, initY, Math.cos(randomAng), Math.sin(randomAng));
			
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("NumFilaments"))  {
			MyosinMiniFilament.numInitialMinifilaments = in.nextInt();
		} else if (tag.equals("discRadius")) {
				discRadius = in.nextDouble();
				if (discRadius > Sim2D.xDimension) { throw new Exception("Unworkable discRadius in MakeDiscOfMinifilaments!!"); }
		}
		else throw new Exception("MakeFixedFilaments.loadParameter(): got bad tag = " + tag);
	}
}

