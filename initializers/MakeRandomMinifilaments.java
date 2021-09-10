package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import io.*;
import main.*;


public class MakeRandomMinifilaments  extends Initializer
{
		
	public void init() {
		double initX;
		double initY;
		double randomAng;
		
		for(int i = 0; i < MyosinMiniFilament.numInitialMinifilaments; i++) {
			initX = Math.random()*Sim2D.xDimension;
			initY = Math.random()*Sim2D.yDimension;
			randomAng = 2*Math.PI*Math.random();
			MyosinMiniFilament.makeMiniFilament(initX, initY, Math.cos(randomAng), Math.sin(randomAng));
			
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("NumMiniFilaments"))  {
			MyosinMiniFilament.numInitialMinifilaments = in.nextInt();
		}
		else throw new Exception("MakeFixedFilaments.loadParameter(): got bad tag = " + tag);
	}
}

