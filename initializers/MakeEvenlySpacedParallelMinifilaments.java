package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import io.*;
import main.*;


public class MakeEvenlySpacedParallelMinifilaments  extends Initializer
{
	int numFilaments;
	double initY;
	double leftEndpoint;
	double rightEndpoint;
	public void init() {
		
		
	
		
//		for(int i = 0; i < numFilaments; i++) {
		if(MyosinMiniFilament.numInitialMinifilaments == 1) {
			double initX = 0.5*(leftEndpoint+rightEndpoint);
			MyosinMiniFilament.makeMiniFilament(initX, initY, 0, 0);
		}
		else {
			for(int i = 0; i < MyosinMiniFilament.numInitialMinifilaments; i++) {
				double initX = leftEndpoint+i*(rightEndpoint-leftEndpoint)/(MyosinMiniFilament.numInitialMinifilaments-1);
				MyosinMiniFilament.makeMiniFilament(initX, initY, 0, 0);
			}
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("numFilaments"))  {
			numFilaments = in.nextInt();
		}
		else if (tag.equals("initY"))  {
			initY = in.nextDouble();
		}
		else if (tag.equals("leftXPosition"))  {
			leftEndpoint = in.nextDouble();
		}
		else if (tag.equals("rightXPosition"))  {
			rightEndpoint = in.nextDouble();
		}
		else throw new Exception("MakeFixedFilaments.loadParameter(): got bad tag = " + tag);
	}
}

