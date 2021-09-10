package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import io.*;
import parameters.*;

abstract public class Initializer
{
	
	abstract public void init() throws Exception ;
	
	public void load(AMInputStream in)  throws Exception {
		String tag = in.nextTag();
		while(!tag.equals("endInitializer"))  {
			loadParameter(tag, in);
			tag = in.nextTag();
		}
	}
	
	abstract public void loadParameter(String tag, AMInputStream in)  throws Exception;
}

