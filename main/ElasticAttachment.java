/**
 * ElasticAttachhment.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;

import parameters.*;
import util.Point2D;
import io.*;

public class ElasticAttachment extends FilamentAttachment
{

	static String className = new String("main.ElasticAttachment");
	static double variableSpringConstant=0;
	double springConstant=0;
	
	public double getElasticResistance() {
		return springConstant;
	}
	
	public void relax() {
		setPosition(theMonomer.getLocation());
	}

	public void getElasticForce() {
//		double pinStretch = Point2D.getDistance(theMonomer.getLocation(), attachmentPoint);
//		pinForce.getVector(theMonomer.getLocation(), attachmentPoint);
		double pinStretch = Point2D.getDistanceIgnorePeriod(theMonomer.getLocation(), attachmentPoint);
		pinForce.getVectorIgnorePeriod(theMonomer.getLocation(), attachmentPoint);
		pinForce.scale(springConstant); //now force vector
		theMonomer.addForce(pinForce);
		forceTrack.registerValue(pinForce);
	}

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equalsIgnoreCase("SpringConstant")) {
			springConstant = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("VariableSpringConstant")) {
			String parName = in.nextString();
			if(!Parameters.isValidParameterName(parName)) {
				throw new Exception("ElasticAttachment.loadParameter(): tried to assign a variable spring constant parameter anme that doesn't exist");
			}
			Parameters.addParameterListener(parName,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						springConstant = p.getDoubleValue();
					}
				}
			);
			springConstant = Parameters.getParameter(parName).getDoubleValue();
		}
		else throw new Exception("ElasticAttachment.loadParameter(): got bad tag = " + tag);
	}

	public FilamentAttachment makeNewInstance()throws Exception {
		ElasticAttachment ea = this.getClass().newInstance();
		ea.springConstant = this.springConstant;
		return ea;
	}
	static void staticInit() {

		Parameters.addParameter(className,"variableSpringConstant",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				variableSpringConstant=p.getDoubleValue();
			}
		}
		);

	}

}

