/**
 * FilamentAttachment.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;


import util.*;
import analysis.*;
import io.*;

import java.text.DecimalFormat;


abstract public class FilamentAttachment  {

		
	/** The monomer this is attached to. */
	Monomer theMonomer;
	
	/** Holds the initial attachment site */
	public Point2D initialAttachmentPoint;
	
	/** Holds the current attachment site */
	public Point2D attachmentPoint;
	
	/** point holds the current pin force */
	Point2D pinForce = new Point2D();

	ValueTracker forceTrack;

	public void getElasticForce() {}
	
	public double getViscousResistance() {
		return 0;
	}

	public double getElasticResistance() {
		return 0;
	}
	
	public void relax() {}

	
	public void setPosition(Point2D p) {
		attachmentPoint = Point2D.Copy(p);
		//attachmentPoint.wrapPoint();
	}
		
	public void setInitialPosition(Point2D p) {
		initialAttachmentPoint = Point2D.Copy(p);
		//initialAttachmentPoint.wrapPoint();
		setPosition(p);
	}
	
	public void restoreInitialPosition() {
		setPosition(initialAttachmentPoint);
		
	}
			
	public void setMonomer(Monomer m) {
		theMonomer = m;
	}
	
	public void detach() {
		forceTrack = null;
		theMonomer = null;
	}
	
	public void initForceTracking() {
		forceTrack = new ValueTracker (100,ValueTracker.POINT_TYPE);
	}
	
	public Point2D getPinForce() {
		return pinForce;
	}
	
	public Point2D getTimeAveragedForce() {
		return forceTrack == null ? new Point2D(0,0) : forceTrack.averagePtVal();
	}

	public void load(AMInputStream in)  throws Exception {
		String tag = in.nextTag();
		while(!tag.equals("endAttachment"))  {
			loadParameter(tag, in);
			tag = in.nextTag();
		}
	}
	
	abstract public void loadParameter(String tag, AMInputStream in)  throws Exception;
	
	
	abstract public FilamentAttachment makeNewInstance()throws Exception;

}

