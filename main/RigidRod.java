package main;


import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.text.DecimalFormat;
import java.util.*;
import java.awt.geom.AffineTransform;
import java.awt.Graphics2D;
import sun.tools.tree.ThisExpression;

import analysis.*;
import util.*;
import collision.*;
import io.*;
import iterators.GlidingAssayEvaluator;
import gui.*;
import parameters.*;

/** A Rigid Rod that can segment itself with respect to periodic boundary conditions, tally up forces
 * and compute its own motions.  Extended by Actin.and BarrierElement.  */
public class RigidRod extends Thing {
	
	
	public int forceTextLocator = 1;  // 1 for above ends, -1 for below
	
	/** radius of the rod in nm. */
	public double radius = 2.0;
	 
	/** Total length of this rod from end to end.Used in calculations of viscosity etc.. */
	public double physicalLength;
		
	Point2D end1;
	
	Point2D end2;
	
	/** The maximum length allowable for a rod to avoid confusion with periodic boundary conditions. */
	static public int maxLength;
	
	/** Unit vector parallel to long axis of filament. */
	public Point2D uVect = new Point2D ();
	
	/** Unit vector normal to long axis of filament */
	public Point2D upVect = new Point2D ();

	/** PERIODIC BOUNDARIES: Number of segments this filament is divided into. */
	public int segCt = 0;
	
	
	/** KEEP TRACK OF DISPLACEMENTS OVER TIME */
	public Point2D lastCM = new Point2D();

	/** Maximum displacement of filament endpoints since last call to getMaxMoveSinceLastCheck() */
	double getMaxMoveSinceLastCheck = 0;
	
	/** Cumulative displacement of center of mass since last call to getMaxMoveSinceLastCheck() */
	Point2D cumDeltaCM = new Point2D();
	
	/** unit vector orthogonal to filament at last call to getMaxMoveSinceLastCheck() */
	Point2D refUPVect = new Point2D();
	
	/** Cumulative roptation of filament since last call to getMaxMoveSinceLastCheck() */
	double cumDeltaAngle = 0;

	/* PERIODIC BOUNDARIES: Endpoints of segments. */
	public Point2D [][] segs;
	
	/** PERIODIC BOUNDARIES: Segment lengths. */
	public double [] segLength = new double[3];
	

	/**********  FORCE CALCULATIONS  ***********/
	
	/* Reused for torque calculations. */
	Point2D torqArm = new Point2D();
	
	/* Reused in translational force calculations. */
	Point2D pTobVec = new Point2D();
	
	/* drag coefficients, needed for calculating diffusion */
	double untweakedActinfilTransGamma = 0;
	double untweakedActinfilPerpGamma = 0;
	double untweakedRotationalGamma = 0;
	
	double tweakDragFactor = 1;
	
	/* drag coefficients, needed for calculating diffusion */
	double actinfilTransGamma = 0;
	double actinfilPerpGamma = 0;
	double rotationalGamma = 0;
	
	/* For calculating diffusion */
	double diffusionParallel = 0;
	double diffusionPerpendicular = 0;
	double diffusionRotational = 0;
	
	double diffusionStepParallel = 0;
	double diffusionStepPerpendicular = 0;
	double diffusionStepRotational = 0;
	
	/* Random number generator to generate diffusive steps. */
	Random generator = new Random(Sim2D.generator.nextLong());
	
	/* For calculating forces due to diffusion */
	Point2D brownianForceParallel=new Point2D();
	Point2D brownianForcePerpendicular=new Point2D();
	Point2D brownianForceRotational=new Point2D();
	
	/* For summing forces */
	Point2D forceSum=new Point2D();
	
	/* Velocity of center of mass. */
	public Point2D velocity=new Point2D();
	
	/** total roational stiffness.  Used each timestep to determine ministep size. */
	double rotK = 0;
	
	/** total translational stiffness parallel to filament. Used each timestep to determine ministep size. */
	double transParK = 0;

	/** total translational stiffness orthogonal to filament. Used each timestep to determine ministep size. */
	double transOrthoK = 0;

	/* Accumulates torque on this rod each timestep. */
	double torque = 0;
	
	/* Accumulates drag on this rod each timestep. */
	double totalDrag = 0;
	
	/* Brownian torque at each timestep. */
	double brownianTorque = 0;
	
	/* change in angle this timestep. */
	double deltaTheta = 0;
	
	/* Current angle. */
	double thetaNOW = 0;
	
	/* Next angle */
	double thetaNEXT = 0;
	
	/* angular velocity this timestep */
	double angularVel = 0;

	Color myColor;
	
	static DecimalFormat lengthFormat = new DecimalFormat ("#000.#; #000.#");	// a clean format nm lengths
	static DecimalFormat forceFormat = new DecimalFormat ("#0.000#; #0.000#");	// a clean format nm lengths



	
	/**
	 * Constructor
	 *
	 * @param    initX               a  double  = initial X position
	 * @param    initY               a  double = initial Y position
	 * @param    initAng             a  double = initial angle
	 *
	 */
	public RigidRod (double initX, double initY, double initAng, double r) {
		this(initX,initY,new Point2D(Math.cos(initAng), Math.sin(initAng)),r);
	}
	
	/**
	 * Constructor
	 *
	 * @param    initX               a  double  = initial X position
	 * @param    initY               a  double = initial Y position
	 * @param    uV             a  double = initial unit Vector
	 *
	 */
	public RigidRod (double initX, double initY, Point2D uV, double r ) {
		super (initX,initY);
		end1 = new Point2D(); end2 = new Point2D();
		radius = r;
		uVect = new Point2D(uV);
		uVect.uVect();
		upVect.orthogonalVector(uVect);
		for (int i=0;i<3;i++) { segLength[i] = 0; }
	}
		
	
	public void setPhysicalLength(double l) {
		physicalLength = l;
	}
	
	/** Determine the positions of pointed and barbed ends from center of mass and orientation (µVect) */
	public void useCenterMassAndUVectToDetermineEnds(){
		end1.add(cm,-physicalLength/2,uVect);
		end2.add(cm,physicalLength/2,uVect);
	}
		
	
	//. FORCES **** FORCES **** FORCES **** FORCES **** FORCES **** FORCES **** FORCES **** FORCES
	
	/** Add the force to the filament's center of mass, thus no torque. */
	public void addForce(Point2D f){
		forceSum.inc(f);
	}
	
	/** Add the force at the specified location, incrementing both the net force on the center of mass, and the net torque.*/
	public void addForce (Point2D f, Monomer m, double orthoK, double parK){
		addForce(f);
		torqArm.getVector(cm,m.cm);
		addTorque(torqArm,f);
		transParK += parK;
		transOrthoK += orthoK;
		rotK += orthoK*f.length();
		
	}

	/** Add the force at the specified location, incrementing both the net force on the center of mass, and the net torque.*/
	public void addForce (Point2D f, Point2D loc){
		addForce(f);
		torqArm.getVector(cm,loc);
		addTorque(torqArm,f);
	}
	
	
	/** Calculate a torque based on force f and torque arm r and add it to net torque for this timestep. */
	public void addTorque(Point2D r, Point2D f) {
		torque += Point2D.Cross(r,f);
	}
	

	
	// INTEGRATION STEP ***** INTEGRATION STEP ***** INTEGRATION STEP ***** INTEGRATION STEP *****
	
	public void step(double dT) {
		forceBasedMotion(dT);
		updateMaxMove(dT);
	}
		
	private double getMaxStep() {
		double maxStep = Math.max(actinfilTransGamma/transParK,actinfilPerpGamma/transOrthoK);
		return 0.5*Math.max(maxStep,rotationalGamma/rotK);
	}
		
	private void updateMaxMove(double dT) {
		cumDeltaCM.inc(dT,velocity);
		cumDeltaAngle += deltaTheta;
		Point2D cumDisp = new Point2D();
		double scale = 0.5*physicalLength*cumDeltaAngle;
		cumDisp.add(cumDeltaCM,scale,refUPVect);
		double moveB = cumDisp.length();
		cumDisp.add(cumDeltaCM,-scale,refUPVect);
		double moveP = cumDisp.length();
		double max = Math.max(moveP,moveB);
		if(max > getMaxMoveSinceLastCheck) getMaxMoveSinceLastCheck = max;
	}


	public double getMaxMoveSinceLastCheck() {
		cumDeltaCM.zero();
		cumDeltaAngle = 0;
		refUPVect.set(upVect);
		double ret = getMaxMoveSinceLastCheck;
		getMaxMoveSinceLastCheck = 0;
		return ret;
	}

	/** all the calculations of movement and rotation due to impinging forces */
	public void forceBasedMotion(double dT) {
		
		//***********************************************************************************************
		// PARALLEL DIFFUSION FORCES - moving in the direction of uVect
		// you are going to diffuse with this step size in the direction of uVect
		diffusionStepParallel = Math.sqrt(2*diffusionParallel*dT)* generator.nextGaussian();
		
		//Convert the step in parallel direction to x,y components of forces in the universal coordinate system
		brownianForceParallel.x = uVect.x * ((actinfilTransGamma / dT) * diffusionStepParallel);
		brownianForceParallel.y = uVect.y * ((actinfilTransGamma / dT) * diffusionStepParallel);
		
		
		
		//***********************************************************************************************
		//PERPENDICULAR DIFFUSION FORCES - moving in the direction of upVect
		diffusionStepPerpendicular = Math.sqrt(2*diffusionPerpendicular*dT)* generator.nextGaussian();
		
		brownianForcePerpendicular.x =upVect.x * ((actinfilPerpGamma / dT) * diffusionStepPerpendicular);
		brownianForcePerpendicular.y = upVect.y * ((actinfilPerpGamma / dT) * diffusionStepPerpendicular);
		
		
		
		//TOTAL TRANSLATIONAL FORCES**********************************************************************
		
		//calculate the total force acting on the center of mass
		//System.out.println("transXforce before Brownian and constant = " + totalTranslationalForce.x);
		double brownianTweak = 1.0;
		if (Sim2D.brownianOn) {
			forceSum.x += brownianTweak*(brownianForceParallel.x + brownianForcePerpendicular.x);
			forceSum.y += brownianTweak*(brownianForceParallel.y + brownianForcePerpendicular.y);
		}
		
		totalDrag += actinfilTransGamma;

		velocity.x = forceSum.x / totalDrag;
		velocity.y = forceSum.y / actinfilTransGamma;
		
		//velocity.x = 0;
		//velocity.y=0; 
		
		// Move the center of mass.
		cm.inc(dT,velocity);
		
		if (!Point2D.pointOK(cm)) {
			FileOps.reportln ("actin cm is NaN: actin = " + this + " v.x = " + velocity.x + " v.y = " + velocity.y);
		}
		
		//ROTATIONAL DIFFUSION FORCES**********************************************************************
		if (Sim2D.brownianOn) {
			diffusionStepRotational = Math.sqrt(2*diffusionRotational*dT)* generator.nextGaussian();
			brownianTorque = rotationalGamma * diffusionStepRotational / dT;
		} else {
			brownianTorque = 0;
		}
		
		torque +=  brownianTweak*brownianTorque;
		
		angularVel = torque / rotationalGamma;
		
		//TOTAL ANGULAR FORCES**************************************************************************

		deltaTheta = angularVel * dT;		//contribution to change in theta from diffusion
					
		thetaNOW = Math.atan2 (uVect.y, uVect.x);
		thetaNEXT = thetaNOW+deltaTheta;

		uVect.x = Math.cos (thetaNEXT);
		uVect.y = Math.sin (thetaNEXT);
		uVect.uVect();
		upVect.orthogonalVector(uVect);
				
		// Finally, calculate the position of the ends
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
	}
	
	public void scaleTweakDragFactor(double s)  {
		tweakDragFactor*=s;
		reTweakViscosities();
	}
	
	public void setTweakDragFactor(double f) {
		tweakDragFactor=f;
		reTweakViscosities();
	}
	
	private void reTweakViscosities()  {
		rotationalGamma = untweakedRotationalGamma*tweakDragFactor;
		actinfilTransGamma = untweakedActinfilTransGamma*tweakDragFactor;
		actinfilPerpGamma = untweakedActinfilPerpGamma*tweakDragFactor;
		
		diffusionParallel = Constants.kT / actinfilTransGamma;
		diffusionPerpendicular = Constants.kT / actinfilPerpGamma;
		diffusionRotational =  Constants.kT / rotationalGamma;
	}
		
	
	/** Compute length dependent drag coefficients. */
	public void evaluateProperties () {
		double naturalLogTerm = (double) Math.log(physicalLength/(2*radius));
		double PiViscosityLength= (double) Math.PI* Constants.viscosity*physicalLength;

		// double pValue = (2*radius)/length;
		// the Sigmas would normally be looked up from a table given the pValue.
		// Jacked from page 107 of Jonathan Howard's Mechanics of Motor Proteins and the Cytoskeleton
		double rotSigma = -0.66;
		double parSigma = -0.20;
		double perpSigma = 0.86;
		
		double pValue = (2*radius)/physicalLength;
		boolean found=false;
		for(int i=0;i<Drag.cylinderDragValues.length&&!found;i++){
			if(pValue>=Drag.cylinderDragValues[i][Drag.pValueIndex]){
				parSigma=Drag.cylinderDragValues[i][Drag.parIndex];
				perpSigma=Drag.cylinderDragValues[i][Drag.perpIndex];
				rotSigma=Drag.cylinderDragValues[i][Drag.rotIndex];
				found=true;
			}
		}
		
		if(!found){
			int defaultIndex=Drag.cylinderDragValues.length-1;
			parSigma=Drag.cylinderDragValues[defaultIndex][Drag.parIndex];
			perpSigma=Drag.cylinderDragValues[defaultIndex][Drag.perpIndex];
			rotSigma=Drag.cylinderDragValues[defaultIndex][Drag.rotIndex];
		}
		untweakedRotationalGamma = (PiViscosityLength*physicalLength*physicalLength)/(3*(naturalLogTerm+rotSigma));
		untweakedActinfilTransGamma = (2*PiViscosityLength)/(naturalLogTerm+parSigma);
		untweakedActinfilPerpGamma = (4*PiViscosityLength)/(naturalLogTerm+perpSigma);
		
		reTweakViscosities();
	}
	
	/** Zero out all of the instance variables that assume new values during each timestep. */
	public void reset() {
		super.reset();
		forceSum.zero();
		torque = 0;
		totalDrag = 0;
		brownianTorque = 0;
		rotK = 0;
		transParK = 0;
		transOrthoK = 0;
		brownianForcePerpendicular.zero();
		brownianForceParallel.zero();
		brownianForceRotational.zero();
		deltaTheta = 0;
		angularVel = 0;
	}
	
	public void zeroForces() {
		forceSum.zero();
		torque = 0;
		totalDrag = 0;
	}
	
	public void jump() {
		cm.set(Math.random()*Sim2D.xDimension,Math.random()*Sim2D.yDimension);
		double ang = Math.random()*2*Math.PI;
		uVect.x = Math.cos (ang);
		uVect.y = Math.sin (ang);
		uVect.uVect();
		upVect.orthogonalVector(uVect);
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
	}
	
// DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING *****
	
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		switch (segCt) {
		case 1:
			Point2D.drawLineTwixtPoints(segs[0][0],segs[0][1],g);
			break;
		case 2:
			Point2D.drawLineTwixtPoints(segs[0][0],segs[0][1],g);
			Point2D.drawLineTwixtPoints(segs[1][0],segs[1][1],g);
			break;
		case 3:
			Point2D.drawLineTwixtPoints(segs[0][0],segs[0][1],g);
			Point2D.drawLineTwixtPoints(segs[1][0],segs[1][1],g);
			Point2D.drawLineTwixtPoints(segs[2][0],segs[2][1],g);
			break;
		}
	}
	
	public void segment() {
		segs = Point2D.segmentLine(end1,end2,uVect);
		segCt = segs.length;
		for(int i = 0; i < 3; i++) {
			if(i < segCt) {
				segLength[i] = Point2D.getDistance(segs[i][0], segs[i][1]);
			}
			else segLength[i] = 0;
		}
	}

	public void wrapCoordinates() {
		cm.wrapPoint();
		end1.wrapPoint();
		end2.wrapPoint();
	}
		
			
}



