package main;

/*
	A barrier other things cannot cross
*/

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Random;

import util.*;
import collision.*;
import gui.*;

/** An barrier elements.  */
public class BarrierElement extends Thing {
	
	/** The list of all barrier elements. */
	static public BarrierElement [] theBarrierElements = new BarrierElement [50000];
		
	/** Current number of barrier elements. */
	static public int barrierElementCt = 0;

	/** Color for drawing. */
	Color barrierColor = Color.green;
	
	/** radius of a barrier. */
	static double radius = 5.0;
	
	/** minimum reflecting distance */
	static double minReflectingDistance = 8;
	
	/** The position of this barrier element in the array of all barrier elements. */
	int myBarrierElementNumber;

	/** True if this barrier element does not move in response to force. */
	boolean isPinned = false;

	/** length from end-to-end */
	public double barrierLength;
	public static double defaultBarrierLength = 100; // nm
			
	/** Unit vector parallel to long axis of barrier element. */
	Point2D uVect = new Point2D ();
	
	/** Unit vector normal to long axis of barrier element */
	Point2D upVect = new Point2D ();

	/** Position of one end (a) of barrier element. */
	public Point2D aEnd = new Point2D();
	
	/** Position of one end (b) of barrier element. */
	public Point2D bEnd = new Point2D();
		


	/** PERIODIC BOUNDARIES: Number of segments this filament is divided into. */
	int segCt = 0;
	
	/** segment endpoints*/
	Point2D [][] segs;
	
	/** PERIODIC BOUNDARIES: segment lengths. */
	double [] segLength = new double[3];
	
	/**********  CALCULATIONS  ***********/

// ??? don't know if I'll need this ??? doubt it, but saving it for now
	/* Reused for torque calculations. */
	Point2D torqArm = new Point2D();

// ??? don't know if I'll need this ??? doubt it, but saving it for now
	/* Reused in translational force calculations. */
	Point2D aTobVec = new Point2D();

// ??? don't know if I'll need this ??? doubt it, but saving it for now
	/* drag coefficients, needed for calculating diffusion */
	double barrierElementTransGamma = 0;
	double barrierElementPerpGamma = 0;
	double rotationalGamma = 0;

// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* For calculating diffusion */
	double diffusionParallel = 0;
	double diffusionPerpendicular = 0;
	double diffusionRotational = 0;
	
	double diffusionStepParallel = 0;
	double diffusionStepPerpendicular = 0;
	double diffusionStepRotational = 0;

// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* Random number generator to generate diffusive steps. */
	Random generator = new Random(Sim2D.generator.nextLong());
	
// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* For calculating forces due to diffusion */
	Point2D brownianForceParallel=new Point2D();
	Point2D brownianForcePerpendicular=new Point2D();
	Point2D brownianForceRotational=new Point2D();
	
// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* For summing forces */
	Point2D forceSum=new Point2D();

// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* Velocity of center of mass. */
	Point2D velocity=new Point2D();						// has x and y component

// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* Accumulates torque on this filament each timestep. */
	double torque = 0;

// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* Brownian torque at each timestep. */
	double brownianTorque = 0;
	
// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* change in angle this timestep. */
	double deltaTheta = 0;
	
// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* Current angle. */
	double thetaNOW = 0;
	
// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* Next angle */
	double thetaNEXT = 0;
	
// ??? don't know if I'll need this ??? doubt it, but saving it for now (might be good to have later)
	/* angular velocity this timestep */
	double angularVel = 0;

	Color myColor;
	
	static public void parametersChanged() {
	}

	// CONSTRUCTORS ***** CONSTRUCTORS ***** CONSTRUCTORS ***** CONSTRUCTORS***** CONSTRUCTORS

	/** Make a new barrier element at (0,0), angle = 0, length = 100, not Pinned */
	public BarrierElement () {
		this(0,0,0,defaultBarrierLength,false);
		barrierLength = defaultBarrierLength;
	}

	/** Make a new barrier element at (0,0), angle = 0, length = specified length, not Pinned */
	public BarrierElement (double length) {
		this(0,0,0,length,false);
		barrierLength = length;
	}
	
	/**
	 * Constructor
	 *
	 * Make a new barrier element from (x1,y1) to (x2,y2), pinned in place
	 *
	 * @param    XONE	               a  double
	 * @param    YONE	               a  double
	 * @param    XTWO	               a  double
	 * @param    YTWO	               a  double
	 */
	public BarrierElement (double xOne, double yOne, double xTwo, double yTwo) {
		this(((xOne + xTwo)/2), ((yOne + yTwo)/2), (Math.asin((yTwo-yOne)/(Math.pow (((xTwo-xOne)*(xTwo-xOne)+(yTwo-yOne)*(yTwo-yOne)), 0.5)))),(Math.pow (((xTwo-xOne)*(xTwo-xOne)+(yTwo-yOne)*(yTwo-yOne)), 0.5)),false);
		barrierLength = Math.pow (((xTwo-xOne)*(xTwo-xOne)+(yTwo-yOne)*(yTwo-yOne)), 0.5);
	}
	
	
	/**
	 * Constructor
	 *
	 * Make a new barrier element at (initX,initY), angle = 0, length = length, not Pinned
	 *
	 * @param    initX               a  double
	 * @param    initY               a  double
	 * @param	 barrierLength		 a  double
	 */
	public BarrierElement (double initX, double initY, double length) {
		this(initX,initY,0,length,false);
		barrierLength = length;
	}
	
//  using another construct with inputs (double, double, double, double) for barrier between two pts
//	/**
//	 * Constructor
//	 * Make a new BarrierElement that is not pinned
//	 *
//	 * @param    initX               a  double  = initial X position
//	 * @param    initY               a  double = initial Y position
//	 * @param    initAng             a  double = initial angle
//	 * @param    barrierLength       a double = length in nm
//	 *
//	 */
//
//	public BarrierElement (double initX, double initY, double initAng) {
//		this(initX,initY,initAng,length,false);
//		barrierLength = length;
//	}
	
	/**
	 * Constructor
	 *
	 * @param    initX               a  double  = initial X position
	 * @param    initY               a  double = initial Y position
	 * @param    initAng             a  double = initial angle
	 * @param    barrierLength       a double = length in nm
	 * @param    pinned              a  boolean = does this filament move?
	 *
	 */
	public BarrierElement (double initX, double initY, double initAng, double length, boolean pinned) {
		super (initX,initY);
		this.isPinned = pinned;
		defineColor();
		uVect.set(Math.cos(initAng), Math.sin(initAng));
		uVect.uVect();

		upVect.orthogonalVector(uVect);
					
		barrierLength = length;
		useCenterMassAndUVectToDetermineEnds();
		segmentBarrierElement();
		wrapBarrierCoordinates();
		addBarrierElement(this);
		evaluateProperties();
	}
	
	
	// INTEGRATION STEP ***** INTEGRATION STEP ***** INTEGRATION STEP ***** INTEGRATION STEP *****
	
	/** Looop through all filaments and call each one's step() method. */
	static public void stepAllBarrierElements(double dT)  {

		for (int i = 0; i < barrierElementCt; i++) {
			if (!theBarrierElements[i].removeMe) { theBarrierElements[i].step(dT); }
		}
	}
		
	/** Register all filaments with the Collision detectors actin mesh. */
	static public void meshBarrierElement()  {
//		System.out.println("BarrierElement.meshBarrierElement: method started...");
		for (int i = 0; i < barrierElementCt; i++) {
			theBarrierElements[i].fillBarrierMesh();
//			System.out.println("BarrierElement.meshBarrierElement: filling mesh for barrier "+i);
		}
//		System.out.println("BarrierElement.meshBarrierElement: method complete.");
	}

	private void fillBarrierMesh() {
//		System.out.println("Mesh.fillBarrierMesh: segment method started... segCt = "+barrier.segCt);
		switch (segCt) {
		case 1:
//			System.out.println("Mesh.fillBarrierMesh: case 1");
			Mesh.BARRIER_MESH.addLineSegmentToMesh(myBarrierElementNumber,segs[0][0],segs[0][1]);
			break;
		case 2:
//			System.out.println("Mesh.fillBarrierMesh: case 2");
			Mesh.BARRIER_MESH.addLineSegmentToMesh(myBarrierElementNumber,segs[0][0],segs[0][1]);
			Mesh.BARRIER_MESH.addLineSegmentToMesh(myBarrierElementNumber,segs[1][0],segs[1][1]);
			break;
		case 3:
			System.out.println("Mesh.fillBarrierMesh: case 3");
			Mesh.BARRIER_MESH.addLineSegmentToMesh(myBarrierElementNumber,segs[0][0],segs[0][1]);
			Mesh.BARRIER_MESH.addLineSegmentToMesh(myBarrierElementNumber,segs[1][0],segs[1][1]);
			Mesh.BARRIER_MESH.addLineSegmentToMesh(myBarrierElementNumber,segs[2][0],segs[2][1]);
			break;
		}
//		System.out.println("Mesh.fillBarrierMesh: segment method complete.");
	}

/** Determine the positions of the ends (a,b) from center of mass and orientation (µVect) */
	public void useCenterMassAndUVectToDetermineEnds(){
		aEnd.add(cm,-barrierLength/2,uVect);
		bEnd.add(cm,barrierLength/2,uVect);
	}
	
	/** Set the new center of mass from the positions of ends (a b). */
	public void setNewCM(){
		aTobVec.getVector(aEnd, cm);
		aTobVec.uVect();
		cm.add(aEnd,barrierLength/2,aTobVec);
		cm.wrapPoint();
	}
	
	static public double getTotalBarrierElementLength() {
		double l = 0;
		for(int i = 0; i < barrierElementCt; i++) {
			l+= theBarrierElements[i].barrierLength;
		}
		return l;
	}
				
	
	/** Melt this entire barrier element, targeting it for destruction.
	 */
	public void meltBarrierElement(){
		setForRemove();
	}
	
	

	/** Calculate the brwonian forces and moments on this barrier element, add them to the forces accumulated due to
	 * collisions(?), then compute and take a translational and rotational step.
	 * Note that this method assumes that all other forces acting on this barrier element have already been
	 * computed and added to the barrier element.
	 */
	public void step(double dT) {
		
		if(!isPinned) {
			
			//***********************************************************************************************
			// PARALLEL DIFFUSION FORCES - moving in the direction of uVect
			// you are going to diffuse with this step size in the direction of uVect
			diffusionStepParallel = Math.sqrt(2*diffusionParallel*dT)* generator.nextGaussian();
			
			//Convert the step in parallel direction to x,y components of forces in the universal coordinate system
			brownianForceParallel.x = uVect.x * ((barrierElementTransGamma / dT) * diffusionStepParallel);
			brownianForceParallel.y = uVect.y * ((barrierElementTransGamma / dT) * diffusionStepParallel);
			
			
			
			//***********************************************************************************************
			//PERPENDICULAR DIFFUSION FORCES - moving in the direction of upVect
			diffusionStepPerpendicular = Math.sqrt(2*diffusionPerpendicular*dT)* generator.nextGaussian();
			
			brownianForcePerpendicular.x =  upVect.x * ((barrierElementPerpGamma / dT) * diffusionStepPerpendicular);
			brownianForcePerpendicular.y =  upVect.y * ((barrierElementPerpGamma / dT) * diffusionStepPerpendicular);
			
			
			
			//TOTAL TRANSLATIONAL FORCES**********************************************************************
			
			//calculate the total force acting on the center of mass
			//System.out.println("transXforce before Brownian and constant = " + totalTranslationalForce.x);
			forceSum.x += (brownianForceParallel.x + brownianForcePerpendicular.x);
			forceSum.y += (brownianForceParallel.y + brownianForcePerpendicular.y);
			
			velocity.x = forceSum.x / barrierElementTransGamma;
			velocity.y = forceSum.y / barrierElementTransGamma;
			
			// Move the center of mass.
			cm.inc(dT,velocity);
			if (!Point2D.pointOK(cm)) {
				System.out.println ("barrier element cm is NaN: barrier element = " + this + " v.x = " + velocity.x + " v.y = " + velocity.y);
			}
			
			//ROTATIONAL DIFFUSION FORCES**********************************************************************
			diffusionStepRotational = Math.sqrt(2*diffusionRotational*dT)* generator.nextGaussian();
			
			double brownianTorque = rotationalGamma * diffusionStepRotational / dT;
			
			torque +=  brownianTorque;
			
			angularVel = torque / rotationalGamma;
			
			
			//TOTAL ANGULAR FORCES**************************************************************************
	
			deltaTheta = angularVel * dT;		//contribution to change in theta from diffusion
			
			thetaNOW = Math.atan2 (uVect.y, uVect.x);
			thetaNEXT = thetaNOW +deltaTheta;
	
			uVect.x = Math.cos (thetaNEXT);
			uVect.y = Math.sin (thetaNEXT);
			uVect.uVect();
			upVect.orthogonalVector(uVect);
		}
				
		// Finally, calculate the position of the ends
		useCenterMassAndUVectToDetermineEnds();
		segmentBarrierElement();
		wrapBarrierCoordinates();
		resetCounters();
	}
	
	/** Zero out all of the instance variables that assume new values during each timestep. */
	public void resetCounters() {
		forceSum.zero();
		torque = 0;
		brownianTorque = 0;
		brownianForcePerpendicular.zero();
		brownianForceParallel.zero();
		brownianForceRotational.zero();
		velocity.zero();
		deltaTheta = 0;
		angularVel = 0;
	}
	
	
	
	/** Add the force to the barrier element's center of mass, thus no torque. */
	public void addForce(Point2D f){
		forceSum.inc(f);
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
				
	/** Compute length dependent drag coefficients. */
	public void evaluateProperties () {
		double naturalLogTerm = (double) Math.log(barrierLength/(2*radius));
		double PiViscosityLength= (double) Math.PI* Constants.viscosity*barrierLength;

		// double pValue = (2*radius)/length;
		// the Sigmas would normally be looked up from a table given the pValue.
		// Jacked from page 107 of Jonathan Howard's Mechanics of Motor Proteins and the Cytoskeleton
		double rotSigma = -0.66;
		double parSigma = -0.20;
		double perpSigma = 0.86;
		
		double pValue = (2*radius)/barrierLength;
		boolean found=false;
		for(int i=0;i<Drag.cylinderDragValues.length&&!found;i++){
			if(pValue>=Drag.cylinderDragValues[i][Drag.pValueIndex]){
				parSigma=Drag.cylinderDragValues[i][Drag.parIndex];
				perpSigma=Drag.cylinderDragValues[i][Drag.perpIndex];
				rotSigma=Drag.cylinderDragValues[i][Drag.rotIndex];
				found=true;
			//	System.out.println("using pValue="+pValue+"  index="+i+"  parSigma="+parSigma+"  perpSigma="+perpSigma+"  rotSigma="+rotSigma);
			}
		}
		
		if(!found){
			int defaultIndex=Drag.cylinderDragValues.length-1;
			parSigma=Drag.cylinderDragValues[defaultIndex][Drag.parIndex];
			perpSigma=Drag.cylinderDragValues[defaultIndex][Drag.perpIndex];
			rotSigma=Drag.cylinderDragValues[defaultIndex][Drag.rotIndex];
			//System.out.println("pValue outside of range so using default");
		}
		
		rotationalGamma = (PiViscosityLength*barrierLength*barrierLength)/(3*(naturalLogTerm+rotSigma));
		barrierElementTransGamma = (2*PiViscosityLength)/(naturalLogTerm+parSigma);
		barrierElementPerpGamma = (4*PiViscosityLength)/(naturalLogTerm+perpSigma);
		
		diffusionParallel = Constants.kT / barrierElementTransGamma;
		diffusionPerpendicular = Constants.kT / barrierElementPerpGamma;
		diffusionRotational =  Constants.kT / rotationalGamma;
	}

	
// DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING *****
	public void defineColor(){
		myColor=new Color(255,0,0);
		double value=Math.random();
		if(value<0.15){
			myColor=new Color(255,0,0);
		}
		else if(value<0.25){
			myColor=new Color(225,0,0);
		}
		else if(value<0.5){
			myColor=new Color(200,0,0);
		}
		else if(value<0.75){
			myColor=new Color(175,0,0);
		}
		else{
			myColor=new Color(150,0,0);
		}
		//System.out.println("MyColor="+myColor.getRed()+","+myColor.getGreen()+","+myColor.getBlue());
	}
	
	
	double angle=0;
	GeneralPath perimeter = new GeneralPath(); // create GeneralPath object
	boolean showBox=false;
	boolean showEnds=true;
	boolean showEndStrings=false;
	boolean showCenterMass=false;
	boolean showLine=true;
	boolean showCollisionStuff=false;
	Point2D tempPoint2DP = new Point2D();
	Point2D tempPoint2DB = new Point2D();
	Point2D tempPoint2Dp1I = new Point2D();
	Point2D tempPoint2Dp2I = new Point2D();
	
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		//if (this instanceof GridActin) { System.out.println ("gonna draw a grid actin that is " + this.length + " nm long centered at " + this.centermass.reportVals()); }
		g.setColor(barrierColor);
		
		int XaEnd = aEnd.getPixX();
		int YaEnd = aEnd.getPixY();
		int XbEnd = bEnd.getPixX();
		int YbEnd = bEnd.getPixY();
		
		// BELOW IS CRAP FOR TESTING THE COLLISON DETECTOR
		if (showCollisionStuff) {
			g.setColor(Color.yellow);
			g.drawOval(tempPoint2DP.getPixX(),tempPoint2DP.getPixY(),5,10);
			g.drawOval(tempPoint2DB.getPixX(),tempPoint2DB.getPixY(),5,10);
			g.drawOval(tempPoint2Dp1I.getPixX(),tempPoint2Dp1I.getPixY(),5,10);
			g.drawOval(tempPoint2Dp2I.getPixX(),tempPoint2Dp2I.getPixY(),5,10);
		}
		// ABOVE IS CRAP FOR TESTING THE COLLISON DETECTOR
		
		int xCenter;
		int yCenter;
		int pixelRadius;
		int pixelDiameter;

		if (showBox) {
			if (segCt >= 1) { perimeterRender(segs[0][0], segs[0][1], g, scale, offset); }
			if (segCt >= 2) { perimeterRender(segs[1][0], segs[1][1], g, scale, offset); }
		
		}
		
		if(showEnds){
			int EndPixelRadius = (int) (3*radius*scale);
			int EndPixelDiameter = (int) (6*radius*scale);
			
			g.setColor(Color.WHITE);
			g.fillOval(XaEnd-EndPixelRadius, YaEnd-EndPixelRadius, EndPixelDiameter, EndPixelDiameter);
			
			g.setColor(Color.green);
			g.fillOval(XbEnd-EndPixelRadius, YbEnd-EndPixelRadius, EndPixelDiameter, EndPixelDiameter);
		}
		
		if(showEndStrings){
			int EndPixelRadius = (int) (radius*scale);
			int EndPixelDiameter = (int) (2*radius*scale);
			int padding=20;
			g.setColor(Color.WHITE);
			g.drawString("A_END",XaEnd, YaEnd+padding);
			
			g.setColor(Color.green);
			g.drawString("B_END",XbEnd, YbEnd+padding);
		}
		
		
		if(showCenterMass){
			g.setColor(Color.MAGENTA);
			xCenter = cm.getPixX();
			yCenter = cm.getPixY();
			pixelRadius = (int) (radius*scale);
			pixelDiameter = (int) (2*radius*scale);
			
			g.fillOval(xCenter-pixelRadius,yCenter-pixelRadius,pixelDiameter,pixelDiameter);
		}
		
		if(showLine){
			g.setColor(barrierColor);
			switch (segCt) {
			case 1:
				drawLineTwixtPoint2Ds(segs[0][0],segs[0][1],g);
				break;
			case 2:
				drawLineTwixtPoint2Ds(segs[0][0],segs[0][1],g);
				drawLineTwixtPoint2Ds(segs[1][0],segs[1][1],g);
				break;
			case 3:
				drawLineTwixtPoint2Ds(segs[0][0],segs[0][1],g);
				drawLineTwixtPoint2Ds(segs[1][0],segs[1][1],g);
				drawLineTwixtPoint2Ds(segs[2][0],segs[2][1],g);
				break;
			}
		}
	}
	
	public void drawLineTwixtPoint2Ds (Point2D p1, Point2D p2, Graphics g) {
		g.drawLine(p1.getPixX(), p1.getPixY(), p2.getPixX(), p2.getPixY());
	}
	
	public void perimeterRender (Point2D end1, Point2D end2, Graphics g, double scale, double [] offset) {
		Graphics2D g2d=(Graphics2D)g;
		double height=radius;
		double tempxT1,tempxT2,tempyT1,tempyT2;
		double tempxB1,tempxB2,tempyB1,tempyB2;
		double ux,uy,pux,puy;
		ux=uVect.x;
		uy=uVect.y;
		pux=upVect.x;
		puy=upVect.y;
		tempxT1=((end1.x-height*pux)-offset[0])*scale;
		tempyT1=((end1.y-height*puy)-offset[0])*scale;
		tempxT2=((end2.x-height*pux)-offset[0])*scale;
		tempyT2=((end2.y-height*puy)-offset[0])*scale;
		
		tempxB1=((end1.x+height*pux)-offset[0])*scale;
		tempyB1=((end1.y+height*puy)-offset[0])*scale;
		tempxB2=((end2.x+height*pux)-offset[0])*scale;
		tempyB2=((end2.y+height*puy)-offset[0])*scale;
		
		g2d.setPaint(Color.RED);
	
		g2d.setPaint(myColor);
		perimeter.reset();
		perimeter.moveTo((float)tempxT1,(float)tempyT1);
		perimeter.lineTo((float)tempxT2,(float)tempyT2);
		perimeter.lineTo((float)tempxB2,(float)tempyB2);
		perimeter.lineTo((float)tempxB1,(float)tempyB1);
		perimeter.closePath();
		
		g2d.fill(perimeter);
	}
	
// END OF DRAWING ***** END OF DRAWING ***** END OF DRAWING ***** END OF DRAWING ***** END OF DRAWING *****
			
	public static void addBarrierElement (BarrierElement newBarrierElement) {
		newBarrierElement.myBarrierElementNumber = barrierElementCt;
		theBarrierElements[barrierElementCt] = newBarrierElement;
		barrierElementCt ++;
	}

	public void die() {
	}
	
		
	
	public void removeMe () {
		super.removeMe();
		int swapId = myBarrierElementNumber;
		theBarrierElements[swapId] = theBarrierElements[barrierElementCt-1];
		theBarrierElements[swapId].myBarrierElementNumber = swapId;
		barrierElementCt --;
	}
	
	public void segmentBarrierElement() {
		segs = Point2D.segmentLine(aEnd,bEnd,uVect);
		segCt = segs.length;
		for(int i = 0; i < 3; i++) {
			if(i < segCt) {
				segLength[i] = Point2D.getDistance(segs[i][0], segs[i][1]);
			}
			else segLength[i] = 0;
		}
		wrapBarrierCoordinates();
	}

	private void wrapBarrierCoordinates() {
		cm.wrapPoint();
		aEnd.wrapPoint();
		bEnd.wrapPoint();
	}
		

	
	public void printMe() {
		System.out.println(this+"myBarrierElementNumber="+myBarrierElementNumber+"  thingNumber="+myThingNumber+"  simTime="+Sim2D.simulationTime+"  barrierLength="+barrierLength +"  removeMe="+removeMe);
	}
	
	static public void writeHistograms() {
	}
	
	public static void printAll () {
		System.out.println("ALL BARRIER ELEMENT PRINT");
		for (int i=0;i<BarrierElement.barrierElementCt;i++) {
			BarrierElement.theBarrierElements[i].printMe();
		}
		System.out.println("END ALL BARRIER ELEMENT PRINT");
	}

	static public void staticInit() {
		
		Sim2D.addProxy(new AgentProxy() {
			public void registerForCollisionDetection() {
				meshBarrierElement();
			}
			public void step(double dT) {
				// Barrier elements move and update positions (moot for now).
				stepAllBarrierElements(dT);
			}
			public void reset() {
				for(int i = 0; i < barrierElementCt; i++) {
					theBarrierElements[i].die();
					theBarrierElements[i] = null;
				}
				barrierElementCt = 0;
			}
		});
	}

}

