package main;

/**
/*	A rectangular array of myosin V tails (coiled coil/oligomerization domains).
	 */

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

import util.*;
import collision.*;
import gui.*;
import parameters.*;

public class MyosinVTail extends Thing {

	static CollisionDetector myCollisionDetector = new MyosinVTailCollisionDetector();
	
	/** array holding all MyosinVTails */
	static public MyosinVTail[] theMyosinVTails = new MyosinVTail[10000];
	
	/** Current number of myosinVtails. */
	static public int myosinVTailCt = 0; // initialize the LongMyosinVTail count
	
//	/** Color for drawing. */
//	Color MyosinVTailColor = Color.green;
//	MyosinVTailColor = 0;
//	setPaintColor();
//	public Color myColor;
	public Color myosinVTailColor = Color.blue;
	
	/** color for jumping. */
	Color jumpingColor = Color.DARK_GRAY;

	/** Unique index of this myosinVtail in the array of all myosinVtails. */
	public int myMyosinVTailNumber;
	/** indicates if this myosinVtail has at least 1 head bound to actin */
//	public boolean boundToActin;
	public int headOneState;
	public int headTwoState;
	
	public int headOnePosition;
	public int headTwoPosition;
	
	public int headOneMonomer;
	public int headTwoMonomer;
	
	public int headOneActinFil;
	public int headTwoActinFil;

	/** define states for myosin V tail, part of a run or not part of a run... */
	static public int
		NOT_RUNNING = 0,
		RUNNING = 1;

	public int runStatus = NOT_RUNNING;
	
	/** descriptors for individual myosin V head positions */
	static public int
		UNBOUND = 0,
		SINGLE = 1,
		LEADING = 2,
		TRAILING = 3,
		INTERSECTION = 4;

	/** histogram of time spent with one or both heads bound, runlengths, velocities */
	static public HistogramPlus runTimeHistogram;
	static public HistogramPlus runlengthHistogram;
	static public HistogramPlus runVelocityHistogram;

	/** run start time and position */
	public double runStartTime;
	public double runStartX;
	public double runStartY;
	
	/** run current time and positions */
	public double runEndTime;
	public double runEndX;
	public double runEndY;
	
	public double runlength;
	public double runtime;
	public double runvelocity;
			
	static public int runCount = 0;
	public int thisRun;

	static public double jCutoff = 300;  // distance cutoff for storing runlengths and velocities
	static public double minRunTime = 0.1;  // time cutoff for storing runlengths and velocities
	        
	/** mean runlength and velocity */
	static public double meanRunlength;
	static public double meanRunVelocity;
	static public double sdRunlength;
	static public double sdRunVelocity;
	static public int nRunlengths;
	static public int nRunVelocities;
	static public double seRunlength;
	
	/** vectors to gather distributions for runlength and velocity */
	static public Vector runlengthsList = new Vector();
	static public Vector runVelocitiesList = new Vector();
	static public int runlengthCt = 0;

	static public double
		meanRunlengthList,
		sdRunlengthList,
		meanRunVelocitiesList,
		sdRunVelocitiesList,
		seRunlengthList;
	
	/** indicates if this myosinVtail spans one or more boundaries */
	boolean wrapped = false;
	
	/** increase in diffusivity out of simulation plane (i.e. off cortex) */
	static double outOfPlaneDiffFactor = 1000;

	/** Out of plane diffusivity. */
	double outOfPlaneDiff;

	/** Z position of A end. */
	int zPosA = 0;

	/** Z position of B end. */
	int zPosB = 0;

	/** Maximum Z position. */
	static int maxZPos = 8;

	/** Probabaility that the myosinVtail will take a step in Z. */
	static double zPosDeltaProb = 0.0;

	/** Proabability that a myosinVtail will release the cortex altogether. */
	static double rdmCortexReleaseProb = 0.0;

	/** increase in diffusivity for myosinVtail not bound to cortex. */
	static double offCortexDiffMult = 1000;

	/** Length of a myosinVtail from end to end. */
	static double length;

	/** Half the length of a myosinVtail rod. */
	static double halfLength;

	/** Width of the myosinVtail. */
	static double width;

	/** Half the width of the myosinVtail. */
	static double halfWidth;

	/**
	 * square of collision radius used for myosinVtail collision calculations.
	 */
	static public double collisionRadSqrd;

	/** number of myosinV heads on each end of tail. */
	static int nMyoVHeadsA;
	static int nMyoVHeadsB;
	
	/** Unit vector along myosinVtail axis. */
	Point2D uVect = new Point2D();

	/** Unit vector perpendicular to myosinVtail axis. */
	Point2D upVect = new Point2D();

	/** Reused for torque calculations. */
	Point2D torqArm = new Point2D();

	/** Reused for collision caclulations. */
	public Point2D tempP = new Point2D();

	/*
	 * First end Point2D of the myosinVtail, calculated based on length, centermass,
	 * and uVect
	 */
	Point2D end1 = new Point2D();

	/*
	 * Second end Point2D of the myosinVtail, calculated based on length, centermass,
	 * and uVect
	 */
	Point2D end2 = new Point2D();

	Point2D A1 = new Point2D();
	Point2D A2 = new Point2D();
	Point2D B1 = new Point2D();
	Point2D B2 = new Point2D();
	
//	Point2D lastMVTailPosition = new Point2D();
//	Point2D nextMVTailPosition = new Point2D();
	public double lastXpos;
	public double lastYpos;
	public double nextXpos;
	public double nextYpos;
	
	/** Array of myosinV heads attached to "A" end of myosinVtail. */
	MyosinV[] myAmyosinVs;

	/** Array of myosinV heads attached to "B" end of myosinVtail. */
	MyosinV[] myBmyosinVs;

	/** Array of all myosinV heads. */
	MyosinV[] myMyosinVs;

	/** drag coefficients, needed for calculating diffusion */
	public double myosinVtailTransGamma, myosinVtailPerpGamma, rotGam;

	/** For calculating diffusion. */
	double diffPar, diffPerp, diffRot;
	double curDPar, curDPerp, curDRot; // normally equal to diffPar, etc, but
										// used in step() to increase off cortex
										// diffusivity
	double diffStepPar, diffStepPerp, diffStepRot;

	Random generator = new Random(Sim2D.generator.nextLong());
	Random jumpGen = new Random(Sim2D.generator.nextLong());

	/** For calculating forces due to diffusion. */
	Point2D brownianForceParallel = new Point2D();
	Point2D brownianForcePerpendicular = new Point2D();
	Point2D brownianForceRotational = new Point2D();

	/** For summing forces. */
	Point2D totalTranslationalForce = new Point2D();
	Point2D velocity = new Point2D();

	double torque = 0;
	double brownianTorque = 0;
	double deltaTheta = 0;
	double thetaNOW = 0;
	double thetaNEXT = 0;
	double angularVel = 0;

	boolean myosinVOnLongAxis = false;

	// MOSHE DIAGNOSTIC
	// for determining if myosinV heads in myosinVtail are active or "dead"
	boolean myoVTailstate = MyosinV.initMyoVTailstate;

	/** for rendering */
	boolean showInfo = false;
	boolean showEndsInfo = false;
	boolean showUVect = true;
	boolean fileMyoVTailOffCortex;
	GeneralPath perimeter = new GeneralPath();


	// CONSTRUCTOR
	// *************************************************************************************************

	public MyosinVTail(double initX, double initY, double uVectX,
			double uVectY) {
		super(initX, initY);

		setMyosinVTailParams();

		end1.add(cm, -halfLength, uVect);
		end2.add(cm, halfLength, uVect);

		uVect.set(uVectX, uVectY);
		uVect.uVect();
		upVect.orthogonalVector(uVect);

		myosinVTailColor = Color.blue;
		
		createMyosinVHeads(nMyoVHeadsA, nMyoVHeadsB);
		addMyosinVTail(this);
		evaluateProperties();
	}

	public void createMyosinVHeads(int a, int b) {
		zPosA = 0;// (int)(Math.random()*maxZPos);
		zPosB = 0;// (int)(Math.random()*maxZPos);

		a = 2 * ((a + 1) / 2); // n should always be even and there should be at
								// least 2 heads at either end
		b = 2 * ((b + 1) / 2); // n should always be even and there should be at
		// least 2 heads at either end

		myAmyosinVs = new MyosinV[a];
		myBmyosinVs = new MyosinV[b];
		myMyosinVs = new MyosinV[myAmyosinVs.length + myBmyosinVs.length];

		for (int i = 0; i < myAmyosinVs.length; i++) {
			myAmyosinVs[i] = new MyosinV(this);
			myMyosinVs[i] = myAmyosinVs[i];
			myAmyosinVs[i].orientation = 1;
			myAmyosinVs[i].iAmAnAMyoV = true;
		}
		for (int i = 0; i < myBmyosinVs.length; i++) {
			myBmyosinVs[i] = new MyosinV(this);
			myMyosinVs[i + myAmyosinVs.length] = myBmyosinVs[i];
			myBmyosinVs[i].orientation = -1;
			myBmyosinVs[i].iAmAnAMyoV = false;
		}

		// A's are the top row, and B's are the bottom row
		double spacing = MyosinV.myoVRadius;
		Point2D temp = new Point2D();
		for (int i = 0; i < myAmyosinVs.length; i += 2) {
			temp.x = halfLength + i * spacing;
			temp.y = halfWidth;
			myAmyosinVs[i].setInitialPositionRelativeToMyosinVTail(temp);
			temp.y = -halfWidth;
			myAmyosinVs[i + 1].setInitialPositionRelativeToMyosinVTail(temp);
		}
		for (int i = 0; i < myBmyosinVs.length; i += 2) {
			temp.x = -halfLength - i * spacing;
			temp.y = halfWidth;
			myBmyosinVs[i].setInitialPositionRelativeToMyosinVTail(temp);
			temp.y = -halfWidth;
			myBmyosinVs[i + 1].setInitialPositionRelativeToMyosinVTail(temp);
		}

		moveMyosinVs();
	}

	public MyosinVTail() {
		this(0, 0, 1, 0);
	}

	public MyosinVTail(double initX, double initY) {
		this(initX, initY, 1, 0);
	}

	static public void fillMVTailMesh() {
		MyosinVTail m;
		for (int i = 0; i < myosinVTailCt; i++) {
			m = theMyosinVTails[i];
			Mesh.MVTAIL_MESH.addPointToMesh(m.myMyosinVTailNumber,m.cm,MyosinV.myoVCollisionRadius);
		}
	}

	static public void stepAllMyosinVtails(double dT) {
		for (int i = 0; i < myosinVTailCt; i++) {
			if (!theMyosinVTails[i].removeMe) {
				theMyosinVTails[i].step(dT);
			}
		}
	}

	// createMyMonomers (numMonomers-2);
	// length = monLength*numMonomers;
	// addLongMiniF(this);
	//
	// uVect.Set(uVectX, uVectY);
	// uVect.uVect();

	public void addForce(Point2D f) {
		totalTranslationalForce.inc(f);
	}

	public void addForce(Point2D f, Point2D loc) {
		addForce(f);
		torqArm.getVector(cm, loc);
		addTorque(torqArm, f);
	}

	public void addTorque(Point2D r, Point2D f) {
		torque += Point2D.Cross(r, f);
	}

	// *************************************************************************************************

	public void step(double dT) {
		
		
//		PositionInfo info = new PositionInfo();
		assignLastMVTailPosition(cm);
		
		setMyosinVTailParams();

		zPosCalculations(dT);

		// System.out.println ("zPosA,zPosB = " + zPosA + "," + zPosB);

		// ***********************************************************************************************
		// PARALLEL DIFFUSION FORCES - moving in the direction of uVect
		// you are going to diffuse with this step size in the direction of
		// uVect
		diffStepPar = Math.sqrt(2 * curDPar * dT)
				* generator.nextGaussian();

		// Convert the step in parallel direction to x,y components of forces in
		// the universal coordinate system
		brownianForceParallel.set((myosinVtailTransGamma * diffStepPar)
				/ dT, uVect);

		// ***********************************************************************************************
		// PERPENDICULAR DIFFUSION FORCES - moving in the direction of upVect
		diffStepPerp = Math.sqrt(2 * curDPerp * dT)
				* generator.nextGaussian();

		brownianForcePerpendicular.set((myosinVtailPerpGamma * diffStepPerp)
				/ dT, upVect);

		// ***********************************************************************************************
		// TOTAL TRANSLATIONAL FORCES
		// calculate the total force acting on the center of mass

		totalTranslationalForce.inc(brownianForceParallel);
		totalTranslationalForce.inc(brownianForcePerpendicular);

		velocity.set(1 / myosinVtailTransGamma, totalTranslationalForce); // WHAT
																		// DRAG
																		// COEFFICIENT
																		// TO
																		// USE?

		// Move the center of mass.

		cm.inc(dT, velocity);
		

		
		// ************************************************************************************************
		// ROTATIONAL DIFFUSION FORCES
		diffStepRot = Math.sqrt(2 * curDRot * dT)
				* generator.nextGaussian();

		brownianTorque = rotGam * diffStepRot / dT;

		torque += brownianTorque;

		angularVel = torque / rotGam;

		// ************************************************************************************************
		// TOTAL ANGULAR FORCES
		deltaTheta = angularVel * dT; // contribution to change in
												// theta from diffusion

		thetaNOW = Math.atan2(uVect.y, uVect.x);
		// deltaTheta=Math.PI*0.0005;
		thetaNEXT = thetaNOW + deltaTheta;

		uVect.x = Math.cos(thetaNEXT);
		uVect.y = Math.sin(thetaNEXT);

		uVect.uVect();

		// also set the perp. unit vector
		upVect.orthogonalVector(uVect);

		// calculate the position of the ends
		end1.add(cm, -halfLength, uVect);
		end2.add(cm, halfLength, uVect);

		// wrap all positions
		cm.wrapPoint();
		boolean end1Wrapped = end1.wrapPoint();
		boolean end2Wrapped = end2.wrapPoint();
		if (end1Wrapped | end2Wrapped) {
			wrapped = true;
		} else {
			wrapped = false;
		}

		// Move the MyosinVs
		moveMyosinVs();

		// Random release from cortex sim
		if (Math.random() < rdmCortexReleaseProb * dT) {
			// System.out.println(this + " is now off cortex");
			releaseAll();
			zPosA = 1;
			zPosB = 1;
		}

		if (!Point2D.pointOK(cm)) {
			System.out.println("myosinVtail cm is NaN");
			System.out.println("velocity is " + velocity.reportVals());
			System.out.println("totalTranslationalForce is "
					+ totalTranslationalForce.reportVals());
			Toolkit.getDefaultToolkit().beep();
			System.exit(0);
		}

//		checkBoundToActin();
		getHeadPositions();
		if (runStatus == NOT_RUNNING && headOnePosition == LEADING) {
			initiateRun();
		}
		if (runStatus == NOT_RUNNING && headOnePosition == TRAILING) {
			initiateRun();
		}
		if (runStatus == RUNNING && headOnePosition == UNBOUND && headTwoPosition == UNBOUND) {
			endRun();
		}

//		System.out.println("mvtail " + myMyosinVTailNumber + ": last x,y = " + lastMVTailPosition.x + ", " + lastMVTailPosition.y);
		
		assignNextMVTailPosition(cm);

//		System.out.println("mvtail " + myMyosinVTailNumber + ": last x,y = "+lastXpos+", "+lastYpos+"; next x,y = " + nextXpos + ", " +nextYpos);
//		OK to this Point2D!
		reset();
	}

	public void setDiffMultipliers(double transScale, double rotScale) {
		curDPar = transScale * diffPar;
		curDPerp = transScale * diffPerp;
		curDRot = rotScale * diffRot;
	}

	public void zPosCalculations(double dT) {
		double dTransF = 1;
		double dRotF = 1;
		if (zPosA > 0 & zPosB > 0) {
			dTransF = offCortexDiffMult;
			dRotF = offCortexDiffMult;
		}
		setDiffMultipliers(dTransF, dRotF);

		if (myosFreeA()) {
			if (Math.random() < dT * zPosDeltaProb) {
				if (Math.random() < 0.5) {
					zPosA--;
					if (zPosA < 0) {
						zPosA = 0;
					}
				} else {
					zPosA++;
					if (zPosA > maxZPos) {
						zPosA = maxZPos;
					}
				}
			}
		}

		if (myosFreeB()) {
			if (Math.random() < dT * zPosDeltaProb) {
				if (Math.random() < 0.5) {
					zPosB--;
					if (zPosB < 0) {
						zPosB = 0;
					}
				} else {
					zPosB++;
					if (zPosB > maxZPos) {
						zPosB = maxZPos;
					}
				}
			}
		}
	}

	public boolean allMyosFree() {
		for (int i = 0; i < myAmyosinVs.length; i++) {
			if (!myAmyosinVs[i].isFree())
				return false;
		}
		for (int i = 0; i < myBmyosinVs.length; i++) {
			if (!myBmyosinVs[i].isFree())
				return false;
		}
		return true;
	}

	public boolean myosFreeA() {
		for (int i = 0; i < myAmyosinVs.length; i++) {
			if (!myAmyosinVs[i].isFree())
				return false;
		}
		return true;
	}

	public boolean myosFreeB() {
		for (int i = 0; i < myBmyosinVs.length; i++) {
			if (!myBmyosinVs[i].isFree())
				return false;
		}
		return true;
	}

	public boolean offCortex() {
		if (zPosA != 0 & zPosB != 0) {
			return true;
		}
		return false;
	}

	public void moveMyosinVs() {
		setMyosinVTailParams();
		Point2D temp = new Point2D();
		for (int i = 0; i < myAmyosinVs.length; i++) {
			temp.x = cm.x + myAmyosinVs[i].mfAttachmentPtLocal.x * uVect.x
					+ myAmyosinVs[i].mfAttachmentPtLocal.y * upVect.x;
			temp.y = cm.y + myAmyosinVs[i].mfAttachmentPtLocal.x * uVect.y
					+ myAmyosinVs[i].mfAttachmentPtLocal.y * upVect.y;

			// wrap
			temp.wrapPoint();
			myAmyosinVs[i].mfAttachmentPt.set(temp);
			if (myAmyosinVs[i].isFree()) {
				myAmyosinVs[i].bindingSite.set(temp);
			}
		}
		for (int i = 0; i < myBmyosinVs.length; i++) {
			temp.x = cm.x + myBmyosinVs[i].mfAttachmentPtLocal.x * uVect.x
					+ myBmyosinVs[i].mfAttachmentPtLocal.y * upVect.x;
			temp.y = cm.y + myBmyosinVs[i].mfAttachmentPtLocal.x * uVect.y
					+ myBmyosinVs[i].mfAttachmentPtLocal.y * upVect.y;

			// wrap
			temp.wrapPoint();
			myBmyosinVs[i].mfAttachmentPt.set(temp);
			if (myBmyosinVs[i].isFree()) {
				myBmyosinVs[i].bindingSite.set(temp);
			}
		}
	}

	public void reset() {
		// reset forces and velocities
		totalTranslationalForce.zero();
		torque = 0;

		brownianForcePerpendicular.zero();
		brownianForceParallel.zero();
		brownianForceRotational.zero();
		velocity.zero();

		deltaTheta = 0;
		angularVel = 0;
	}

	public void evaluateProperties() {
		setMyosinVTailParams();
		double naturalLogTerm = (double) Math.log(halfLength
				/ (MyosinV.myoVRadius));
		double PiViscosityLength = (double) Math.PI * Constants.viscosity
				* length;

		// double pValue = (2*radius)/length;
		// the Sigmas would normally be looked up from a table given the pValue.
		// Jacked from page 107 of Jonathan Howard's Mechanics of Motor Proteins
		// and the Cytoskeleton
		double rotSigma = -0.66;
		double parSigma = -0.20;
		double perpSigma = 0.86;

		double pValue = MyosinV.myoVRadius / halfLength;

		// System.out.println();
		// System.out.println("pValue="+pValue+"  radius="+radius+"  length="+length);
		boolean found = false;
		for (int i = 0; i < Drag.cylinderDragValues.length && !found; i++) {
			if (pValue >= Drag.cylinderDragValues[i][Drag.pValueIndex]) {
				parSigma = Drag.cylinderDragValues[i][Drag.parIndex];
				perpSigma = Drag.cylinderDragValues[i][Drag.perpIndex];
				rotSigma = Drag.cylinderDragValues[i][Drag.rotIndex];
				found = true;
				// System.out.println("using pValue="+pValue+"  index="+i+"  parSigma="+parSigma+"  perpSigma="+perpSigma+"  rotSigma="+rotSigma);
			}
		}

		if (!found) {
			int defaultIndex = Drag.cylinderDragValues.length - 1;
			parSigma = Drag.cylinderDragValues[defaultIndex][Drag.parIndex];
			perpSigma = Drag.cylinderDragValues[defaultIndex][Drag.perpIndex];
			rotSigma = Drag.cylinderDragValues[defaultIndex][Drag.rotIndex];
			// System.out.println("pValue outside of range so using default");
		}

		rotGam = (PiViscosityLength * length * length)
				/ (3 * (naturalLogTerm + rotSigma));
		myosinVtailTransGamma = (2 * PiViscosityLength)
				/ (naturalLogTerm + parSigma);
		myosinVtailPerpGamma = (4 * PiViscosityLength)
				/ (naturalLogTerm + perpSigma);
		/*
		 * System.out.println("rotationalGamma="+rotationalGamma);
		 * System.out.println
		 * ("translationalGammaParallel="+translationalGammaParallel);
		 * System.out.println("translationalGammaPerpendicular="+
		 * translationalGammaPerpendicular);
		 */

		diffPar = Constants.kT / myosinVtailTransGamma;
		diffPerp = Constants.kT / myosinVtailPerpGamma;
		diffRot = Constants.kT / rotGam;

		setDiffMultipliers(1, 1); // default to no multiplier for diffusion
									// terms

		if (diffPerp > diffPar) {
			outOfPlaneDiff = outOfPlaneDiffFactor * diffPerp;
		} else {
			outOfPlaneDiff = outOfPlaneDiffFactor * diffPar;
		}
	}

//	public void drawYourself (Graphics g, double scale, double [] offset) {
//		// put the code here to draw the object on "g"
//	}
	
//	public void setPaintColor(Color mmm) {
//		myColor = mmm;
//	}
		
	public void drawYourself(Graphics g, double scale, double[] offset) {
		setMyosinVTailParams();
		if (offCortex()) {
			System.out.println("off!");
			return;
		} // no rendering if off cortex
		Graphics2D g2d = (Graphics2D) g;
		g2d.setPaint(myosinVTailColor);

		// draw from end1
		A1.add(end1, -halfWidth, upVect);
		A2.add(end1, length, uVect, -halfWidth, upVect);
		B1.add(end1, halfWidth, upVect);
		B2.add(end1, length, uVect, halfWidth, upVect);
		drawPerimeter(g2d);

		if (wrapped) { // if spanning a boundary we need to draw from both ends
			// draw from end2
			A1.add(end2, -length, uVect, -halfWidth, upVect);
			A2.add(end2, -halfWidth, upVect);
			B1.add(end2, -length, uVect, halfWidth, upVect);
			B2.add(end2, halfWidth, upVect);
			drawPerimeter(g2d);

		if (showInfo) {
			float END1X = (float) ((end1.x - offset[0]) * scale);
			float END1Y = (float) ((end1.y - offset[0]) * scale);

			float END2X = (float) ((end2.x - offset[0]) * scale);
			float END2Y = (float) ((end2.y - offset[0]) * scale);

			float CMX = (float) ((cm.x - offset[0]) * scale);
			float CMY = (float) ((cm.y - offset[0]) * scale);

			// g2d.drawString("A1",(int)A1X,(int)A1Y);
			// g2d.drawString("A2",(int)A2X,(int)A2Y);
			// g2d.drawString("B1",(int)B1X,(int)B1Y);
			// g2d.drawString("B2",(int)B2X,(int)B2Y);

			g2d.drawString("end1", (int) END1X, (int) END1Y);
			g2d.drawString("end2", (int) END2X, (int) END2Y);

			int EndPixelRadius = (int) (10 * scale);
			int EndPixelDiameter = (int) (2 * 10 * scale);

			if (showEndsInfo) {
				int A1X = A1.getPixX();
				int A1Y = A1.getPixY();

				int A2X = A2.getPixX();
				int A2Y = A2.getPixY();

				int B1X = B1.getPixX();
				int B1Y = B1.getPixY();

				int B2X = B2.getPixX();
				int B2Y = B2.getPixY();
				g2d.setPaint(Color.CYAN);
				g.fillOval((int) (A1X - EndPixelRadius),
						(int) (A1Y - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);
				g2d.setPaint(Color.WHITE);
				g.fillOval((int) (A2X - EndPixelRadius),
						(int) (A2Y - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);
				g2d.setPaint(Color.RED);
				g.fillOval((int) (B1X - EndPixelRadius),
						(int) (B1Y - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);
				g2d.setPaint(Color.YELLOW);
				g.fillOval((int) (B2X - EndPixelRadius),
						(int) (B2Y - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);

				g2d.setPaint(Color.MAGENTA);
				g.fillOval((int) (CMX - EndPixelRadius),
						(int) (CMY - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);

				g2d.setPaint(Color.YELLOW);
				g.fillOval((int) (END1X - EndPixelRadius),
						(int) (END1Y - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);

				g2d.setPaint(Color.CYAN);
				g.fillOval((int) (END2X - EndPixelRadius),
						(int) (END2Y - EndPixelRadius), EndPixelDiameter,
						EndPixelDiameter);
			}

			if (showUVect) {
				int U1X = (int) ((cm.x - offset[0]) * scale);
				int U1Y = (int) ((cm.y - offset[0]) * scale);

				double lengthOfUVectDrawing = 20;
				double newx = cm.x + lengthOfUVectDrawing * uVect.x;
				double newy = cm.y + lengthOfUVectDrawing * uVect.y;
				int U2X = (int) ((newx - offset[0]) * scale);
				int U2Y = (int) ((newy - offset[0]) * scale);

				g2d.setPaint(Color.cyan);
				g2d.drawLine(U1X, U1Y, U2X, U2Y);

				int UP1X = cm.getPixX();
				int UP1Y = cm.getPixY();

				newx = cm.x + lengthOfUVectDrawing * upVect.x;
				newy = cm.y + lengthOfUVectDrawing * upVect.y;
				int UP2X = (int) ((newx - offset[0]) * scale);
				int UP2Y = (int) ((newy - offset[0]) * scale);

				g2d.setPaint(Color.magenta);
				g2d.drawLine(UP1X, UP1Y, UP2X, UP2Y);
				}
			}
		}
	}
		/*
		 * int U1X=(int)((centermass.x-offset[0])*scale); int
		 * U1Y=(int)((centermass.y-offset[0])*scale);
		 *
		 * g2d.drawString(this+" myID= "+Integer.toString(myMyosinVTailamentNumber)+"  "
		 * +Double.toString(brownianTorque),U1X,50+myMyosinVTailamentNumber*20);
		 */

	public void drawPerimeter(Graphics2D g2d) {
		setMyosinVTailParams();
		int A1X = A1.getPixX();
		int A1Y = A1.getPixY();

		int A2X = A2.getPixX();
		int A2Y = A2.getPixY();

		int B1X = B1.getPixX();
		int B1Y = B1.getPixY();

		int B2X = B2.getPixX();
		int B2Y = B2.getPixY();

		perimeter.reset();
		perimeter.moveTo(A1X, A1Y);
		perimeter.lineTo(A2X, A2Y);
		perimeter.lineTo(B2X, B2Y);
		perimeter.lineTo(B1X, B1Y);
		perimeter.closePath();

		g2d.draw(perimeter);
		Composite originalComposite = g2d.getComposite();
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,
				0.5f));
		g2d.fill(perimeter);
		g2d.setComposite(originalComposite);
	}

	public void setMyosinVTailParams() {
	}

	public void releaseAll() {
		for (int i = 0; i < myMyosinVs.length; i++) {
			myMyosinVs[i].releaseMon();
		}
	}

	public boolean allMyosinVsFree() {
		for (int i = 0; i < myMyosinVs.length; i++) {
			if (!myMyosinVs[i].isFree()) {
				return false;
			}
		}
		return true;
	}


	static public void makeMyosinVTail(double initX, double initY,
			double uVectX, double uVectY) {
		MyosinVTail mt = new MyosinVTail(initX, initY, uVectX,
				uVectY);
	}

	public static void makeTestBlob(int numMyosinV, double blobRadius) {
		for (int i = 0; i < numMyosinV; i++) {
			double theta = 2 * Math.PI * Math.random();
			double initX = Math.random() * blobRadius * Math.cos(theta)
					+ (Sim2D.xDimension) / 2;
			double initY = Math.random() * blobRadius * Math.sin(theta)
					+ (Sim2D.yDimension) / 2;
			MyosinVTail mt = new MyosinVTail(initX, initY,
					(2 * Math.random() - 1), (2 * Math.random() - 1));
		}
	}

	public static void addMyosinVTail(MyosinVTail newMyosinVTail) {
		newMyosinVTail.myMyosinVTailNumber = myosinVTailCt;
		theMyosinVTails[myosinVTailCt] = newMyosinVTail;
		myosinVTailCt++;
	}

	public void die() {
		myAmyosinVs = null;
		myBmyosinVs = null;
		myMyosinVs = null;
	}

	public void removeMe() {
		super.removeMe();
		int swapId = myMyosinVTailNumber;
		theMyosinVTails[swapId] = theMyosinVTails[myosinVTailCt - 1];
		theMyosinVTails[swapId].myMyosinVTailNumber = swapId;
		myosinVTailCt--;
	}
	
	public MyosinV getFreeMyosinV() {
		for (int i = 0; i < myMyosinVs.length; i++) {
			if (myMyosinVs[i].isFree()) {
				return myMyosinVs[i];
			}
		}
		return null;
	}
	
//	public void checkBoundToActin() {
//		headOneState = myAmyosinVs[0].state;
//		headTwoState = myAmyosinVs[1].state;
//		if (!(headOneState == 0 && headTwoState == 0))
//			boundToActin = true;
//		else
//			boundToActin = false;
//		return;
//	}
	
	public void getHeadPositions() {
		headOneState = myAmyosinVs[0].state;
		headTwoState = myAmyosinVs[1].state;
		// first checking the unbound and single-head bound states...
		if (headOneState == 0 && headTwoState == 0) {
			headOnePosition = UNBOUND;
			myAmyosinVs[0].setHeadPosition(UNBOUND);
			headTwoPosition = UNBOUND;
			myAmyosinVs[1].setHeadPosition(UNBOUND);
		}
		if (headOneState == 0 && headTwoState != 0) {
			headOnePosition = UNBOUND;
			myAmyosinVs[0].setHeadPosition(UNBOUND);
			headTwoPosition = SINGLE;
			myAmyosinVs[1].setHeadPosition(SINGLE);
		}
		if (headOneState != 0 && headTwoState == 0) {
			headOnePosition = SINGLE;
			myAmyosinVs[0].setHeadPosition(SINGLE);
			headTwoPosition = UNBOUND;
			myAmyosinVs[1].setHeadPosition(UNBOUND);
		}
		// next check for spanning across different filaments
		if (headOneState != 0 && headTwoState != 0 && myAmyosinVs[0].myBoundActinFil != myAmyosinVs[1].myBoundActinFil) {
				headOnePosition = INTERSECTION;
				myAmyosinVs[0].setHeadPosition(INTERSECTION);
				headTwoPosition = INTERSECTION;
				myAmyosinVs[1].setHeadPosition(INTERSECTION);
		}
		if (headOneState != 0 && headTwoState != 0 && myAmyosinVs[0].myBoundActinFil == myAmyosinVs[1].myBoundActinFil) {
		// finally decide which is leading head, which is trailing head
			if (myAmyosinVs[0].myBoundMonomer > myAmyosinVs[1].myBoundMonomer) {
				headOnePosition = LEADING;
				myAmyosinVs[0].setHeadPosition(LEADING);
				headTwoPosition = TRAILING;
				myAmyosinVs[1].setHeadPosition(TRAILING);
			}
			if (myAmyosinVs[0].myBoundMonomer < myAmyosinVs[1].myBoundMonomer) {
				headOnePosition = TRAILING;
				myAmyosinVs[0].setHeadPosition(TRAILING);
				headTwoPosition = LEADING;
				myAmyosinVs[1].setHeadPosition(LEADING);
			}
		}
	}
	
	
	public void setRunStatus(int st) {
		runStatus = st;
	}

	public void initiateRun() {
		runStartTime = Sim2D.simulationTime;
		runStartX = cm.x;
		runStartY = cm.y;
		runCount ++;
		thisRun = runCount;
		setRunStatus(RUNNING);
// diagnostic: start run position and time, output to console while code runs
//		System.out.println("run " + thisRun + " initiated at (" + runStartX + "," + runStartY + ")");
	}

	public void endRun() {
		runEndTime = Sim2D.simulationTime;
		runEndX = cm.x;
		runEndY = cm.y;
		runtime = runEndTime - runStartTime;
		runlength = Math.pow (((Math.pow ((runEndX - runStartX), 2)) + (Math.pow ((runEndY - runStartY), 2))), 0.5);
		runvelocity = (runlength)/(runtime);
		setRunStatus(NOT_RUNNING);
		if (runtime >= minRunTime && runlength >= jCutoff) {
			storeRunValues();
		}
// diagnostic: end run position and time, output to console while code runs
//		System.out.println("run " + thisRun + " ended at (" + runEndX + "," + runEndY + ")");
//		System.out.println("run " + thisRun + ": length = " + runlength + ", time = " + runtime + ", velocity = " + runvelocity);
//		System.out.println("runlength = " + runlength);// + "; velocity = " + runvelocity);

//		HERE! THIS IS WHERE I OUTPUT TO CONSOLE FOR COPY-PASTE TO MATLAB!
//		System.out.println(runlength + ", " + runvelocity + ";");// + "; velocity = " + runvelocity);
	}
	
	public void storeRunValues () {
		runlengthHistogram.addValue(runlength);
		runVelocityHistogram.addValue(runvelocity);
		runTimeHistogram.addValue(runtime);
		runlengthsList.add(runlength);
		runVelocitiesList.add(runvelocity);
		runlengthCt++;
//		System.out.println("values added! "+runlength);
//		calculateOutputValues();
	}
	
	static public void writeHistograms() {
		runTimeHistogram.writeToFile();
		runlengthHistogram.writeToFile();
		runVelocityHistogram.writeToFile();
	}
	
	//this is a ghetto output method, printing into the console frame...
	static public void calculateOutputValues() {
		meanRunlength = runlengthHistogram.getMean() - jCutoff;
		meanRunVelocity = runVelocityHistogram.getMean();
		sdRunlength = runlengthHistogram.getStdev();
		sdRunVelocity = runVelocityHistogram.getStdev();
		nRunlengths = runlengthHistogram.getPointCt();
		nRunVelocities = runVelocityHistogram.getPointCt();
		seRunlength = calculateSEcountFromHist(runlengthHistogram,meanRunlength);
		meanRunlengthList = calculateMean(runlengthsList) - jCutoff;
		sdRunlengthList = calculateSD(runlengthsList,meanRunlengthList);
		meanRunVelocitiesList = calculateMean(runVelocitiesList);
		sdRunVelocitiesList = calculateSD(runVelocitiesList,meanRunVelocitiesList);
		seRunlengthList = calculateSEcount(runlengthsList,meanRunlengthList);
		
//		for(int i=0;i<runlengthHistogram.binCt;i++) {
//			nRunlengths += runlengthHistogram.binSums[i];
//		}
//		for(int i=0;i<runVelocityHistogram.binCt;i++) {
//			nRunVelocities += runVelocityHistogram.binSums[i];
//		}
		System.out.println("X_hist = " + meanRunlength + " ± " + seRunlength + " (mean ± se); n = " + nRunlengths);
		System.out.println("V_hist = " + meanRunVelocity + " ± " + sdRunVelocity + " (mean ± sd); n = " + nRunVelocities);
		System.out.println("X_list = " + meanRunlengthList +" ± "+ seRunlengthList + " (mean ± se); n = " +runlengthCt);
		System.out.println("V_list = " + meanRunVelocitiesList +" ± "+ sdRunVelocitiesList + " (mean ± sd); n = " +runlengthCt);
	}
	
	static public double calculateMean(Vector v_in) {
		int n = v_in.size();
		double d;
		Double[] vArray = new Double[n];
		v_in.toArray(vArray);
		double m_out;
				
		double sum = 0;
		for(int j = 0; j < n; j++) {
			sum = sum + vArray[j];
		}
		m_out = sum/(n+1);
		return m_out;
	}

	static public double calculateSD(Vector v_in, double m_in) {
		double sd_out;
		int n = v_in.size();
		double d;
		Double[] vArray = new Double[n];
		v_in.toArray(vArray);
		
		double sumdiffsq = 0;
		for(int k = 0; k < n; k++) {
			sumdiffsq = sumdiffsq + Math.pow((vArray[k]-m_in),2);
		}
		sd_out = Math.pow((sumdiffsq/(n+1)),0.5);
		return sd_out;
	}
	
	static public double calculateSEcount(Vector v_in, double m_in) {
		double se_out;
		int n = v_in.size();
		se_out = m_in / (Math.pow((n+1),0.5));
		return se_out;
	}
	
	static public double calculateSEcountFromHist(HistogramPlus hist, double m_in) {
		double se_out;
		int n = hist.getPointCt();
		se_out = m_in / (Math.pow(n,0.5));
		return se_out;
	}

	public double getBindingProb(Actin a) {
		if (a instanceof Bundle) {
			if (headOnePosition == UNBOUND && headTwoPosition == UNBOUND) {
				return MyosinV.myoVBundleInitBindingProb;
			} else {
				return MyosinV.myoVBundleSecondHeadBindingProb;
			}
		}
		if (!(a instanceof Bundle)) {
			if (headOnePosition == UNBOUND && headTwoPosition == UNBOUND) {
				return MyosinV.myoVInitBindingProb;
			} else {
				return MyosinV.myoVSecondHeadBindingProb;
			}
		} else {
			System.out.println("can't determine myoVBindingProb");
			return 0;
//			System.out.println("can't determine myoVBindingProb");
		}
	}
	
	public void assignLastMVTailPosition(Point2D sourcePt) {
		lastXpos = sourcePt.x;
		lastYpos = sourcePt.y;
//		System.out.println("MyosinVTail.assignMVTailPosition: xin = "+sourcePt.x+", xout = "+lastXpos);
	}
	
	public void assignNextMVTailPosition(Point2D sourcePt) {
		nextXpos = sourcePt.x;
		nextYpos = sourcePt.y;
//		System.out.println("MyosinVTail.assignMVTailPosition: xin = "+sourcePt.x+", xout = "+nextXpos);
	}

	public void reflectPositionOffBarrier(double xCorr, double yCorr) {
		cm.x = cm.x + xCorr;
		cm.y = cm.y + yCorr;
		myAmyosinVs[0].reflectOffBarrier(xCorr,yCorr);
		myAmyosinVs[0].reflectOffBarrier(xCorr,yCorr);
//		System.out.println("MyosinVTail.reflectPositionOffBarrier: run!");
	}

	
//	public class PositionInfo {
//		Point2D lastMVTailPosition = new Point2D ();
//		Point2D nextMVTailPosition = new Point2D ();
//	}
		static public void staticInit() {

		Parameters.addParameterSetListener (
			new ParameterSetListener() {
				public void parametersChanged() {
					length = MyosinV.myosinVTailLength;
					halfLength = 0.5 * length;
					width = MyosinV.myosinVTailWidth;
					halfWidth = 0.5 * width;
					collisionRadSqrd = MyosinV.myoVCollisionRadius*MyosinV.myoVCollisionRadius;
					nMyoVHeadsA = MyosinV.nMyosinVHeadsA;
					nMyoVHeadsB = MyosinV.nMyosinVHeadsB;
				}
			}
		);

		Sim2D.addProxy(new AgentProxy() {
			public void registerForCollisionDetection() {
				fillMVTailMesh();
			}
			public void checkCollisions(double simTime) {
				myCollisionDetector.checkCollisions(simTime);
			}
			public void init() {
				myCollisionDetector.init();
			}
			public void step(double dT) {
				//  MyosinVTails move and update the positions of their heads
				stepAllMyosinVtails(dT);
			}
			public void checkBarrierCrossings() {
				myCollisionDetector.checkBarrierCrossings();
			}
			public void initDiagnostics() {
				runTimeHistogram = new HistogramPlus(20,0,10,Sim2D.basePath,"runTimes",null,true,false,false,false);
				runlengthHistogram = new HistogramPlus(20,0,1000,Sim2D.basePath,"runlengths",null,true,false,false,false);
				runVelocityHistogram = new HistogramPlus(20,0,1000,Sim2D.basePath,"runVelocities",null,true,false,false,false);
		//		if (Sim2D.writeState) {
		//			writeHistograms();
		//		}
			}
		
			public void reset() {
				for (int i = 0; i < myosinVTailCt; i++) {
					theMyosinVTails[i].die();
					theMyosinVTails[i] = null;
				}
				myosinVTailCt = 0;
			}

		});
	}


}


