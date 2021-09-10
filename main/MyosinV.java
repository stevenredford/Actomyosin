package main;


import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics;
import java.awt.Graphics2D;

import collision.*;
import util.*;
import parameters.*;
import gui.*;

public class MyosinV extends Thing {
	
	static CollisionDetector myCollisionDetector = new MyoVCollisionDetector();
	
	/** histograms containing runlengths and velocities */
	static public HistogramPlus runlengths;
	static public HistogramPlus velocities;
	static public HistogramPlus lifetimeHistogram; // time bound to actin for each myosin

	/** timestamp when myosin V head binds to actin, amount of time bound */
	double bindStartTime;
	
	double dwellTime;
	
	/** An array containing all the MyosinVs */
	static public MyosinV [] theMyosinVs = new MyosinV [50000];
	
	/** Current number of MyosinVs. */
	static public int myosinVCt = 0;
		
	/** Diameter of a myosinV head in nm */
	static double headDiameter;

	/** How close an actin filament must be to the stalk's base for binding to occur. */
	static public double minBindingDistance;
	
	static public int
	
		FREE = 0,
		BOUND = 1,
		BOUND_PS = 2;

	public int state = FREE;

	static public int
		UNBOUND = 0,
		SINGLE = 1,
		LEADING = 2,
		TRAILING = 3,
		INTERSECTION = 4;
	
	public int headPos = UNBOUND;
		
	/** Whether this myosinV head is at the "A" end or the "B" end of the filament. */
	double orientation;
		
	/** Position of this myosinV in the array of all myosinVs. */
	int myMyosinVNumber;

	/** The myosinVtail this myosinV is attached to. */
	MyosinVTail myMyosinVTail;

	/** True if this myosinV head is at "A" end of myosinVTail. */
	boolean iAmAnAMyoV;
	
	/** The actin monomer this myosinV head is bound to. */
	Monomer boundMon;
	
	public int myBoundMonomer;
	public int myBoundActinFil;
		
	double bindingProb;
	
	/** The Point2D at which the myosinV stalk attaches to the myosinVtail expressed
	 * with respect to the myosinVTails local coordinate system. */
	Point2D mfAttachmentPtLocal=new Point2D();

	/** The Point2D at which the myosinV stalk attaches in the myosinVtail in global coordinates. */
	public Point2D mfAttachmentPt=new Point2D();
		
	/** The Point2D at which the myosinV head binds to the myosinVtail. */
	public Point2D bindingSite = new Point2D();
	
	Point2D forceSite = new Point2D();
		
	/** A unit vector Point2Ding from the myosinVTail center to the attachment site for this myosinV. */
	Point2D oriUVect = new Point2D();
	
	/** Force this myosinV head exerts on actin filament. Reused in force calculations. */
	Point2D forceVector = new Point2D();
	
	/** Magnitude of the force this myosinV head exerts on actin filament. Reused in force calculations. */
	double forceMag = 0;
	
	/** Time-averaged force on this myosinV head - used for force-dependent binding kinetics. */
	double forceAv;
	
	/** How far the myosinV-actin bind is stretched defined as the distance between the
	 * myosinV/myosinVTail attachment site and the myosinV/Factin binding site minus the stalk length.
	 */
	double bondStretch = 0;

	/** Store some values for when and where a run starts */
	
	double runStartT = -1;
	double runStartX = -1;
	double runStartY = -1;
	
	/** Color to draw this myosinV. */
	Color myosinVColor = Color.blue;
	
	/** Used during rendering. */
	Point2D theLine = new Point2D();			// reused in rendering
	
	/** flag used in rendering from files */
	boolean fileMyoVOffCortex = false;
		
	public MyosinV (MyosinVTail myf) {
		super(0,0);
		myMyosinVTail = myf;
		addMyosinV(this);
	}
	
	static public void stepAllMyosinVs(double dT) {
		for (int i = 0; i < myosinVCt; i++) {
			theMyosinVs[i].step(dT);
		}
	}
	
	public void setInitialPositionRelativeToMyosinVTail(Point2D relMyosinVTailPosition){
		this.mfAttachmentPtLocal.set(relMyosinVTailPosition);
	}
	
	public void setState(int st) {
		state = st;
	}

	public void setHeadPosition(int hp) {
		headPos = hp;
	}

	public void step (double dT) {
		
		doStateChangeKinetics(dT);
		updateSites();
//		System.out.println("myosinV head number " + this.myMyosinVNumber + " state is " + this.state + " at time step " + Sim2D.counter);
		doForces();
	}
	
	private void doStateChangeKinetics(double dT) {
		if(state == FREE) {
//			System.out.println("myosinV head number " + this.myMyosinVNumber + " is free");
			bindingSite.zero();
			return;
		}
		else if(state == BOUND) {
//			System.out.println("myosinV head number " + this.myMyosinVNumber + " is bound to monomer number " + boundMon.myMonNumber + " at time step " + Sim2D.counter);
			// test for release transition
			if(Math.random() < getReleaseProb()*dT) {
//			if(Math.random() < Constants.myoVUniformReleaseProb*Sim2D.deltaT) {
				releaseMon();
//				System.out.println("myosinV head number " + this.myMyosinVNumber + " is now free");
				return;
			}
			// test for powerstroke transition
			if(Math.random() < getPowerStrokeProb()*dT) {
				setState(BOUND_PS);
//				System.out.println("myosinV head number " + this.myMyosinVNumber + " has performed a power stroke");
			}
			
			return;
		}
			
		else if(state == BOUND_PS) {
			// test for release transition
			if(Math.random() < getReleaseProb()*dT) {
//			if(Math.random() < Constants.myoVUniformReleaseProb*Sim2D.deltaT) {
				releaseMon();
				return;
			}
			takeProcessiveStep(dT);
			return;
		}
	}
		
// don't need this method anymore...
	private void takeProcessiveStep(double dT) {
		if(Math.random() < myoVUniformRebindingProb*dT) {
			Actin tmpFilament = boundMon.myFilament ;
			double nxtMonPosition = boundMon.myFilIndex*Actin.monLength + myoVStepSize;
			Monomer tmpMon = tmpFilament.getMonomerAt(nxtMonPosition);
			if(tmpMon != null && tmpMon.isFreeMyosinV()) { // tmpMon == null -> step beyond barbed end
//				System.out.println("Monomer at next position myfilindex= " + tmpMon.myFilIndex + " and boundmon.myfilindex is= " + boundMon.myFilIndex);
				//System.out.println("for myoVhead num " + myMyosinVNumber + ", original monomer filament index is " + boundMon.myFilIndex);
				releaseMon();
				bindMon(tmpMon);
			}
		}
	}

	/* Calculate monomer binding sites and force sites with respect to new monomer location. */
	private void updateSites() {
	
		if(isFree()) return;
		
		bindingSite.set(boundMon.getLocation());
		if(state == BOUND) {
			forceSite.set(bindingSite);
		}
		else { // state == BOUND_PS
			forceSite.add(bindingSite, myoVStepSize,boundMon.myFilament.uVect);
		}
	}
	
	/* calculate and apply motor forces. */
	private void doForces() {
		if(state == FREE) return;
		forceVector.getVector(mfAttachmentPt,forceSite);
		forceVector.uVect();
		bondStretch = Point2D.getDistance(forceSite, mfAttachmentPt) - minBindingDistance;  // calculate distance between center of particles
		if (bondStretch < 0) { bondStretch = 0; }
		if (!Point2D.pointOK(forceVector)) {
			System.out.println ("forceVector in MyosinV.applyForce() is NaN");
			System.out.println ("bondStretch = " + bondStretch);
			System.out.println ("mFAttachmentSite is " + mfAttachmentPt.reportVals());
			System.out.println ("forceSite is " + forceSite.reportVals());
			System.out.println ("boundMon is " + boundMon);
			System.out.println ("boundMon.myFilament is " + boundMon.myFilament + " with cm at " + boundMon.myFilament.cm.reportVals());
		}
		forceMag = myoVSpringConstant*bondStretch;
		forceAv = (1-myoVForceAvFactor)*forceAv + myoVForceAvFactor*forceMag;
		forceVector.scale(forceMag);
		myMyosinVTail.addForce(forceVector, mfAttachmentPt);
		forceVector.scale(-1);
		boundMon.addForce(forceVector);
	}
	
	public double getAngle(){
		oriUVect.x = orientation*myMyosinVTail.uVect.x;
		oriUVect.y = orientation*myMyosinVTail.uVect.y;
		
		return Point2D.getAngle(boundMon.myFilament.uVect, oriUVect);
	}
	
	public Point2D getAim(){
		oriUVect.x = orientation*myMyosinVTail.uVect.x;
		oriUVect.y = orientation*myMyosinVTail.uVect.y;
		return oriUVect;
	}

	public double getReleaseProb () {
	if (boundMon.myFilament instanceof Bundle) {
		return getReleaseProbBundle ();
		} else {
			return getReleaseProbSingleFil ();
		}
	}
	
	public double getReleaseProbSingleFil () {
		if (myoVUseForceBasedRelease && headPos == LEADING) {
			return myoVLeadForceBasedReleaseBase*Math.exp(forceMag*myoVLeadForceBasedReleaseExp);
		}
		if (myoVUseForceBasedRelease && headPos == TRAILING) {
			return myoVTrailForceBasedReleaseBase*Math.exp(forceMag*myoVTrailForceBasedReleaseExp);
		}
//		if (Constants.myoVUseForceBasedRelease && headPos != LEADING && headPos != TRAILING) {
//			return Constants.myoVForceBasedReleaseBase*Math.exp(forceMag*Constants.myoVForceBasedReleaseExp)*Sim2D.deltaT;
		if (!myoVUseForceBasedRelease && headPos == LEADING) {
			return myoVLeadReleaseProb;
		}
		if (!myoVUseForceBasedRelease && headPos == TRAILING) {
			return myoVTrailReleaseProb;
		}
		if (!myoVUseForceBasedRelease && headPos == INTERSECTION) {
			return myoVIntersectReleaseProb;
		}
		if (!myoVUseForceBasedRelease && headPos == SINGLE) {
			return myoVSingleReleaseProb;
		} else {
			return myoVUniformReleaseProb;
		}
	}
	
	public double getReleaseProbBundle () {
		if (myoVUseForceBasedRelease && headPos == LEADING) {
			return myoVBundleLeadForceBasedReleaseBase*Math.exp(forceMag*myoVBundleLeadForceBasedReleaseExp);
		}
		if (myoVUseForceBasedRelease && headPos == TRAILING) {
			return myoVBundleTrailForceBasedReleaseBase*Math.exp(forceMag*myoVBundleTrailForceBasedReleaseExp);
		}
//		if (Constants.myoVUseForceBasedRelease && headPos != LEADING && headPos != TRAILING) {
//			return Constants.myoVForceBasedReleaseBase*Math.exp(forceMag*Constants.myoVForceBasedReleaseExp)*Sim2D.deltaT;
		if (!myoVUseForceBasedRelease && headPos == LEADING) {
			return myoVBundleLeadReleaseProb;
		}
		if (!myoVUseForceBasedRelease && headPos == TRAILING) {
			return myoVBundleTrailReleaseProb;
		}
		if (!myoVUseForceBasedRelease && headPos == INTERSECTION) {
			return myoVIntersectReleaseProb;
		}
		if (!myoVUseForceBasedRelease && headPos == SINGLE) {
			return myoVBundleSingleReleaseProb;
		} else {
			return myoVBundleUniformReleaseProb;
		}
	}

	public double getPowerStrokeProb () {
		if (myoVUseForceBasedRelease && headPos == LEADING) {
			return myoVLeadForceBasedPowerStrokeProb;
		}
		if (myoVUseForceBasedRelease && headPos == TRAILING) {
			return myoVTrailForceBasedPowerStrokeProb;
		}
		if (!myoVUseForceBasedRelease && headPos == LEADING) {
			return myoVLeadPowerStrokeProb;
		}
		if (!myoVUseForceBasedRelease && headPos == TRAILING) {
			return myoVTrailPowerStrokeProb;
		}
		if (!myoVUseForceBasedRelease && headPos == INTERSECTION) {
			return myoVIntersectPowerStrokeProb;
		}
		if (!myoVUseForceBasedRelease && headPos == SINGLE) {
			return myoVSinglePowerStrokeProb;
		} else {
			return myoVPowerStrokeProb;
		}
	}
		
	public void drawYourself (Graphics g, double scale, double [] offset) {
		if (myMyosinVTail.offCortex()) { return; }  // don't render if my myosinV tail is off-cortex
		Graphics2D g2d=(Graphics2D)g;
		Composite originalComposite = g2d.getComposite();
   	 	
		int xOrigin = mfAttachmentPt.getPixX();
		int yOrigin = mfAttachmentPt.getPixY();
		int pixelDiameter = (int) (headDiameter*scale);
		int pixelRadius = pixelDiameter/2;
		
		g2d.setPaint(myosinVColor);
		g2d.drawOval(xOrigin-pixelRadius, yOrigin-pixelRadius, pixelDiameter, pixelDiameter);
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
		g2d.fillOval(xOrigin-pixelRadius, yOrigin-pixelRadius, pixelDiameter, pixelDiameter);
		g2d.setComposite(originalComposite);

		if (!isFree()){
			g2d.setPaint(Color.cyan);
			int xBSite = bindingSite.getPixX();
			int yBSite = bindingSite.getPixY();
			g2d.drawOval(xBSite-pixelRadius, yBSite-pixelRadius, pixelDiameter, pixelDiameter);
			g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
			g2d.fillOval(xBSite-pixelRadius, yBSite-pixelRadius, pixelDiameter, pixelDiameter);
			g2d.setComposite(originalComposite);
			
			// line rendering that works with crossing boundaries
			theLine.getVector(bindingSite,mfAttachmentPt);
			int theLineX = theLine.getPixX();
			int theLineY = theLine.getPixY();
			
			g2d.setPaint(Color.cyan);
			g2d.drawLine(xOrigin, yOrigin, xOrigin-theLineX, yOrigin-theLineY);
			g2d.drawLine(xBSite,yBSite,xBSite+theLineX,yBSite+theLineY);
		}
	}
		
	public boolean isFree (){
		//monLinkSanityCk();
		return boundMon == null;
	}
	
	public boolean canBind () {
		if (iAmAnAMyoV) {
			if (myMyosinVTail.zPosA != 0) { return false; }
		} else {
			if (myMyosinVTail.zPosB != 0) { return false; }
		}
		return isFree();
	}
	

 	public void bindMon (Monomer mon) {
 		releaseMon();
 		mon.bindMyoV(this);
 		setState(BOUND);
 		myBoundMonomer = boundMon.myMonNumber;
 		myBoundActinFil = boundMon.myFilament.myActinNumber;
 	}
 	
 	public void releaseMon () {
		if(boundMon!=null){
			boundMon.releaseMyoV();
		}
		setState(FREE);
		boundMon = null;
	}
 	
 	public void monLinkSanityCk () {
 		if (boundMon != null) {
 			if (boundMon.boundMyoV != this) { releaseMon(); }
 		}
 	}
 	
	public static void reportBoundMyoVs () {
		int boundMyoVs = 0;
		for (int i=0;i<myosinVCt;i++) {
			if (!theMyosinVs[i].isFree()) { boundMyoVs++; }
		}
		System.out.println ("MyosinVs bound to monomers = " + boundMyoVs);
	}
 	
 	public static void addMyosinV (MyosinV newMyosinV) {
		newMyosinV.myMyosinVNumber = myosinVCt;
		theMyosinVs[myosinVCt] = newMyosinV;
		myosinVCt ++;
	}
	
	public void die() {
		boundMon = null;
		myMyosinVTail = null;
	}
	
	public void removeMe () {
		super.removeMe();
		int swapId = myMyosinVNumber;
		theMyosinVs[swapId] = theMyosinVs[myosinVCt-1];
		theMyosinVs[swapId].myMyosinVNumber = swapId;
		myosinVCt --;
	}
	
	public double getMyosinVBindingProb(Actin a) {
		MyosinVTail tail = myMyosinVTail;
		return tail.getBindingProb(a);
	}

	public void reflectOffBarrier(double xCorr, double yCorr) {
		this.cm.x = this.cm.x + xCorr;
		this.cm.y = this.cm.y + yCorr;
	}

	static boolean myoVUseForceBasedRelease = false;
	static boolean initMyoVTailstate = true;
	static int nMyosinVHeadsA = 2;
	static int nMyosinVHeadsB = 0;
	static int nMyosinVtails = 1;
	static double myosinVTailLength = 80;  // MV ~80, MX ~21
	static double myosinVTailWidth = 3;
	static double myoVSpringConstant = 1;
	static double myoVRadius = 8;
	static double myoVStalkLength = 24;  // MV ~24, MX ~13
	static double myoVStepSize = 72.0;  // MV ~72.0, MX ~36.0
	static double myoVForceBasedReleaseBase = 50;
	static double myoVLeadForceBasedReleaseBase = 10;
	static double myoVTrailForceBasedReleaseBase = 10;
	static public double myoVBundleLeadForceBasedReleaseBase = 10;
	static public double myoVBundleTrailForceBasedReleaseBase = 10;
	static double myoVForceBasedReleaseExp = 0.2;
	static double myoVLeadForceBasedReleaseExp = 0.1;
	static double myoVTrailForceBasedReleaseExp = 0.1;
	static public double myoVBundleLeadForceBasedReleaseExp = 0.1;
	static public double myoVBundleTrailForceBasedReleaseExp = 0.1;
	static double myoVInitBindingProb = 900; // MV ~400
	static public double myoVBundleInitBindingProb = 900; // MV ~200
	static double myoVSecondHeadBindingProb = 900;  // MV ~200, MX ~400
	static public double myoVBundleSecondHeadBindingProb = 900;  // MV ~200, MX ~400
	static double myoVUniformReleaseProb = 100;
	static double myoVBundleUniformReleaseProb = 100;
	static double myoVLeadReleaseProb = 0;
	static double myoVTrailReleaseProb = 15;  // MV ~15, MX ~30
	static public double myoVBundleLeadReleaseProb = 0;
	static public double myoVBundleTrailReleaseProb = 15; // MV ~15, MX ~30
	static double myoVIntersectReleaseProb = 50;
	static double myoVSingleReleaseProb = 300;  // MV ~18, MX ~108 ?
	static public double myoVBundleSingleReleaseProb = 300;   //MV ~18, MX ~18 ?
	static double myoVUniformRebindingProb = 0;
	static double myoVPowerStrokeProb = 50;
	static double myoVLeadPowerStrokeProb = 50;
	static double myoVTrailPowerStrokeProb = 5;
	static double myoVIntersectPowerStrokeProb = 0;
	static double myoVSinglePowerStrokeProb = 0;
	static double myoVLeadForceBasedPowerStrokeProb = 800;
	static double myoVTrailForceBasedPowerStrokeProb = 0;
	static public double myoVAngleBias = -1.1;
	static double myoVForceAvFactor = 1;
	static public double myoVCollisionRadius = 10;
	static String className = new String("main.MyosinV");

	static void staticInit() {

		Parameters.addParameterSetListener (
			new ParameterSetListener() {
				public void parametersChanged() {
					headDiameter = 2.0*myoVRadius;
					minBindingDistance = myoVStalkLength + headDiameter;
				}
			}
		);
		
		
		
		Parameters.addParameter(className, "myoVUseForceBasedRelease",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVUseForceBasedRelease = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className, "initMyoVTailstate",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					initMyoVTailstate = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className, "nMyosinVHeadsA",2,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					nMyosinVHeadsA = (int)p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "nMyosinVHeadsB",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					nMyosinVHeadsB = (int)p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "nMyosinVtails",1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					nMyosinVtails = (int)p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myosinVTailLength",80,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myosinVTailLength = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myosinVTailWidth",3,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myosinVTailWidth = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVSpringConstant",1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVSpringConstant = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVRadius",8,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVRadius = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVStalkLength",24,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVStalkLength = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVStepSize",72,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVStepSize = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVForceBasedReleaseBase",50,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVForceBasedReleaseBase = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVLeadForceBasedReleaseBase",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVLeadForceBasedReleaseBase = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVTrailForceBasedReleaseBase",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVTrailForceBasedReleaseBase = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleLeadForceBasedReleaseBase",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleLeadForceBasedReleaseBase = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleTrailForceBasedReleaseBase",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleTrailForceBasedReleaseBase = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVForceBasedReleaseExp",0.2,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVForceBasedReleaseExp = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVLeadForceBasedReleaseExp",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVLeadForceBasedReleaseExp = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVTrailForceBasedReleaseExp",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVTrailForceBasedReleaseExp = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleLeadForceBasedReleaseExp",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleLeadForceBasedReleaseExp = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleTrailForceBasedReleaseExp",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleTrailForceBasedReleaseExp = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVInitBindingProb",900,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVInitBindingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleInitBindingProb",900,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVInitBindingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVSecondHeadBindingProb",900,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVSecondHeadBindingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleSecondHeadBindingProb",900,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleSecondHeadBindingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVUniformReleaseProb",100,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVUniformReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleUniformReleaseProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleUniformReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVLeadReleaseProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVLeadReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVTrailReleaseProb",15,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVTrailReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleLeadReleaseProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleLeadReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleTrailReleaseProb",15,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleTrailReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVIntersectReleaseProb",50,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVIntersectReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVSingleReleaseProb",300,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVSingleReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVBundleSingleReleaseProb",300,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVBundleSingleReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVUniformRebindingProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVUniformRebindingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVPowerStrokeProb",50,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVPowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVLeadPowerStrokeProb",50,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVLeadPowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVTrailPowerStrokeProb",5,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVTrailPowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVIntersectPowerStrokeProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVIntersectPowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVSinglePowerStrokeProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVSinglePowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVLeadForceBasedPowerStrokeProb",800,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVLeadForceBasedPowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVTrailForceBasedPowerStrokeProb",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVTrailForceBasedPowerStrokeProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVAngleBias",-1.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVAngleBias = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVForceAvFactor",1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVForceAvFactor = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoVCollisionRadius",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoVCollisionRadius = p.getDoubleValue();
				}
			}
		);
		
		
		Sim2D.addProxy(new AgentProxy() {
			public void checkCollisions(double simTime) {
				myCollisionDetector.checkCollisions(simTime);
			}
			public void init() {
				myCollisionDetector.init();
			}
			public void doForces(double dT) {
				// Myosin heads exert forces.
				stepAllMyosinVs(dT);
			}
			public void initDiagnostics() {
				runlengths = new HistogramPlus(20,0,1000,Sim2D.basePath,"runlengths",null,false,true,false,false);
				velocities = new HistogramPlus(20,0,500,Sim2D.basePath,"velocities",null,false,true,false,false);
			}
			
			public void reset() {
				for(int i = 0; i < myosinVCt; i++) {
					theMyosinVs[i].die();
					theMyosinVs[i] = null;
				}
				myosinVCt = 0;
				runlengths.clearBins();
				velocities.clearBins();
			}
		});
	}
}

/*
 *
 *
 *
 
 if (Actin.useBundling&&myBoundPartner.myFilament.bundled && Point2D.getDistance(myBoundPartner.centermass, myBoundPartner.myFilament.bEnd) < myBoundPartner.myFilament.length/2) {
				if(angle >= Math.PI-MyosinV.minBindAngle){
					if (myBoundPartner.prev != null){
						if (myBoundPartner.prev.isFree()){
							Monomer m = myBoundPartner.prev;
							myBoundPartner.release();
							m.addBoundPartner(this);
							addBoundPartner(m);
						}
					}
				}
				else myBoundPartner.release();
			}
 
 
 */



