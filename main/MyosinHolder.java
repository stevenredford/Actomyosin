package main;

/*
	Defines a point of attachment for myosin heads.  Extended to make e.g. Minifilaments
*/

import java.awt.Graphics;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Random;

import util.*;

public class MyosinHolder extends Thing {
	
	/** array holding all MyosinHolders */
	static MyosinHolder [] theMyoHolders = new MyosinHolder [10000];
	
	/** Current number of MyosinHolders. */
	static int myoHolderCt = 0;								// initialize the LongMiniF count
	
	/** I am "The One" for connectivity rendering/calculating */
	static MyosinHolder rootMyoHolder;
		
	/** Unique index of this MyosinHolder in the array of all MyosinHolders. */
	int myMyoHolderNumber;
	
	/** indicates if this MyosinHolder spans one or more boundaries */
	boolean wrapped = false;
	
	/** Maximum Z position. */
	static int maxZPos = 8;
	
	/** Reused for torque calculations. */
	Point2D torqArm = new Point2D();
	
	/** Reused for collision caclulations. */
	public Point2D tempP = new Point2D();
	
	/** Unit vector along minifilament axis. */
	public Point2D uVect = new Point2D ();

	/** Unit vector perpendicular to minifilament axis. */
	public Point2D upVect = new Point2D ();
	
	/** Array of all myosin heads. */
	Myosin [] myMyosins;
		
	/** drag coefficients, needed for calculating diffusion */
	public double myosinfilTransGamma = 0;
	public double myosinfilPerpGamma = 0;
	public double myosinRotGam = 0;
	
	/** For calculating diffusion. */
	double diffPar,diffPerp,diffRot;
	double curDPar,curDPerp,curDRot;  // normally equal to diffPar, etc, but used in step() to increase off cortex diffusivity
	double diffStepPar,diffStepPerp,diffStepRot;
	
	Random generator = new Random(Sim2D.generator.nextLong());
	Random jumpGen = new Random(Sim2D.generator.nextLong());

	/** For calculating forces due to diffusion. */
	Point2D brownianForceParallel=new Point2D();
	Point2D brownianForcePerpendicular=new Point2D();
	Point2D brownianForceRotational=new Point2D();
	
	/** For summing forces. */
	Point2D totalTranslationalForce=new Point2D();
	Point2D velocity=new Point2D();
	
	double torque = 0;
	double brownianTorque = 0;
	double deltaTheta = 0;
	double thetaNOW = 0;
	double thetaNEXT = 0;
	double angularVel = 0;
			
	
	/** for rendering */
	boolean showInfo=false;
	boolean showEndsInfo=false;
	boolean showUVect=true;
	boolean fileMyoFilOffCortex;
	GeneralPath perimeter = new GeneralPath();
	
	
	
	// CONSTRUCTOR *************************************************************************************************
	
	public MyosinHolder (double initX, double initY, double uVectX, double uVectY) {
		super(initX,initY);
				
	}
			
	public void createMyosinHeads(int n){
		
	}

	public MyosinHolder () {
		this(0,0,1,0);
	}
	
	public MyosinHolder (double initX, double initY) {
		this(initX,initY,1,0);
	}

	public void addForce(Point2D f){
		totalTranslationalForce.inc(f);
	}
	
	public void addForce (Point2D f, Point2D loc){
		addForce(f);
		torqArm.getVector(cm,loc);
		addTorque(torqArm,f);
	}
	
	public void addTorque(Point2D r, Point2D f) {
		torque += Point2D.Cross(r,f);
	}
	
	//*************************************************************************************************
	

	public void step (double dT) {
		
		
	}
	
	public void setDiffMultipliers (double transScale, double rotScale) {
		curDPar = transScale*diffPar;
		curDPerp = transScale*diffPerp;
		curDRot = rotScale*diffRot;
	}
	
	public boolean offCortex () {
		return false;
	}
	
	
	public void moveMyosins () {
		
	}
	

		public void reset () {
		// reset forces and velocities
		super.reset();
		totalTranslationalForce.zero();
		torque = 0;
		
		brownianForcePerpendicular.zero();
		brownianForceParallel.zero();
		brownianForceRotational.zero();
		velocity.zero();
		
		deltaTheta = 0;
		angularVel = 0;
	}
	
	public void evaluateProperties () {
		
	}
		
	public void drawYourself (Graphics g, double scale, double [] offset) {
		
	}
	
	public static void startConnectedRecursion () {
		rootMyoHolder.connectedLevel = 0;
		rootMyoHolder.setConnected(rootMyoHolder.connectedLevel);
	}
	
	/** Call "setConnected" for all actin filaments this MyosinHolder contains... those filaments will, in turn
	 *  call "setConnected" for all MyosinHolders they are attached to.  This recursion is broken by setting the
	 *  "isConnected" flag.
	 */
	public void setConnected (int conLevel) {
		connectedLevel = Math.min(conLevel, connectedLevel); //always select most closely connected path as connectedLevel
		setConnectedMyMyosins();
		if (isConnected) return;		// if already connected
		for (int i=0;i<myMyosins.length;i++) {
			if (!myMyosins[i].isFree()) {
				myMyosins[i].boundMon.myFilament.setConnected(connectedLevel+1);
			}
		}
		
		isConnected = true;
	}
	
	public void setConnected (Cluster myCluster) {
		if (clusterTime == Sim2D.simulationTime) return;	// this myosinholder already a member
		myCluster.addElement(this,false);
	
		setToClusterMyMyosins();
		
		for (int i=0;i<myMyosins.length;i++) {
			if (!myMyosins[i].isFree()) {
				myMyosins[i].boundMon.myFilament.setConnected(myCluster);
			}
		}
	}
	
	public void setConnectedMyMyosins () {
		for (int i=0;i<myMyosins.length;i++) {
			myMyosins[i].setConnectedFromHolder();
		}
	}
	
	public void setToClusterMyMyosins () {
		for (int i=0;i<myMyosins.length;i++) {
			myMyosins[i].setToHolderCluster();
		}
	}
	
	public void setToTendrilMyMyosins () {
		for (int i=0;i<myMyosins.length;i++) {
			myMyosins[i].setToHolderTendril();
		}
	}
	
	public void traceTendril (ForceTendril myTendril) {
		if (tendrilTime == Sim2D.simulationTime) return;	// this myosinholder already a member
		myTendril.addElement(this);
	
		setToTendrilMyMyosins();
		
		for (int i=0;i<myMyosins.length;i++) {
			if (!myMyosins[i].isFree()) {
				myMyosins[i].boundMon.myFilament.traceTendril(myTendril);
			}
		}
	}
	
	public void releaseAll() {
		for (int i=0;i < myMyosins.length; i++) {
			myMyosins[i].releaseMon();
		}
	}
	
	public boolean allMyosinsFree() {
		for (int i=0;i < myMyosins.length; i++) {
			if(!myMyosins[i].isFree()) {
				return false;
			}
		}
		return true;
	}
	
	public static void addMyosinHolder (MyosinHolder newHolder) {
		newHolder.myMyoHolderNumber = myoHolderCt;
		theMyoHolders[myoHolderCt] = newHolder;
		myoHolderCt ++;
	}
	
	public void die()  {
		myMyosins = null;
	}
		
	
	public void removeMe () {
		super.removeMe();
		int swapId = myMyoHolderNumber;
		theMyoHolders[swapId] = theMyoHolders[myoHolderCt-1];
		theMyoHolders[swapId].myMyoHolderNumber = swapId;
		myoHolderCt --;
	}
	
	public Myosin getFreeMyosin() {
		for (int i=0;i<myMyosins.length;i++) {
			if (myMyosins[i].isFree()){
				return myMyosins[i];
			}
		}
		return null;
	}
	
	static public void staticInit() {
		Sim2D.addProxy(new AgentProxy() {
			public void reset() {
				for(int i = 0; i < myoHolderCt; i++) {
					theMyoHolders[i].die();
					theMyoHolders[i] = null;
				}
				myoHolderCt = 0;
			}
		});
	}
	

}

