package main;

/*
	A rectangular array of myosin filaments. Based loosely on Verkhovsky AB, Svitkina TM, Borisy GG JCB 1995.
*/

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Random;

import util.*;
import parameters.*;


public class MyosinSurface extends MyosinHolder {
	
	static public String SURFACE_DENSITY = "SurfaceDensity",
							LINE_DENSITY = "LineDensity";
	
	/** array holding all MyosinSurfaces */
	public static MyosinSurface [] theMyoSurfaces = new MyosinSurface [10000];
	
	/** Current number of MyosinSurfaces. */
	static int myoSurfaceCt = 0;

	/** Unique index of this MyosinSurface in the array of all MyosinSurfaces. */
	int myMyoSurfaceNumber;
	
	String densityMode = SURFACE_DENSITY;
	

	/** drag coefficients, needed for calculating diffusion */
	double myosinfilTransGamma,myosinfilPerpGamma,rotGam;
	
	/** For calculating diffusion. */
	double diffPar,diffPerp,diffRot;
	double curDPar,curDPerp,curDRot;  // normally equal to diffPar, etc, but used in step() to increase off cortex diffusivity
	double diffStepPar,diffStepPerp,diffStepRot;
	
	Random generator = new Random(Sim2D.generator.nextLong());
	Random jumpGen = new Random(Sim2D.generator.nextLong());
	
	double xDim = 0;  // dimensions of the myosin surface
	double yDim = 0;
	
	
	// CONSTRUCTOR *************************************************************************************************
	/** specify a myosin surface by center point (initX,initY), dimensions (xDim,yDim), and myosin density */
	public MyosinSurface (double initX, double initY, double xDim, double yDim) {
		super(initX,initY);
		
		this.xDim = xDim;
		this.yDim = yDim;
		this.densityMode = SURFACE_DENSITY;
		double area = xDim*yDim;
		int numMyosToCreate = (int)(Math.ceil(myoDensity*area));
		createMyosinHeads (numMyosToCreate);
		addMyoSurface(this);
	}

	public MyosinSurface () {
		this(0,0,1,1);
	}
	
	public MyosinSurface (double initX, double initY) {
		this(initX,initY,1,1);
	}
	
	public MyosinSurface (Point2D startPt, Point2D stopPt, boolean uniform) {
		double lineDist = Point2D.getDistanceIgnorePeriod(startPt, stopPt);
		this.densityMode = densityMode;
		Point2D myoPt = new Point2D();
		Point2D lineVec = new Point2D();
		lineVec.getVectorIgnorePeriod(startPt, stopPt);
		lineVec.uVect();
		int numMyos = (int)(lineDist*myoLineDensity);
		
		myMyosins = new Myosin[numMyos+1];
		
		if (uniform) {
			double spacing = lineDist/numMyos;
			for (int i=0;i<numMyos;i++) {
				myoPt.add(startPt,i*spacing,lineVec);
				myMyosins[i] = createMyosinAtPoint(myoPt);
			}
			myMyosins[numMyos] = createMyosinAtPoint(stopPt);
		}
		addMyoSurface(this);
	}
	
	
	/** specify a myosin surface with a single myosin at the specified point */
	public MyosinSurface (Point2D myoCenter) {
		myMyosins = new Myosin[1];
		myMyosins[0] = createMyosinAtPoint (myoCenter);
		addMyoSurface(this);
	}
	
	public void createMyosinHeads(int n){
		myMyosins = new Myosin[n];
		Myosin nuMyo;
		for (int i=0;i<n;i++) {
			double rdmX = cm.x + (2*Math.random()-1)*xDim/2.0;
			double rdmY = cm.y + (2*Math.random()-1)*yDim/2.0;
			nuMyo = new Myosin(this);
			nuMyo.cm.set(rdmX,rdmY);
			nuMyo.attachPt.set(rdmX,rdmY);
			nuMyo.bindingSite.set(rdmX,rdmY);
			myMyosins[i] = nuMyo;
		}
	}
	
	public Myosin createMyosinAtPoint (Point2D myoCenter) {
		Myosin nuMyo = new Myosin(this);
		nuMyo.cm.set(myoCenter.x,myoCenter.y);
		nuMyo.attachPt.set(myoCenter.x,myoCenter.y);
		nuMyo.bindingSite.set(myoCenter.x,myoCenter.y);
		return nuMyo;
	}

	public void addForce(Point f){
	}
	
	public void addForce (Point f, Point loc){
	}
	
	public void addTorque(Point r, Point f) {
	}
	
	//*************************************************************************************************
	
	static public void stepAllMyosinSurfaces(double dT) {
		for (int i = 0; i < myoSurfaceCt; i++) {
			if (!theMyoSurfaces[i].removeMe) { theMyoSurfaces[i].step(dT); }
		}
	}

	public void step (double dT) {
		
		moveMyosins (); // all this does is reset head location if not free
		
	}
	
	public void moveMyosins () {
		for (int i=0;i < myMyosins.length; i++) {
			if (myMyosins[i].isFree()) {
				myMyosins[i].bindingSite.set(myMyosins[i].attachPt);
			}
		}
	}

	public void reset () {
	
	}
	
	public void evaluateProperties () {
		
	}
	
	public void drawYourselfFromFile () {
	
	}
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		
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
	
	
	public static void addMyoSurface (MyosinSurface newMyoSurface) {
		newMyoSurface.myMyoSurfaceNumber = myoSurfaceCt;
		theMyoSurfaces[myoSurfaceCt] = newMyoSurface;
		myoSurfaceCt ++;
	}
	
	public void die()  {
		myMyosins = null;
	}
		
	
	public void removeMe () {
		super.removeMe();
		int swapId = myMyoSurfaceNumber;
		theMyoSurfaces[swapId] = theMyoSurfaces[myoSurfaceCt-1];
		theMyoSurfaces[swapId].myMyoSurfaceNumber = swapId;
		myoSurfaceCt --;
	}
	
	public Myosin getFreeMyosin() {
		for (int i=0;i<myMyosins.length;i++) {
			if (myMyosins[i].isFree()){
				return myMyosins[i];
			}
		}
		return null;
	}
	
	static String className = new String("main.MyosinSurface");
	
	static public double myoDensity;
	
	static public double myoLineDensity;
	
	static public void staticInit() {
		
		Sim2D.addProxy(new AgentProxy() {
			public void step(double dT) {
				// MyosinSurfaces update positions of their heads
				stepAllMyosinSurfaces(dT);
			}
			public void reset() {
				for(int i = 0; i < myoSurfaceCt; i++) {
					theMyoSurfaces[i].die();
					theMyoSurfaces[i] = null;
				}
				myoSurfaceCt = 0;
			}
		});

		
		Parameters.addParameterSetListener (
			new ParameterSetListener() {
				public void parametersChanged() {
				}
			}
		);
		Parameters.addParameter(className, "myoSurfaceDensity",0.0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoDensity = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "myoLineDensity",0.0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					myoLineDensity = p.getDoubleValue();
				}
			}
		);
	}

}

