package main;

/*
  It's a crosslinker
 */



import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.util.*;


import util.*;
import collision.*;
import iterators.*;
import gui.*;
import io.*;
import parameters.*;



public class Crosslinker extends Thing{
	
	static Vector XLVariants = new Vector();
	
	XLVariant myXLVariant;
	
	static CollisionDetector myCollisionDetector = new XLCollisionDetector();
	
	
	static Vector bondListeners = new Vector();

	/** An array containing all the crosslinkers */
	static public Crosslinker [] theCrosslinkers = new Crosslinker [20000];
	
	/** Current number of crosslinkers. */
	static public int crosslinkerCt = 0;
	
	/** How far a linker bound to one filament can reach out to bind a second filament. */
	static public double reach;
	
	/** position of this crosslinker in allCrossLinkers */
	public int myCrosslinkerNumber;
	
	/** Whether first end is bound. */
	public Monomer bound1 = null;
	
	/** Whether second end is bound. */
	public Monomer bound2 = null;
	
	double bond1Starts = -1;
	
	double bond2Starts = -1;
	
	double linkStarts = -1;
	
	/** Position of binding site 1. */
	public Point2D site1 = new Point2D();
	
	/** Position of binding site 2. */
	public Point2D site2 = new Point2D();

	/** Current length of doubly bound linker. */
	double length = 0;
	
	/** Unit vector pointing from site1 to site2. */
	public Point2D uVect = new Point2D();
	
	/**  Used each timestep for force calculation. */
	Point2D forceVector = new Point2D();
	
	/** Reused to render cross-linker. */
	Point2D theLine = new Point2D();
	
	/** radius of drawn cross-linker. (*/
	double renderRadius=5;
	
	/** moving average of cross link force, controlled by:
		forceAv(t) = (1-lam)*forceAv(t-1) + lam*force(t) : 0 <= lam <= 1
	*/
	public double forceAv;
	
	/** The crosslink force computed during each step */
	double forceMag;
	
	Random generator = new Random(Sim2D.generator.nextLong());
	
	
	//*******************************************************************************************************
	//CONSTRUCTOR
	
	/** Default constructor. */
	public Crosslinker() {
		this(null);
	}
	
	/** Make a cross-linker and bind it to monomer m. */
	public Crosslinker(Monomer m) {
		super(new Point2D());
		addCrosslinker(this);
		setBound1(m);
	}
	
	/**Make a cross-linker that links two monomers*/
	public Crosslinker(Monomer m1, Monomer m2){
		super(new Point2D());
		addCrosslinker(this);
		setBound1(m1);
		setBound2(m2);
	}
	
	/** Register a CrosslinkListener object that registers whenever
	 * a singly or doubly bound crosslinker releases a bond.
	 */
	static public void addBondListener(CrosslinkListener l) {
		bondListeners.add(l);
	}

	public void miniStep() {
		updateBindingSites();
		doForces();
	}
	
	public Point2D getSite2IgnorePeriod() {
		Point2D p = new Point2D();
		p.add(site1,length,uVect);
		return p;
	}
		
	public void doForces() {

		//calculate the distance between the binding sites;
		length = Point2D.getDistance(site1, site2);
		// calculate forces based on distance between binding sites
		forceVector.getVector(site1, site2);	// make vector between centers
		forceVector.scale(1.0/length);			// convert to unit vector
		uVect.copy(forceVector);				// store this value in uVect...
		forceMag = xlSpringConstant*(length - getRestLength());	// magnitude of repulsion force
		forceVector.scale(-forceMag);
			
		forceAv = (1-Sim2D.forceAverageFraction)*forceAv + Sim2D.forceAverageFraction*forceMag;
				
		if (!Point2D.pointOK(forceVector)) { forceVector.zero(); FileOps.reportln("zeroing xl force"); }
		bound2.addXLForce(forceVector, this, 2);											// force acting on iParticle
		forceVector.scale(-1);
		bound1.addXLForce(forceVector, this, 1);
	}
		

	//*******************************************************************************************************
	//STEP FUNCTION
			
	/** Bound cross-linker calculates spring forces and adds them to attached filament.  All cross-linkers update
	 * their positions relative to attached cross-linkers. */
	public void step (double dT) {
		
		updateBindingSites();
		
		if (isBound()){
			doForces();
			site1.wrapPoint();
			site2.wrapPoint();
			doForceBasedRelease(dT);
		}
		
		else if(site1Bound()) {
			if(Math.random() < xlReleaseProb*dT) {
				release1();
				return;
			}
			uVect.orthogonalVector((bound1.getFilament().uVect));
		}
		else  if(site2Bound()) {
			if(Math.random() < xlReleaseProb*dT) {
				release2();
				return;
			}
			uVect.orthogonalVector((bound2.getFilament().uVect));
		}
		
		reset();
	}
	
	private void doForceBasedRelease(double dT) {
	
		double releaseProb = xlReleaseProb*dT*Math.exp(xlForceBasedReleaseFac*forceAv);
		if(Math.random() < releaseProb) {
			if(Math.random() < 0.5) {
				release1();
			}
			else {
				release2();
			}
		}
	}
	
	public void setConnected (int conLevel, int fromEnd) {
		connectedLevel = Math.min(conLevel, connectedLevel); //always select most closely connected path as connectedLevel
		if (isConnected) return;
		isConnected = true;
		if (fromEnd == 1) {
			bound2.myFilament.setConnected(connectedLevel+1);
		} else {
			bound1.myFilament.setConnected(connectedLevel+1);
		}
	}
	
	public void setConnected (Cluster myCluster, int fromEnd) {
		if (clusterTime == Sim2D.simulationTime) return;	// this xl already a member
		myCluster.addElement(this,false);
		
		if (fromEnd == 1) {
			try {
				bound2.myFilament.setConnected(myCluster);
			} catch (NullPointerException npe) {}
		} else {
			try {
				bound1.myFilament.setConnected(myCluster);
			} catch (NullPointerException npe) {}
		}
	}
	
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) return;
		Graphics2D g2d=(Graphics2D)g;
		
		int rad = (int) (renderRadius*scale);
		int diam = (int) (2*renderRadius*scale);

		int xBindingSite1, yBindingSite1;
		int xBindingSite2, yBindingSite2;
		g2d.setPaint(Color.cyan);
		if(isFree()) {
			//g2d.setPaint(getColorFree());
			if(site2Bound()) {
				xBindingSite2 = site2.getPixX();
				yBindingSite2 = site2.getPixY();
				//g2d.setPaint(getColorBound());
				g2d.drawOval(xBindingSite2-rad, yBindingSite2-rad, diam, diam);
			}
			if(site1Bound()) {
				xBindingSite1 = site1.getPixX();
				yBindingSite1 = site1.getPixY();
				//g2d.setPaint(getColorBound());
				g2d.drawOval(xBindingSite1-rad, yBindingSite1-rad, diam, diam);
			}
		}
		
		else {
			xBindingSite2 = site2.getPixX();
			yBindingSite2 = site2.getPixY();
			xBindingSite1 = site1.getPixX();
			yBindingSite1 = site1.getPixY();
			g2d.setPaint(Color.green);

			// line rendering that works with periodic boundary conditions
			theLine.getVector(site1, site2);
			theLine.uVect();
			theLine.scale(length);
			int theLineX = theLine.getPixX();
			int theLineY = theLine.getPixY();
			g2d.setPaint(getColorRod());
			g2d.drawLine(xBindingSite1, yBindingSite1, xBindingSite1+theLineX, yBindingSite1+theLineY);
			g2d.drawLine(xBindingSite2, yBindingSite2, xBindingSite2-theLineX, yBindingSite2-theLineY);

			g2d.setPaint(getColorBound());
			g2d.drawOval(xBindingSite2-rad, yBindingSite2-rad, diam, diam);
			g2d.drawOval(xBindingSite1-rad, yBindingSite1-rad, diam, diam);
		}
	}
	
	public Point2D getCentroid()  {
		Point2D p = null;
		if(site1!= null && site2!= null) {
			p = new Point2D(0.5*(site1.x+site2.x),0.5*(site1.y+site2.y));
		}
		else if(site1 != null) p = new Point2D(site1);
		else if(site2 != null) p = new Point2D(site2);
		return p;
	}

	public double getCentroidX()  {
	
		if(site1!= null && site2!= null) {
			return 0.5*(site1.x + site2.x);
		}
		else if(site1 != null) return site1.x;
		else if(site2 != null) return site2.x;
		else return -1;
	}
			
	public double getCentroidY()  {
	
		if(site1!= null && site2!= null) {
			return 0.5*(site1.y + site2.y);
		}
		else if(site1 != null) return site1.y;
		else if(site2 != null) return site2.y;
		else return -1;
	}
			
	

	//*************************************************************************************************************
	/** True if site1 is bound. */
	public boolean site1Bound(){
		return bound1 != null;
	}
	
	/** True if site2 is bound. */
	public boolean site2Bound(){
		return bound2 != null;
	}
	
	public Point2D getSinglyBoundSite() {
		return site1Bound() ? site1 : site2Bound() ? site2 : null;
	}
	
	/** Calls releaseMon() for site1. */
	public void release1() {
		releaseMon(bound1);
	}
	
	/** Calls releaseMon() for site2. */
	public void release2() {
		releaseMon(bound2);
	}
	
	/** Release the monomer.  If both sites are now free, schedule for removal. */
	public void releaseMon (Monomer mon) {
		if (bound1 == mon){
			bound1.boundLinker = null;
			bound1 = null;
			if(bound2 == null) {  // singly bound
				setForRemove();
				for(Enumeration en = bondListeners.elements(); en.hasMoreElements(); ) {
					((CrosslinkListener)en.nextElement()).singlyBoundLinkerReleased(Sim2D.simulationTime-bond1Starts);
				}
				bond1Starts = -1;
			}
			else { // doubly bound
				for(Enumeration en = bondListeners.elements(); en.hasMoreElements(); ) {
					((CrosslinkListener)en.nextElement()).doublyBoundLinkerReleased(Sim2D.simulationTime-bond1Starts);
				}
				linkStarts = -1;
			}
			return;
		}
		else if (bound2 == mon){
			bound2.boundLinker = null;
			bound2 = null;
			if(bound1 == null) { // singly bound
				setForRemove();
				for(Enumeration en = bondListeners.elements(); en.hasMoreElements(); ) {
					((CrosslinkListener)en.nextElement()).singlyBoundLinkerReleased(Sim2D.simulationTime-bond1Starts);
				}
				bond2Starts = -1;
			}
			else { // doubly bound
				for(Enumeration en = bondListeners.elements(); en.hasMoreElements(); ) {
					((CrosslinkListener)en.nextElement()).doublyBoundLinkerReleased(Sim2D.simulationTime-bond1Starts);
				}
				linkStarts = -1;
			}
			return;
		}
		else {
			FileOps.reportln("tried to release a monomer = " + mon + " that doesn;t exist! My bound1 = " + bound1 + " and bound2 = " + bound2);
		}
	}
	
	/** sets monomer as bound monomer for site1 and sets monomer.boundLinker to this. */
 	public void setBound1 (Monomer boundMonomer){
		if(bound1 != null) {
			bound1.boundLinker = null;
		}
 		bound1 = boundMonomer;
		if(bound1 != null) {
			bound1.boundLinker = this;
			bond1Starts = Sim2D.simulationTime;
			if(bound2 != null) { // now doubly bound
				linkStarts = Sim2D.simulationTime;
			}
		}
	}

	/** sets monomer as bound monomer for site2 and sets monomer.boundLinker to this. */
 	public void setBound2 (Monomer boundMonomer){
 		
 		
 		if(bound2 != null) {
			bound2.boundLinker = null;
		}
 		bound2 = boundMonomer;
		if(bound2 != null) {
			bound2.boundLinker = this;
			bond2Starts = Sim2D.simulationTime;
			if(bound1 != null) { // now doubly bound
				linkStarts = Sim2D.simulationTime;
			}
		}
	}
	
	public void updateBindingSites() {
		if(site1Bound()) {
			site1.set(bound1.getLocation());
		}
		
		if(site2Bound()){
			site2.set(bound2.getLocation());
		}
	}
			
	public double getInternalForce(){
		if (site1Bound() && site2Bound()){
			return (length-getRestLength())*xlSpringConstant;
		}
		else {
			return 0;
		}
	}
	
	public double getLength() {
		return length;
	}
	
	/** True if either sites unbound. */
	public boolean isFree() {
		return bound1 == null || bound2 == null;
	}
 	
	/** True iff both sites are bound */
 	public boolean isBound(){
 		return bound1 != null && bound2 != null;
 	}
 	
	
	/** Add crosslinker to list of all linkers. */
 	public static void addCrosslinker (Crosslinker newCrosslinker) {
		newCrosslinker.myCrosslinkerNumber = crosslinkerCt;
		theCrosslinkers[crosslinkerCt] = newCrosslinker;
		crosslinkerCt++;
	}
	
	/** True if this crosslink is bound to filament. */
 	public boolean isBoundTo(Actin actin) {
 		boolean isBound = false;
 		if (bound1 != null){
 			isBound = bound1.myFilament.myActinNumber == actin.myActinNumber;
 		}
 		if (!isBound){
 			if (bound2 != null){
 				isBound = bound2.myFilament.myActinNumber == actin.myActinNumber;
 			}
 		}
 		return isBound;
 	}
	
	public void die() {
		bound1 = bound2 = null;
	}
 	
	/** Remove this linker from list of all linkers. */
	public void removeMe() {
		super.removeMe();
		int swapId = myCrosslinkerNumber;
		theCrosslinkers[swapId] = theCrosslinkers[crosslinkerCt-1];
		theCrosslinkers[swapId].myCrosslinkerNumber = swapId;
		crosslinkerCt --;
	}

	/*************  DIAGNOSTICS *********************/
	static public void checkGhostBonds(int ct) {
		Crosslinker cl;
		Monomer m;
		for(int i = 0; i < crosslinkerCt; i++) {
			cl = theCrosslinkers[i];
			if(cl.bound1 != null) {
				m = cl.bound1;
				if(m.boundLinker != cl) {
					FileOps.reportln(ct + "   crosslinker = " + cl.myCrosslinkerNumber + " binds unreciprocating monomer = " + m.myMonNumber);
				}
			}
			else if(cl.bound2 != null) {
				m = cl.bound2;
				if(m.boundLinker != cl) {
					FileOps.reportln(ct + "   crosslinker = " + cl.myCrosslinkerNumber + " binds unreciprocating monomer = " + m.myMonNumber);
				}
			}
		}
	}
	
	static public XLVariant getXLVariant(String name) throws Exception {
		if(name == null) {
			if(XLVariants.size() > 1) {
				throw new Exception("There are multiple Crosslinker variants; you must specify one by name");
			}
			if(XLVariants.isEmpty()) {
				XLVariant v = new XLVariant("default");
				XLVariants.add(v);
				return v;
			}
			else return (XLVariant) XLVariants.elementAt(0);
		}
		
		XLVariant v;
			
		for(Enumeration en = XLVariants.elements(); en.hasMoreElements(); ) {
			v = (XLVariant) en.nextElement();
			if(name.equals(v.getName())) {
				return v;
			}
		}
		throw new Exception("Crosslinker.getXLVariant(): there is no XLVariant named: " + name);
	}
			
	
	public XLVariant getXLVariant() {
		
		return myXLVariant;
	}
	
	public void setXLVariant(XLVariant v) {
		myXLVariant = v;
	}
	
	static public void  loadXLVariant(AMInputStream in)	throws Exception {
		String name;
		try {
			name = in.nextString();
		} catch(Exception e) {
			System.out.println("Crosslinker.loadXLVariant():XLVariant tag must be followed by a variant name");
			throw(e);
		}
			
		try {
			XLVariant v = new XLVariant(name);
			v.load(in);
			XLVariants.add(v);
		} catch(Exception e) {
			System.out.println("Problem loading XLVariant " + e.toString());
			throw(e);
		}
	}
	
	// CROSSLINKER BINDING ***** CROSSLINKER BINDING ***** CROSSLINKER BINDING ***** CROSSLINKER BINDING *****
	
	/** Loop through all filaments and allow each to bind linkers from the pool. */
	static public void doCrossLinkerBinding(double dT) {
		double [] probDist = Actin.getWeightedMonomerChoices();
		for(Enumeration en = XLVariants.elements(); en.hasMoreElements(); ) {
			((XLVariant)en.nextElement()).doCrossLinkerBinding(probDist,dT);
		}
	}
	

	static public void checkCrossBinding(int ct)  {
		Crosslinker cl;
		Monomer m;
		for(int j = 0; j < crosslinkerCt; j++) {
			cl = theCrosslinkers[j];
			if(cl.bound1 != null) {
				m = cl.bound1;
				for(int i = 0; i < crosslinkerCt; i++) {
					if(i!=j && theCrosslinkers[i].bound1 == m) {
						FileOps.reportln(ct + "   crosslinker = " + cl.myCrosslinkerNumber + " binds monomer = " + m.myMonNumber + " orphan-bound to linker = " + theCrosslinkers[i].myCrosslinkerNumber);
					}
					if(i!=j && theCrosslinkers[i].bound2 == m) {
						FileOps.reportln(ct + "   crosslinker = " + cl.myCrosslinkerNumber + " binds monomer = " + m.myMonNumber + " orphan-bound to linker = " + theCrosslinkers[i].myCrosslinkerNumber);
					}
				}
			}
			else if(cl.bound2 != null) {
				m = cl.bound2;
				for(int i = 0; i < crosslinkerCt; i++) {
					if(i!=j && theCrosslinkers[i].bound1 == m) {
						FileOps.reportln(ct + "   crosslinker = " + cl.myCrosslinkerNumber + " binds monomer = " + m.myMonNumber + " orphan-bound to linker = " + theCrosslinkers[i].myCrosslinkerNumber);
					}
					if(i!=j && theCrosslinkers[i].bound2 == m) {
						FileOps.reportln(ct + "   crosslinker = " + cl.myCrosslinkerNumber + " binds monomer = " + m.myMonNumber + " orphan-bound to linker = " + theCrosslinkers[i].myCrosslinkerNumber);
					}
				}
			}
		}
	}
	
	static public void linkerReport(int id) {
		Crosslinker cl;
		for (int i = 0; i < crosslinkerCt; i++) {
			cl = theCrosslinkers[i];
			if(cl.site1Bound() && cl.bound1.myFilIndex < cl.bound1.myFilament.indexOffset) {
				FileOps.reportln("Crosslinker ID = " + id);
				FileOps.report(" i = " + cl.myCrosslinkerNumber);
				FileOps.reportln(" site 1: # = " + cl.bound1.myMonNumber + " index = " + cl.bound1.myFilIndex + " offset = " + cl.bound1.myFilament.indexOffset);
			}
			if(cl.site2Bound() && cl.bound2.myFilIndex < cl.bound2.myFilament.indexOffset) {
 				FileOps.report("Crosslinker ID = " + id);
				FileOps.report(" i = " + cl.myCrosslinkerNumber);
				FileOps.reportln(" site 1: # = " + cl.bound2.myMonNumber + " index = " + cl.bound2.myFilIndex + " offset = " + cl.bound2.myFilament.indexOffset);
			}
		}
	}
	
		
		
		/******* CROSSLINKER STATIC CONSTANT PARAMETERS *********/
	
	/** How far a linker bound to one filament can reach out to bind a second filament. */
	static public double xlMaxBindingDistance;
	

	/** The maximaum allowable interval between crosslinker binding checks.  This
	 * limits the interval by adaptive mechanism.. */
	static public double xlBindCheckMaxInterval;
	
	/** For adaptive adjustment of intervals at which we compute crosslinker bond formation.
	 * Defines the TARGET max probability of forming a new crosslink in the adaptive
	 * interval calculated over all crosslinks. */
	static public double xlTargetBindingProbPerTimestep;
	
	/** For adaptive adjustment of intervals at which we compute crosslinker bond formation.
	 * Defines the max ALLOWABLE probability of forming a new crosslink in the adaptive
	 * interval calculated over all crosslinks. */
	static public double xlMaxBindingProbPerTimestep;
	
	/** For adaptive adjustment of intervals at which we compute crosslinker bond formation.
	 * Defines the TARGET move for a filament monomer during
	 * one adaptive timestep interval
	 **/
	static public double targetMonomerMoveFraction;
	
	/** For adaptive adjustment of intervals at which we compute crosslinker bond formation.
	 * Defines the MAX move allowed for a filament monomer during
	 * one adaptive timestep interval
	 **/
	static public double maxMonomerMoveFraction;
		
	/** Total number of linkers available to bind filaments. */
	static public double totalLinkerPool;
	

	/******* CROSSLINKER PARAMETERS THAT CAN DIFFER ACROSS XLVARIANTS *********/

	/** multiplier for force-based crosslinker release. */
	static private double xlForceBasedReleaseFac;
	
	/** linkers only bind filament if contact angle exceeds this angle. */
	static private double xlMinBindingAngle;
	
	/** Maximum stretch before a linker releases. */
	static private double xlMaxStretch;

	/** Rest length for force calculations. */
	static private double xlRestLength;
	
	/** Stiffness of a crosslinker. */
	static private double xlSpringConstant;	// pN/nm

	/** Probabality that a bound linker will release a bond. */
	static private double xlReleaseProb;
	
	/** Probability a soluble crosslinker will bind a filament from bulk solution. */
	static private double xlRecruitmentProb;
	
	/** Probability a singly-bound linker will bind a second filament if it is within reach. */
	static private double xlBindingProb;
	
	static private Color ColorRod;
	static private Color ColorFree;
	static private Color ColorBound;
	
	static String className = new String("main.Crosslinker");
	
	/*********  access methods for properties that vary in XLVariants. ******/
	
	
	public double getRestLength() {
		if(myXLVariant != null && myXLVariant.hasRestLength()) {
			return myXLVariant.getRestLength();
		}
		return xlRestLength;
	}

	
	public double getSpringConstant() {
		if(myXLVariant != null && myXLVariant.hasSpringConstant()) {
			return myXLVariant.getSpringConstant();
		}
		return xlSpringConstant;
	}

	
	public double getRecruitmentProb() {
		if(myXLVariant != null && myXLVariant.hasRecruitmentProb()) {
			return myXLVariant.getRecruitmentProb();
		}
		return xlRecruitmentProb;
	}
	public double getBindingProb() {
		if(myXLVariant != null && myXLVariant.hasBindingProb()) {
			return myXLVariant.getBindingProb();
		}
		return xlBindingProb;
	}
	public double getReleaseProb() {
		if(myXLVariant != null && myXLVariant.hasReleaseProb()) {
			return myXLVariant.getReleaseProb();
		}
		return xlReleaseProb;
	}

	public double getForceBasedReleaseFac() {
		if(myXLVariant != null && myXLVariant.hasForceBasedReleaseFac()) {
			return myXLVariant.getForceBasedReleaseFac();
		}
		return xlForceBasedReleaseFac;
	}
	public double getMinBindingAngle() {
		if(myXLVariant != null && myXLVariant.hasMinBindingAngle()) {
			return myXLVariant.getMinBindingAngle();
		}
		return xlMinBindingAngle;
	}
	public double getXLMaxStretch() {
		if(myXLVariant != null && myXLVariant.hasMaxStretch()) {
			return myXLVariant.getMaxStretch();
		}
		return xlMaxStretch;
	}
	public Color getColorRod() {
		if(myXLVariant != null && myXLVariant.hasColorRod()) {
			return myXLVariant.getColorRod();
		}
		return ColorRod;
	}
	public Color getColorFree() {
		if(myXLVariant != null && myXLVariant.hasColorFree()) {
			return myXLVariant.getColorFree();
		}
		return ColorFree;
	}
	public Color getColorBound() {
		if(myXLVariant != null && myXLVariant.hasColorBound()) {
			return myXLVariant.getColorBound();
		}
		return ColorBound;
	}
	
	static public void staticInit() {
		
	
		Parameters.addParameterSetListener (
			new ParameterSetListener() {
				public void parametersChanged() {
					reach = xlMaxBindingDistance;
				}
			}
		);
		Parameters.addParameter(className, "xlMaxBindingDistance",30,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlMaxBindingDistance = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlTargetBindingProbPerTimestep",0.2,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlTargetBindingProbPerTimestep = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlMaxBindingProbPerTimestep",0.25,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlMaxBindingProbPerTimestep = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlBindCheckMaxInterval",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlBindCheckMaxInterval = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "targetMonomerMoveFraction",0.5,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					targetMonomerMoveFraction = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "maxMonomerMoveFraction",0.6,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					maxMonomerMoveFraction = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlRestLength",20,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlRestLength = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlSpringConstant",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlSpringConstant = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlReleaseProb",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlReleaseProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlRecruitmentProb",0.001,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlRecruitmentProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlBindingProb",500,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlBindingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlMaxStretch",30,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlMaxStretch = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "totalLinkerPool",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					totalLinkerPool = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlMinBindingAngle",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlMinBindingAngle = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "xlForceBasedReleaseFac",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xlForceBasedReleaseFac = p.getDoubleValue();
				}
			}
		);
		// modified Parameter.java to include getPaintValue
		Parameters.addParameter(className, "ColorRod","Color.cyan",
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						ColorRod = p.getPaintValue();
					}
				}
		);
		Parameters.addParameter(className, "ColorFree","Color.cyan",
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						ColorFree = p.getPaintValue();
					}
				}
		);
		Parameters.addParameter(className, "ColorBound","Color.cyan",
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						ColorBound = p.getPaintValue();
					}
				}
		);
		
	
		Sim2D.addProxy(new AgentProxy() {
			public void checkCollisions(double simTime) {
				myCollisionDetector.checkCollisions(simTime);
			}
			public void init() {
				myCollisionDetector.init();
				if(XLVariants.isEmpty()) {
					XLVariants.add(new XLVariant("default"));
				}
				for(Enumeration en = XLVariants.elements(); en.hasMoreElements(); ) {
					((XLVariant)en.nextElement()).init(totalLinkerPool,xlRestLength,xlSpringConstant,xlBindingProb,xlRecruitmentProb);
	
				}
			}
			public void doForces(double dT) {
				for (int i = 0; i < crosslinkerCt; i++) {
					if (!theCrosslinkers[i].removeMe) { theCrosslinkers[i].step(dT); }
				}
			}
			public void step(double dT) {
				// Crosslinkers bind actin filaments from bulk pool..
				doCrossLinkerBinding(dT);
			}
				
		});
	}

}
	

