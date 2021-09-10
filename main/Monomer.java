package main;

/*
	a spherical particle with Brownian diffusion and elastic collisions with other particles
*/

import java.awt.*;
import java.util.*;

import analysis.ValueTracker;

import util.*;
import analysis.*;

public class Monomer extends Thing implements ForceSubstrate {
	
	/** List of all monomers. */
	static public Monomer[] theMonomers=new Monomer[1000000];
	
	/** Total number of monomers */
	static public int monomerCt = 0;
	
	/** radius of a monomer. */
	static double radius = Actin.actinRadius;
	
	/** ID of this monomer. */
	int myMonNumber;
	
	/** True if monomers should be explicitly drawn */
	static public boolean drawMonomers=false;
	
	/** The filament this monomer belongs to. */
	public Actin myFilament = null;
	
	/** integer location of this monomer in the linked list.. need for fast position updating */
	int myFilIndex;
	
	/** The next monomer in the linked list. */
	public Monomer prev = null;

	/** The previous monomer in the linked list. */
	public Monomer next = null;
	
	/** A Link attached to this Monomer */
	public Link theLink = null;
	
	/** The crosslinker bound to this monomer. */
	public Crosslinker boundLinker = null;
	
	/** The myosin head bound to this monomer. */
	public Myosin boundMyo = null;
	
	/** Flag declaring this monomer as a myosin binding site (or not) */
	public boolean isBindingSite = false;
	
	/** The force applied by this monomer on its filament at the last time-step */
	Point2D myForce = new Point2D();
	
	/** Track the force-state (compression/tension) of this monomer */
	// **** Note... relate the time-scale of this ValueTracker to the buckling time-constant
	public double avgMonStress = 0;
	
	public double buckleProgress = 0;
	
	int monStressSeg = 0;	// marks the stress segment this mon is in... can change every time-step
	
	/** The myosin V head bound to this monomer. */
	MyosinV boundMyoV = null;

	public double birthTime;
	
	int locationLastUpdate = -1;
	
	//*****************************************************************************************

	static public double getMonomerDensity()  {
		double denom = Sim2D.xDimension*Sim2D.yDimension*Sim2D.zDimension*Constants.avagadro;
		return monomerCt*1e30/denom;
	}
			
	
	/** Make a monomer at (0,0). */
	public Monomer (){
		this(0,0);
	}
	
	/** Make a monomer at (initX,initY). */
	public Monomer (double initX, double initY) {
		super(initX,initY);
		birthTime = Sim2D.simulationTime;
		addMonomer(this);
	}
	
	//*****************************************************************************************
	/** Set the actin filament this monomer belongs to and the index of this monomer along filament. */
	public void setActinFilament (Actin fil, int ind){
		myFilament = fil;
		myFilIndex = ind;
	}
	
	/** Set the actin filament this monomer belongs to, but not the filament index.. */
	public void setActinFilament (Actin fil){
		myFilament = fil;
	}

	public void drawYourself (Graphics g, double scale, double [] offset) {

		g.setColor(Color.white);
				
		int xCenter = (int) ((cm.x-offset[0])*scale);
		int yCenter = (int) ((cm.y-offset[1])*scale);
		int pixelRadius = (int) (radius*scale);
		int pixelDiameter = (int) (2*radius*scale);
		
		g.fillOval(xCenter-pixelRadius,yCenter-pixelRadius,pixelDiameter,pixelDiameter);
		
	}
	
	/** Add a force to this monomer (which is transmitted to its filament */
	public void addForce(Point2D tempForce) {
		myForce.copy(tempForce);
		myFilament.addForce(tempForce, cm);
	}
	
	/** Part of ForceSubstrate implementation */
	public void addForce(Point2D p, Point2D tempForce) {
		addForce(tempForce);
	}

	
	public double registerStress(double s, double buckleF, double buckleR) {
		avgMonStress = (1-Sim2D.forceAverageFraction)*avgMonStress + Sim2D.forceAverageFraction*s;
		if(buckleF > 1) {
			buckleProgress = buckleProgress*(1-buckleR) + buckleR;
			if(buckleProgress > 1) {
				//System.out.println("mon# = " + myMonNumber + " fil# = " + myFilament.myActinNumber + " buckleProgress = " + buckleProgress);
				buckleProgress = 1;
			}
		}
		else {
			buckleProgress = buckleProgress*(1-buckleR);
			if(buckleProgress < 0) buckleProgress = 0;
		}
		return buckleProgress;
	}

	
	/** Add a force specifically from a crosslinker... used for tracking crosslinkers on a filament, etc */
	public void addXLForce(Point2D tempForce, Crosslinker xl, int endBound) {
		myForce.copy(tempForce);
		myFilament.addXLForce(tempForce, cm, xl, endBound);
	}
	
	/** Add a force specifically from a myosin... used for counting myosins on a filament, etc */
	public void addMyoForce(Point2D tempForce, Myosin myo) {
		myForce.copy(tempForce);
		myFilament.addMyoForce(tempForce, cm, myo);
	}
	
	/** Get this monomers filament */
	public Actin getFilament() {
		return myFilament;
	}
	
	/** Get the location of this monomer */
	public Point2D getLocation(){
		if (locationLastUpdate == Sim2D.counter) { return cm; } // so this costly update will only happen once per time-step, even if called multiple times
		myFilament.setMonomerPosition(this);
		locationLastUpdate = Sim2D.counter;
		return cm;
	}
	
	/** True iff this monomer has no binding partners. */
	public boolean isFree(){
		return boundLinker == null && boundMyo == null && boundMyoV == null && theLink == null;
	}
	
	/** True iff this monomer is not bound to a crosslinker. */
	public boolean isFreeLinker (){
		return boundLinker == null;
	}
	
	/** This monomer releases its bound linker by calling boundLinker.releaseMon(). */
	public void releaseLinker () {
		if(boundLinker!=null){
			boundLinker.releaseMon(this);
			boundLinker = null;
		}
	}
	
	/** True iff this monomer not currently bound to a myosin. */
	public boolean isFreeMyosin (){
		return boundMyo == null;
	}
		
	/** Bind a myosin and set its boundMon to this monomer */
	public void bindMyo (Myosin myo) {
		if (boundMyo != null) { return; }
		boundMyo = myo;
		boundMyo.boundMon = this;
	}
		
	/** Release the bound a myosin and set its boundMon to this */
	public void releaseMyo () {
		if(boundMyo==null) { return; }
		boundMyo.boundMon = null;
		boundMyo.releaseMon();
		boundMyo = null;
		myForce.zero();
	}
	
	/** True iff this monomer not currently bound to a myosin. */
	public boolean isFreeLink (){
		return theLink == null;
	}
	
	/** Bind a myosin and set its boundMon to this monomer */
	public void bindLink (Link l) {
		if (theLink != null) { return; }
		theLink = l;
	}
	
	/** Release the bound a myosin and set its boundMon to this */
	public void releaseLink () {
		if(theLink==null) { return; }
		theLink.breakLink();
		theLink = null;
	}
	
	public void linkDetached(Link l) {
		theLink = null;
	}
	
	public Point2D getGlobalCoords(double l) {
		return getLocation();
	}
	
	/** True iff this monomer not currently bound to a myosinV. */
	public boolean isFreeMyosinV (){
		return boundMyoV == null;
	}
	
	/** Bind a myosinV and set its boundMon to this monomer */
	public void bindMyoV (MyosinV myoV) {
		if (boundMyoV != null) {
			return;
		}
		boundMyoV = myoV;
		boundMyoV.boundMon = this;
	}
	
	/** Release the bound a myosinV and set its boundMon to this */
	public void releaseMyoV () {
		if(boundMyoV==null) { return; }
		boundMyoV.boundMon = null;
		boundMyoV = null;
	}
	
	/** Release all bound entities. */
	public void releaseAll () {
		releaseLinker();
		releaseMyo();
		releaseMyoV();
		releaseLink();
	}
	
	public static void reportBoundMonomers () {
		int boundMyos = 0;
		int boundMyoVs = 0;
		int boundLinkers = 0;
		int boundMyosWithNullLink = 0;
		int boundMyosWithBadLink = 0;
		int boundMyoVsWithNullLink = 0;
		int boundMyoVsWithBadLink = 0;
		for (int i=0;i<monomerCt;i++) {
			if (!theMonomers[i].isFreeLinker()) { boundLinkers++; }
			if (!theMonomers[i].isFreeMyosin()) { boundMyos++; }
			if (!theMonomers[i].isFreeMyosinV()) { boundMyoVs++; }
			if (!theMonomers[i].isFreeMyosin()) {
				if (theMonomers[i].boundMyo.boundMon != null) {
					if (theMonomers[i].boundMyo.boundMon != theMonomers[i]) {
						boundMyosWithBadLink++;
						System.out.println("bad link" + theMonomers[i].boundMyo);
					}
				} else { boundMyosWithNullLink++; }
			}
			if (!theMonomers[i].isFreeMyosinV()) {
				if (theMonomers[i].boundMyoV.boundMon != null) {
					if (theMonomers[i].boundMyoV.boundMon != theMonomers[i]) {
						boundMyoVsWithBadLink++;
						System.out.println("bad link" + theMonomers[i].boundMyoV);
					}
				} else { boundMyoVsWithNullLink++; }
			}
		}
		System.out.println("Monomers bound to linkers = " + boundLinkers);
		System.out.println("Monomers bound to myos = " + boundMyos);
		System.out.println("Monomers bound to myos with null link back = " + boundMyosWithNullLink);
		System.out.println("Monomers bound to myos with bad link back = " + boundMyosWithBadLink);
		System.out.println("Monomers bound to myoVs = " + boundMyoVs);
		System.out.println("Monomers bound to myoVs with null link back = " + boundMyoVsWithNullLink);
		System.out.println("Monomers bound to myoVs with bad link back = " + boundMyoVsWithBadLink);
	}
	
	/** Add a monomer to list of all monomers. */
	public static void addMonomer (Monomer newMonomer) {
		newMonomer.myMonNumber = monomerCt;
		theMonomers[monomerCt] = newMonomer;
		monomerCt ++;
	}
	
	/** This monomer releases all references  to prepare for death.No need to call remove me:  It
	 * will be removed by the filament.
	 * */
	public void die() {
		releaseAll();
		prev = next = null;
		myFilament = null;
	}
	
	/** Remove this monomer from list of all monomers. */
	public void removeMe () {
		super.removeMe();
		int swapId = myMonNumber;
		theMonomers[swapId] = theMonomers[monomerCt-1];
		theMonomers[swapId].myMonNumber = swapId;
		monomerCt --;
	}
	
	
	static void staticInit() {
	
		Sim2D.addProxy(new AgentProxy() {
			/** Reset everything in before starting next simulation run. */
			public void reset() {
				for(int i =0; i < monomerCt; i++) {
					theMonomers[i].die();
					theMonomers[i] = null;
				}
				monomerCt = 0;
			}
		});

	}
	
	
	Vector markers = new Vector();
	
	public void addMarker(String owner) {
		markers.add(new Mark(owner));
	}
	
	public void clearMarks() {
		markers.clear();
	}
	
	public void setMark(String owner) {
		Mark m;
		for(Enumeration en = markers.elements(); en.hasMoreElements(); ) {
			m = (Mark) en.nextElement();
			if(owner.equals(m.getOwner())) {
				m.setP(getLocation());
			}
		}
	}
	
	public Point2D getMark(String owner) {
		Mark m;
		for(Enumeration en = markers.elements(); en.hasMoreElements(); ) {
			m = (Mark) en.nextElement();
			if(owner.equals(m.getOwner())) {
				return (m.getP());
			}
		}
		return null;
	}

}

