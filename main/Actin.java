package main;


/*
	A filamentous particle with Brownian diffusion and elastic collisions with other particles
 */

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.text.DecimalFormat;
import java.util.*;
import java.awt.Graphics2D;
import sun.tools.tree.ThisExpression;

import analysis.*;
import util.*;
import collision.*;
import io.*;
import iterators.GlidingAssayEvaluator;
import gui.*;
import parameters.*;

/** An actin filament.  */
public class Actin extends RigidRod {
	
	static Vector severCandidates = new Vector();

	static ActinListener monomerListener;

	/** The list of all actin filaments. */
	static public Actin [] theActins = new Actin [50000];

	/** Current number of actin filaments. */
	static public int actinCt = 0;

	/** Barbed end monomer addition probability */
	static double actinBOnProb;
	
	/** Barbed end monomer removal probability */
	static double actinBOffProb;
	
	/** Pointed end monomer addition probability */
	static double actinPOnProb;

	/** Pointed end monomer removal probability */
	static double actinPOffProb;
	
	/** Current concentration of free monomers ((total monomer - polymerized monomer)/volume) */
	static double freeMonomerConc;
		
	static double numberOfMonomers;
	
	/** Length of one monomer in nm. */
	static public double monLength = 2.7;
	
	/** Normal actin color. */
	Color actinColor = Color.red;
	
	/** faded color for connectoivoty drawing. */
	Color fadedColor = new Color (0.1f,0.1f,0.1f);
	
	public int forceTextLocator = 1;  // 1 for above ends, -1 for below
	
	/** radius of an actin filament. */
	static public double actinRadius = 2.0;

	/** The position of this filament in the array of all filaments. */
	public int myActinNumber;
	
	/**for use when simulating polar or apolar bundles*/
	public double initPolarity;
	
	/** True if this actin does not add or lose monomers. */
	boolean isStatic = false;
		
	FilamentAttachment pAttachment = null;
	FilamentAttachment bAttachment = null;
	
	static DecimalFormat pinForceFormat = new DecimalFormat ("#0.00#; #0.00#");

	Point2D aLoad = new Point2D();  // resistive load
		
	InternalStressTracker myStressTracker = null;
		

	/** Connectivity tracking */
	public double connectTime = 0;
	public int myoHolderConnectCt = 0;
	public MyosinHolder [] myoHolderConnects;
	
	/** Crosslinker tracking */
	public double xlTime = 0;
	public Crosslinker [] xlList;
	public int [] xlEndBoundList;
	public int xlCt = 0;
	
	/** Myosin tracking */
	public double myoTime = 0;
	public Myosin [] myoList;
	public int myoCt = 0;
	
	/** Myosin residency */
	int sumMyos = 0;
	double startTime = 0;
	double finishTime = 0;
	
	/** Total length of this filament in nm = nMonomers**/
	public double length;
	
	/** Default filament length in monomers. */
	static int defaultMonomerLength = 200;
	
	/** final length in monomers to which this filament will grow */
	int finalMonomerLength=defaultMonomerLength;
	
	/** The maximum length in monomers allowable for a filament (to avoid confusion with periodic boundary conditions. */
	static public int maxMon;
	
	/** The minimum length of a new filament in monomers. */
	static public int minMon = 3;
	
	/** Current number of monomers for this filament. */
	public int nMonomers=0;
		
	/** Used to offset indexed position of each monomer along filament as pointed end monomers disappear. */
	int indexOffset = 0;
	
	/** Position of barbed end. */
	Point2D bEnd;
	
	/** Position of pointed end. */
	Point2D pEnd;
		
	/** Monomer at pointed end. */
	public Monomer pEndMonomer;

	/** Monomer at barbed end. */
	public Monomer bEndMonomer;
	
	int severingPlace = 0;

	boolean capped = false;
	
	/* used for writing position change in sliding filament assay */
	
	/** Random number generator to generate initial filament lengths. */
	static Random lengthGen = new Random((long)Math.random());
	Color myColor;
	
	/* force analysis in gliding assay */
	static double posMyoForce = 0;
	static double negMyoForce = 0;
	
	static double avgFilamentLifetime = 0;
	
	static DecimalFormat lengthFormat = new DecimalFormat ("#000.#; #000.#");	// a clean format nm lengths
	static DecimalFormat forceFormat = new DecimalFormat ("#0.000#; #0.000#");	// a clean format nm lengths


	// CONSTRUCTORS ***** CONSTRUCTORS ***** CONSTRUCTORS ***** CONSTRUCTORS***** CONSTRUCTORS

	/** Make a new filament at (0,0), angle = 0, length = defaultMonomerLength, not static, not Pinned */
	public Actin () {
		this(0,0,0,defaultMonomerLength,defaultMonomerLength,false,null,null);
	}
	
	/**
	 * Constructor
	 *
	 * Make a new filament at (initX,initY), angle = 0, length = defaultMonomerLength, not static, not Pinned
	 *
	 * @param    initX               a  double
	 * @param    initY               a  double
	 *
	 */
	public Actin (double initX, double initY) {
		this(initX,initY,0, defaultMonomerLength,defaultMonomerLength,false,null,null);
	}
	
	/**
	 * Constructor
	 * Make a new filament that is not static and not pinned
	 *
	 * @param    initX               a  double  = initial X position
	 * @param    initY               a  double = initial Y position
	 * @param    initAng             a  double = initial angle
	 * @param    numMons             an int  = initial number of monomers
	 * @param    finalMons           an int = final number of monomers before capping
	 *
	 */
	public Actin (double initX, double initY, double initAng, int numMons, int finalMons) {
		this(initX,initY,initAng, numMons,finalMons,false,null,null);
	}

	/**
	 * Constructor
	 *
	 * Make a new filament with the specified monomer as its barbed end.
	 *
	 * @param    newBEndMonomer      a  Monomer
	 *
	 */
	public Actin (Monomer newB, Actin parent) {
		super(0.0,0.0,parent.uVect,actinRadius);
		pEnd = end1; bEnd = end2;
		indexOffset = parent.indexOffset;
		isStatic = parent.isStatic;
		bEndMonomer = newB;
		bEndMonomer.next = null;
		pEndMonomer = parent.pEndMonomer;
		Monomer myMonomer = pEndMonomer;
		while (myMonomer != null){
			myMonomer.setActinFilament(this);
			myMonomer = myMonomer.next;
		}
		
		pEnd.set(pEndMonomer.cm);
		bEnd.set( bEndMonomer.cm);
		uVect.set(parent.uVect);
		upVect.orthogonalVector(uVect);
		refUPVect.set(upVect);
		countMonomers();
		lengthFromMonomerCount(nMonomers);
		cm.add(pEnd,length/2,uVect);
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
		evaluateProperties();
	
		capped = false;
		addActin(this);
		initForceTracking();
		setActinColor();
	}
	
	/**
	 * Constructor
	 *
	 * @param    initX               a  double  = initial X position
	 * @param    initY               a  double = initial Y position
	 * @param    initAng             a  double = initial angle
	 * @param    numMons             an int  = initial number of monomers
	 * @param    finalMons           an int = final number of monomers before capping
	 * @param    isStatic            a  boolean = does this filament grow and shrink
	 * @param    pinned              a  boolean = does this filament move?
	 *
	 */
	public Actin (double initX, double initY, double initAng, int numMons, int finalMons, boolean isStatic, FilamentAttachment minusAttachment, FilamentAttachment plusAttachment) {
		super (initX,initY,initAng,actinRadius);
		this.isStatic = isStatic;
		this.pAttachment = minusAttachment;
		this.bAttachment = plusAttachment;

		pEnd = end1; bEnd = end2;
					
		pEndMonomer = new Monomer();
		bEndMonomer = new Monomer();

		pEndMonomer.setActinFilament(this,0);
		pEndMonomer.prev = null;
		pEndMonomer.next = bEndMonomer;
		
		bEndMonomer.setActinFilament(this,1);
		bEndMonomer.prev = pEndMonomer;
		bEndMonomer.next = null;
		
		nMonomers=2;
		finalMonomerLength = finalMons;
		
		createMyMonomers (numMons-2);
		lengthFromMonomerCount(nMonomers);
		updateAllMonomerPositions();
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
		evaluateProperties();
		
		if (isPPinned()) { pAttachment.setInitialPosition(pEnd); pAttachment.setMonomer(pEndMonomer); }
		if (isBPinned()) { bAttachment.setInitialPosition(bEnd); bAttachment.setMonomer(bEndMonomer); }
		
		capped = false;
		
		addActin(this);
		initForceTracking();
		lastCM.copy(cm);  // set initial position
		setActinColor();
		markBindingSites();
	}
	
	
	static public void setMonomerListener(ActinListener l) {
		monomerListener = l;
	}
	
	
	private void initForceTracking()  {
		if(filamentsBuckle || monitorStress) {
			if(myStressTracker == null) {
				myStressTracker = new InternalStressTracker();
			}
			myStressTracker.initForceTracking();
		}
		myoHolderConnects = new MyosinHolder[200];
		xlList = new Crosslinker[200];
		xlEndBoundList = new int [200];
		myoList = new Myosin[200];
		if(trackPinForces) {
			if(pAttachment != null)
				pAttachment.initForceTracking();
			if(bAttachment != null)
				bAttachment.initForceTracking();
		}
	}
	
	/** Create all the monomers. */
	void createMyMonomers(double NumMonomers){
		Monomer newMon;
		for (int i = 0; i < (NumMonomers); i++) {
			//Monomer will be a number representing the number of monomers away from the pointed end
			newMon = addMonomerToBarbedEnd();
			newMon.myFilIndex = i+2;	// already created pEndMonomer which is index 0
		}
	}
	
	void markBindingSites() {
		int stepsToNextBSite = actinBSiteInterval;
		int i = 0;
		Monomer curM = pEndMonomer;
		curM.isBindingSite = true;
		int totalSites = 1;
		while (curM != bEndMonomer) {
			i++;
			curM = curM.next;
			if (i == stepsToNextBSite) {
				i = 0;
				curM.isBindingSite = true;
				totalSites++;
			}
		}
		//System.out.println("Total binding sites marked = " + totalSites);
	}
		
	
	// COLLISION DETECTION ***** COLLISION DETECTION ***** COLLISION DETECTION ***** COLLISION DETECTION *****
	
	/** Register all filaments with the Collision detectors actin mesh. */
	static public void meshActin()  {
		for (int i = 0; i < actinCt; i++) {
			theActins[i].fillActinMesh();
		}
	}

	private void fillActinMesh() {
		for(int i = 0; i < segCt; i++) {
			Mesh.ACTIN_MESH.addLineSegmentToMesh(myActinNumber,segs[i][0],segs[i][1]);
		}
	}

	
	// INTEGRATION STEP ***** INTEGRATION STEP ***** INTEGRATION STEP ***** INTEGRATION STEP *****
	
	static public double getMaxFilamentMoveSinceLastCheck() {
		double val, max = 0;
		for(int i = 0; i < actinCt; i++) {
			val = theActins[i].getMaxMoveSinceLastCheck();
			if(val > max) max = val;
		}
		return max;
	}
	
	/** Looop through all filaments and call each one's step() method. */
	static public void stepAllFilaments(double dT)  {
		
		updateFreeMonomerConc();
		
		if(actinDynamicsMode.equals(FULL_DYNAMICS)) {
			nucleate();
		}

		for (int i = 0; i < actinCt; i++) {
			if (!theActins[i].removeMe) {
				theActins[i].step(dT);
			}
		}
		severFilaments();
		for (int i = 0; i < actinCt; i++) {
			if (!theActins[i].removeMe) {
				theActins[i].checkRelease();
			}
		}
	}
	
	static public void severFilaments()  {
		Monomer m;
		for(Enumeration e = severCandidates.elements(); e.hasMoreElements();) {
			m = (Monomer) e.nextElement();
			m.myFilament.severFilamentAt(m);
		}
		severCandidates.clear();
	}
	
	public void treadmill(int nMonomers) {
				
		if(nMonomers > 0) {
			pAttachment.setMonomer(null);
			bAttachment.setMonomer(null);
			for(int i = 0; i < nMonomers; i++) {
				removeMonomerFromPointedEnd();
				addMonomerToBarbedEnd();
			}
			pAttachment.setMonomer(pEndMonomer);
			bAttachment.setMonomer(bEndMonomer);
			pAttachment.setPosition(new Point2D(pAttachment.attachmentPoint.x + uVect.x*nMonomers*monLength,pAttachment.attachmentPoint.y));
			bAttachment.setPosition(new Point2D(bAttachment.attachmentPoint.x + uVect.x*nMonomers*monLength,bAttachment.attachmentPoint.y));
		}
		else {
			pAttachment.setMonomer(null);
			bAttachment.setMonomer(null);
			for(int i = 0; i < -nMonomers; i++) {
				removeMonomerFromBarbedEnd();
				addMonomerToPointedEnd();
			}
			pAttachment.setMonomer(pEndMonomer);
			bAttachment.setMonomer(bEndMonomer);
			pAttachment.setPosition(new Point2D(pAttachment.attachmentPoint.x - nMonomers*monLength,pAttachment.attachmentPoint.y));
			bAttachment.setPosition(new Point2D(bAttachment.attachmentPoint.x - nMonomers*monLength,bAttachment.attachmentPoint.y));
		}
		cmFromBarbedPointed();
		evaluateProperties();
		markBindingSites();
	}

	
	public void step(double dT) {
		//Determine whether the filament grew or shrank and from what end
		if(actinDynamicsMode.equals(FULL_DYNAMICS) && !isStatic) {
			dynamics();
		}
		
		pinForces();  // add forces to keep one or both ends fixed
		
		resistiveLoad(); // add a constant resistive load
		
		/** Force-based motions of rigid rod */
		super.step(dT);
		
		bucklingCalcs(); // calculate internal stress, sever if force greater than buckling criteria
		
		if(actinDynamicsMode.equals(FULL_DYNAMICS) && !isStatic) {
			cap();
			sever();
		}
		
		if(actinDynamicsMode.equals(RECYCLED_FILAMENTS)) {
			checkJump();
		}
		
		incSumMyos();  // for tracking myosin residency
				
		reset();
	}

	static public void updateFreeMonomerConc() {
		double dT = Sim2D.getDeltaT();
		freeMonomerConc = (numberOfMonomers - Monomer.monomerCt)/(Constants.microMolePerSquareNM*Sim2D.xDimension*Sim2D.yDimension*Sim2D.zDimension);
		double barbedEndRate = barbedEndKon*freeMonomerConc - barbedEndKoff;

		if(barbedEndRate > 0) {
			actinBOnProb = barbedEndRate*dT;
			actinBOffProb = 0;
		}
		else {
			actinBOnProb = 0;
			actinBOffProb = -barbedEndRate*dT;
		}
		double pointedEndRate = pointedEndKon*freeMonomerConc - pointedEndKoff;
		if(pointedEndRate > 0) {
			actinPOnProb = pointedEndRate*dT;
			actinPOffProb = 0;
		}
		else {
			actinPOnProb = 0;
			actinPOffProb = -pointedEndRate*dT;
		}
	}
	
	public void incSumMyos () {
		sumMyos += myoCt;
	}
	
	public double aveMyosAttached () {
		finishTime = Sim2D.simulationTime;
		double timeStepsToAve = (finishTime-startTime)/Sim2D.getDeltaT();
		double aveMyos = sumMyos/timeStepsToAve;
		sumMyos = 0;
		startTime = finishTime;
		return aveMyos;
	}
	static public double aveMyoCt () {
		double totalMyos=0;
		for (int i=0;i<Actin.actinCt;i++) {
			totalMyos+=theActins[i].myoCt;
		}
		return totalMyos/Actin.actinCt;
	}
	
	private Monomer findRareUnboundMonomer() {
		Monomer [] candidates = new Monomer[actinCt];
		int cnt = 0;
		for(Monomer m = bEndMonomer; m != null; m = m.next) {
			if(m.isFreeLinker()) {
				candidates[cnt++] = m;
			}
		}
		if(cnt == 0) return null;
		
		int fnd = (int)(Math.floor(Math.random()*cnt));
		if(fnd == cnt) fnd --;

		return candidates[fnd];
	}
						
	
	
	
	// FILAMENT DYNAMICS ***** FILAMENT DYNAMICS ***** FILAMENT DYNAMICS ***** FILAMENT DYNAMICS *****
	
	/** Check to see if we nucleate a new filament. */
	public static void nucleate(){
		if (nucleating && Math.random() < filamentNucleationProb*Sim2D.getDeltaT()) {
			makeRandomSeed();
		}
	}
	
	/** Create a new 3 monomer "seed" */
	public static void makeRandomSeed () {
		double initX = Math.random()*maxX;
		double initY = Math.random()*maxY;
		double randomAng = 2*Math.PI*Math.random();
		int monomers = 3;
		Actin a = new Actin (initX, initY, randomAng, monomers,monomers);
	}
	
	public void checkJump() {
		if(Math.random() < filamentJumpProb*Sim2D.getDeltaT()) {
			releaseAndJump();
		}
	}
	
	
	/** Check if this filament has less than 3 monomers.  If so, it dies. */
	public void checkRelease() {
		if(nMonomers<3){
			releaseAndDie();
		}
	}
	
	/** Decide whether to add a cap to the barbed end of this filament */
	public void cap () {
		if(!capping) {
			return;
		}
		if(nMonomers >= maxMon) {
			if(!capped)
				FileOps.reportln("filament " + this + " capped at maximum length");
			capped = true;
		}
		else if(Math.random() < filamentCappingProb*Sim2D.getDeltaT()) {
			capped = true;
		}
			
	}
	
	/** Test if this filament was severed.   If so, choose severing location, create the new filament and
	 * reallocate the monomers between the two filaments
	 */
	public void sever () {
	
		if (!severing) { return; }
		
		if (Math.random() > nMonomers*monLength*filamentSeveringProb*Sim2D.getDeltaT()) { return; }
		severingPlace = (int)(((nMonomers-2)* Math.random())+1);
		if(severingPlace<1) severingPlace = 1;
		if(severingPlace > nMonomers-1) severingPlace = nMonomers-1;
				
		Monomer severMon = pEndMonomer;
		for (int i = 0; i < severingPlace; i++) {
			severMon = severMon.next;
		}
		severCandidates.add(severMon);
	}
		
	public void severFilamentAt(Monomer severMon) {

		//System.out.println("sever");
		Monomer newBEndMonomer = severMon.prev;
		
		// reset positions of the endpoint monomers for both new filaments
		bEndMonomer.getLocation();
		pEndMonomer.getLocation();
		severMon.getLocation();
		newBEndMonomer.getLocation();
		
		Actin nuFil = new Actin(newBEndMonomer,this);
				
		severMon.prev = null;
		pEndMonomer = severMon;	//now the severed monomer is the new p end
		pEnd.set(pEndMonomer.getLocation());		//reset p end of filament
		bEnd.set(bEndMonomer.getLocation());
		this.reIndexMonomers();
		
		// reset length, drags, etc for old severed filament
		countMonomers();
		lengthFromMonomerCount(nMonomers);
		cm.add(pEnd,length/2,uVect);
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
		evaluateProperties();
		
		//transfer any endpoint pins
		if (isPPinned()) {
			nuFil.pAttachment = this.pAttachment;
			this.pAttachment = null;
		}
	}
	
	/** Calculate barbed end monomer addition and pointed end removal. */
	public void dynamics (){
		boolean lengthChanged = false;
		
		//  Barbed end dynamics
		if (!capped) {
			if(actinBOnProb > 0) {
				if(Math.random() < actinBOnProb){
					addMonomerToBarbedEnd();
					lengthChanged = true;
				}
			}
			else if(Math.random() < actinBOffProb) {
				if (nMonomers >= 3) {
					removeMonomerFromBarbedEnd();
					lengthChanged = true;
				}
			}
		}
		
		// Pointed end dynamics
		if(actinPOnProb > 0) {
			if(Math.random() < actinPOnProb){
				addMonomerToPointedEnd();
				lengthChanged = true;
			}
		}
		else if(Math.random() < actinPOffProb) {
			if (nMonomers >= 3) {
				removeMonomerFromPointedEnd();
				lengthChanged = true;
			}
		}
		if (lengthChanged) {
			cmFromBarbedPointed();
			evaluateProperties();
		}
	}
	
	
	/** Remove a Monomer from the pointed end. */
	public void removeMonomerFromPointedEnd(){
		
		if(monomerListener != null) {
			monomerListener.monomerReleased(pEndMonomer);
		}
		Monomer tempMonomer = pEndMonomer.next;
		pEndMonomer.next = null;					//break links, about to delete pEnd Monomer
		pEndMonomer.releaseAll();
		pEndMonomer.setForRemove();
		
		tempMonomer.prev = null;
		pEndMonomer = tempMonomer;
		
		pEnd.inc(monLength,uVect);
		
		nMonomers--;
		length -= monLength;
	
		indexOffset++;
	}
	
	/** Delete a monomer from the barbed end. */
	public void removeMonomerFromBarbedEnd(){
		
		if(monomerListener != null) {
			monomerListener.monomerReleased(bEndMonomer);
		}
		Monomer tempMonomer = bEndMonomer.prev;
		bEndMonomer.prev = null;
		bEndMonomer.releaseAll();
		bEndMonomer.setForRemove();
		
		tempMonomer.next = null;
		bEndMonomer = tempMonomer;
		
		nMonomers--;
		length -= monLength;
			
		bEnd.inc(-monLength,uVect);
	}
	
	
	/** Add a monomer to the pointed end of this filament. */
	public void addMonomerToPointedEnd(){
		//add one monomer to the pointed end
		Monomer tempMonomer = new Monomer();
		
		tempMonomer.setActinFilament(this,0); //*** BROKEN CODE
		tempMonomer.prev = null;
		tempMonomer.next = pEndMonomer;
		
		pEndMonomer.prev = tempMonomer;
		
		pEndMonomer = tempMonomer;
		
		pEnd.inc(-monLength,uVect);
		nMonomers++;
		length += monLength;
	}
	
	/** Add a monomer to the barbed end of this filament. */
	public Monomer addMonomerToBarbedEnd(){

		//add one monomer to the barbed end
		Monomer newBEndMon = new Monomer();
		Monomer oldBEndMon = bEndMonomer;
		
		newBEndMon.setActinFilament(this,oldBEndMon.myFilIndex+1);
		newBEndMon.prev = oldBEndMon;
		newBEndMon.next = null;
		
		oldBEndMon.next = newBEndMon;
		bEndMonomer = newBEndMon;
		bEnd.inc(monLength,uVect);
		
		nMonomers++;
		length += monLength;
		
		return bEndMonomer;
	}
	
	/** Melt this entire filament, releasing all of it's monomers binding partners and scheduling
	 * all monomers for destruction.
	 */
	public void meltFilament(){

		nMonomers=0;
		pEndMonomer.next=null;
		pEndMonomer.releaseAll();
		pEndMonomer.setForRemove();
		
		bEndMonomer.prev=null;
		bEndMonomer.releaseAll();
		bEndMonomer.setForRemove();
		
		pEndMonomer=null;
		bEndMonomer=null;
		setForRemove();
	}
		
	//. FORCES **** FORCES **** FORCES **** FORCES **** FORCES **** FORCES **** FORCES **** FORCES
	
	
	/** Add the force Ñfrom a myosinÑ at the specified location, incrementing both the net force on the center of mass, and the net torque.*/
	public void addMyoForce (Point2D f, Point2D loc, Myosin myo){
		addForce(f, loc);
		registerMyo(myo);
		
		/* remove this.. just a check for gliding analysis */
		double fmag = Point2D.dot(f, uVect);
		if (fmag > 0) { posMyoForce += fmag; } else { negMyoForce += fmag; }
	}
	
	public void addXLForce (Point2D f, Point2D loc, Crosslinker xl, int endBound) {
		addForce(f, loc);
		registerXL(xl, endBound);
	}
	
	/** Register crosslinkers that are adding force to this filament */
	public void registerXL(Crosslinker xl, int endBound) {
		// reset counters if this is the first entry for this time-step
		if (xlTime != Sim2D.simulationTime) {
			xlTime = Sim2D.simulationTime;
			xlCt = 0;
		}
		xlList[xlCt] = xl;
		xlEndBoundList[xlCt] = endBound;
		xlCt++;
	}
	
	/** Register myosins that are adding force to this filament */
	public void registerMyo(Myosin myo) {
		// reset counters if this is the first entry for this time-step
		if (myoTime != Sim2D.simulationTime) {
			myoTime = Sim2D.simulationTime;
			myoCt = 0;
		}
		if(myoCt >= myoList.length) {
			Myosin [] tmp = new Myosin[myoList.length+100];
			System.arraycopy(myoList,0,tmp,0,myoList.length);
			myoList = tmp;
		}
		myoList[myoCt] = myo;
		myoCt++;
	}
	
	/** Add MyosinMiniFilament to connected list, if not already present */
	public void addMyoConnect (MyosinHolder myoHolder) {
		// reset counters if this is the first entry for this time-step
		if (connectTime != Sim2D.simulationTime) {
			connectTime = Sim2D.simulationTime;
			myoCt = 0;
			myoHolderConnectCt = 0;
		}
		
		boolean alreadyConnected = false;
		for (int i=0; i<myoHolderConnectCt; i++) {
			if (myoHolder == myoHolderConnects[i]) {
				alreadyConnected = true;
				break;
			}
		}
		if (!alreadyConnected) {
			myoHolderConnects[myoHolderConnectCt] = myoHolder;
			myoHolderConnectCt++;
		}
	}
	
	public void makeMyoFilConnectList () {
		for (int i=0;i<myoCt;i++) {
			if (myoList[i].myHolder != null) {
				addMyoConnect (myoList[i].myHolder);
			}
		}
	}
	
	public static void traceAllClusters () {
		Cluster.resetClusters();
		Actin curA;
		for (int i=0;i<Actin.actinCt;i++) {
			curA = Actin.theActins[i];
			if (curA.clusterTime != Sim2D.simulationTime) {
				Cluster.startCluster(curA);
			}
		}
	}
	
	public void setConnected (Cluster myCluster) {
		if (clusterTime == Sim2D.simulationTime) return;	// this actin already a member
		myCluster.addElement(this,true);
				
		if (Cluster.followMyos && myoTime == Sim2D.simulationTime-Sim2D.getDeltaT()) {

			makeMyoFilConnectList();  // compact list of myosins into unique list of myosin holders
			for (int i=0;i<myoHolderConnectCt;i++) {
				myoHolderConnects[i].setConnected(myCluster);
			}
		}
		
		if (Cluster.followXLs && xlTime == Sim2D.simulationTime-Sim2D.getDeltaT()) {

			for (int i=0;i<xlCt;i++) {
				xlList[i].setConnected(myCluster, xlEndBoundList[i]);
			}
		}
	}
	
	public void setConnected (int conLevel) {
		connectedLevel = Math.min(conLevel, connectedLevel); //always select most closely connected path as connectedLevel
		if (isConnected) return;
		isConnected = true;
		
		if (myoTime == Sim2D.simulationTime) {
			makeMyoFilConnectList();  // compact list of myosins into unique list of myosin holders
			for (int i=0;i<myoHolderConnectCt;i++) {
				myoHolderConnects[i].setConnected(connectedLevel+1);
			}
		}
		
		if (xlTime == Sim2D.simulationTime) {
			for (int i=0;i<xlCt;i++) {
				xlList[i].setConnected(connectedLevel, xlEndBoundList[i]);
			}
		}
	}
	
	public void traceTendril (ForceTendril myTendril) { }
	
	/** Resistive load... simulates loads to reproduce Howard Fig 16.5 */
	public void resistiveLoad () {
		if(applyConstantFractionStallForce){
			int nHeadsEquil = MyosinMiniFilament.predictNumHeadsAtStall();
			double equilStallF = nHeadsEquil*Myosin.getCrossbridgeForce();
			resistiveLoad=fractionStallForceApplied*equilStallF;
		}
		aLoad.copy(uVect);
		aLoad.scale(resistiveLoad);
		addForce(aLoad);
	}
	
	public void relaxAllAttachments() {
		if(pAttachment!=null) {
			pAttachment.relax();
		}
		if(bAttachment!=null) {
			bAttachment.relax();
		}
	}
	
	public void restoreAttachmentPositions() {
		pAttachment.restoreInitialPosition();
		bAttachment.restoreInitialPosition();
		cm.set(0.5*(pAttachment.attachmentPoint.x+bAttachment.attachmentPoint.x),
									0.5*(pAttachment.attachmentPoint.y+bAttachment.attachmentPoint.y));
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
	}
				
	public boolean isPPinned() {
		return pAttachment!= null;
	}
	
	public boolean isBPinned() {
		return bAttachment!= null;
	}
	
	public FilamentAttachment  getPEndAttachment() {
		return pAttachment;
	}
	
	public FilamentAttachment  getBEndAttachment() {
		return bAttachment;
	}
	
	public void stretchPEndAttachment(double f) {
		double d = f/pAttachment.getElasticResistance();
		pAttachment.setPosition(Point2D.sum(pEnd,Point2D.scaledV(d,uVect)));
	}
	
	public void stretchBEndAttachment(double f) {
		double d = f/bAttachment.getElasticResistance();
		bAttachment.setPosition(Point2D.sum(bEnd,Point2D.scaledV(d,uVect)));
	}
	
	/* Pin forces, if one or more ends are pinned  */
	public void pinForces () {
		
		if (isPPinned()) {
			pAttachment.getElasticForce();
		}
		
		if (isBPinned()) {
			bAttachment.getElasticForce();
		}
		
	}
	
	// Determination of buckling of internal segments
	public void bucklingCalcs () {
		if(monitorStress || filamentsBuckle) {
			myStressTracker.trackStresses(this);
		}
		if(filamentsBuckle) {
			myStressTracker.severFromBuckling(this);
			}
		}
		
	public void addSeverCandidate(Monomer m) {
		severCandidates.add(m);
	}

	public void lengthFromMonomerCount(int nMonomers) {
		length = monLength*(nMonomers-1.0);
		setPhysicalLength(monLength*nMonomers);
	}
	
	/** Determine the positions of pointed and barbed ends from center of mass and orientation (µVect) */
	public void useCenterMassAndUVectToDetermineEnds(){
		end1.add(cm,-length/2,uVect);
		end2.add(cm,length/2,uVect);
	}
	
	/** update monomer count and then let RigidRod do the all the physical property calcs. */
	public void evaluateProperties () {
		lengthFromMonomerCount(nMonomers);
		super.evaluateProperties();
	}
	
	/** A calculation of the ratios of baseline (i.e. no force-dependent kinetics) myosin crossbridge
	 * lifetime to the relaxtion time of crossbridge-viscously moving filament relaxation time
	 */
	public static double getNoNameRatio() {
		// If this ratio is << 1 then the motion of filaments is determined by motor dynamics
		// not by the time required to move the filament through a viscous media.
		double crossBridgeLifetime = Myosin.getCrossBridgeLifetime();
		double crossBridgeRelax = theActins[0].actinfilTransGamma/Myosin.myoStrokeSpringConstant;
		System.out.println("crossBridgeRelax = " + crossBridgeRelax + " gamma = " + theActins[0].actinfilTransGamma);
		return crossBridgeRelax/crossBridgeLifetime;
	}

	static public double getTotalFilamentLength() {
		double l = 0;
		for(int i = 0; i < actinCt; i++) {
			l+= theActins[i].length;
		}
		return l;
	}
	
	static public double [] getWeightedMonomerChoices() {
		double l = 0, totalL = getTotalFilamentLength();
		double [] dist = new double[actinCt];
		for(int i = 0; i < actinCt; i++) {
			l+= theActins[i].length;
			dist[i] = l/totalL;
		}
		return dist;
	}
			
	
	/** Bind crosslinker from the "cytoplasmic" pool, A new Crosslinker object is created de novo as it binds. */
	public Monomer bindCrossLinker() {
		int cnt = 0;
		Monomer m = null;
		while(m == null) {
			double arcL = Math.random()*length;
			m = getMonomerAt(arcL);
			if (!m.isFreeLinker()) {
				m = null;
				if(++cnt > nMonomers) {
					m = findRareUnboundMonomer();
					return m;
				}
			}
		}
		if(m.myFilIndex < indexOffset) {
			FileOps.reportln("crosslinker born with attachement to bad monomer index");
		}
		return m;
	}


	/** Set the new center of mass from the positions of barbed and pointed ends. */
	public void cmFromBarbedPointed(){
		pTobVec.getVector(pEnd, cm);
		pTobVec.uVect();
		cm.add(pEnd,length/2,pTobVec);
		cm.wrapPoint();
	}

		
	/** Set the position of Monomer relative to this filament's position and orientation. */
	public void setMonomerPosition(Monomer m){
		m.cm.add(pEnd,(m.myFilIndex-indexOffset)*monLength,uVect);
		m.cm.wrapPoint();
	}
	
	/** Set the position of the Monomer relative to this fliament's position and orientation
	 * asumming that it's position is given by index.
	 */
	public void setMonomerPosition(Monomer m, int index){
		m.cm.add(pEnd,index*monLength,uVect);
		m.cm.wrapPoint();
	}
	
	// Return the monomer at the specified arclength along this filament starting from pointed end. */
	public Monomer getMonomerAt (double arcL) {
		int monIndex=(int)(Math.round(arcL/monLength));
		int ct=0;
		Monomer tempMonomer=pEndMonomer;
		while(ct<monIndex&&tempMonomer!=null){
			tempMonomer=tempMonomer.next;
			ct++;
		}
		if(tempMonomer != null) {
			setMonomerPosition(tempMonomer,ct);
		}
		return tempMonomer;
	}

	/** Returns the current numbr of monomers in this filament. */
	public void countMonomers(){
		nMonomers=0;
		Monomer tempMonomer = pEndMonomer;
		while (tempMonomer != null){
			nMonomers++;
			tempMonomer = tempMonomer.next;
		}
	}
	
	/** Print info about monomers to the console. */
	public void printMonomerInfo(){
		int i=0;
		Monomer tempMonomer = pEndMonomer;
		while (tempMonomer != null){
			System.out.println("i="+Double.toString(i)+" mon="+tempMonomer+"  on "+this+"  with "+nMonomers+" totalMons");
			
			tempMonomer = tempMonomer.next;
			i++;
		}
	}
	
	/** Updates the positions of all monomers relative to this filament's position and orientation */
	public void updateAllMonomerPositions(){
		double i=0;
		Monomer tempMonomer = pEndMonomer;
		while (tempMonomer != null){
			if (!tempMonomer.isFree()) {
				tempMonomer.cm.add(pEnd,i*monLength,uVect);
				tempMonomer.cm.wrapPoint();
			}
			tempMonomer = tempMonomer.next;
			i++;
		}
	}
	
	/** Reassign monomer indices */
	public void reIndexMonomers() {
		int newIndex = 0;
		Monomer tM=pEndMonomer;
		while(tM !=null){
			tM.myFilIndex = newIndex;
			newIndex++;
			tM=tM.next;
		}
	}
	

	/** Cycle through all monomers and call each one's releaseAll() method to relaease all binding partners. */
	public void releaseAll () {
		releaseMonomers();
		zeroForces();
		indexOffset = 0;
	}
	
	/** Release all bound elements, but leave filament and associated monomers in place. */
	public void releaseBoundElements() {
		Monomer tM=pEndMonomer;
		while(tM !=null){
			tM.releaseAll();
			tM=tM.next;
		}
	}
	
	/** Diissolve the filament, releasing all monomers and bound elements and scheduling monomers for removal. */
	public void releaseMonomers() {
		Monomer tM=pEndMonomer;
		while(tM !=null){
			tM.releaseAll();
			tM.prev = null;
			tM.setForRemove();
			tM=tM.next;
		}
		pEndMonomer = null;
		bEndMonomer = null;
	}
	
	
	/** Release all binding partners and then jump to a new random location. */
	public void releaseAndJump () {
		releaseBoundElements();
		reIndexMonomers();
		jump();
	}
	
	/** Dissolve filament, releasing all binding partners and then die. */
	public void releaseAndDie () {
		releaseAll();
		setForRemove();
	}
	
	

	
// DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING ***** DRAWING *****
	
	double angle=0;
	GeneralPath perimeter = new GeneralPath(); // create GeneralPath object
	boolean showBox=false;
	boolean showEnds=true;
	boolean showEndStrings=false;
	boolean showPinForces=false;
	boolean showCenterMass=false;
	boolean showCollisionStuff=false;
	Point2D tempPointP = new Point2D();
	Point2D tempPointB = new Point2D();
	Point2D tempPointp1I = new Point2D();
	Point2D tempPointp2I = new Point2D();
	
	
	public void setActinColor() {
		actinColor = Color.red;
	}
	
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		
		g.setColor(actinColor);
		super.drawYourself(g,scale,offset);
		
		int XpEnd = pEnd.getPixX();
		int YpEnd = pEnd.getPixY();
		int XbEnd = bEnd.getPixX();
		int YbEnd = bEnd.getPixY();
		
		// BELOW IS CRAP FOR TESTING THE COLLISON DETECTOR
		if (showCollisionStuff) {
			g.setColor(Color.yellow);
			g.drawOval(tempPointP.getPixX(),tempPointP.getPixY(),5,10);
			g.drawOval(tempPointB.getPixX(),tempPointB.getPixY(),5,10);
			g.drawOval(tempPointp1I.getPixX(),tempPointp1I.getPixY(),5,10);
			g.drawOval(tempPointp2I.getPixX(),tempPointp2I.getPixY(),5,10);
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
			int EndPixelRadius = (int) (4*actinRadius*scale);
			int EndPixelDiameter = (int) (8*actinRadius*scale);
			
			g.setColor(Color.YELLOW);
			if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) g.setColor(fadedColor);
			g.fillOval(XpEnd-EndPixelRadius, YpEnd-EndPixelRadius, EndPixelDiameter, EndPixelDiameter);
			
			if(capped) {
				g.setColor(Color.GREEN);
				if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) g.setColor(fadedColor);
				EndPixelRadius = (int) (5*actinRadius*scale);
				EndPixelDiameter = (int) (10*actinRadius*scale);
			} else {
				g.setColor(Color.CYAN);
				if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) g.setColor(fadedColor);
			}
			
			g.fillOval(XbEnd-EndPixelRadius, YbEnd-EndPixelRadius, EndPixelDiameter, EndPixelDiameter);
		}
		
		if(showEndStrings){
			int padding=20;
			g.setColor(Color.WHITE);
			g.drawString("P_END",XpEnd, YpEnd+padding);
			
			g.setColor(Color.CYAN);
			g.drawString("B_END",XbEnd, YbEnd+padding);
		}
		
		if (showPinForces) {
			int padding=20;
			g.setColor(Color.WHITE);
			
			if (isPPinned()) {
				//g.drawString("F="+pinForceFormat.format(Point.getMag(pPinForceTrack.averagePtVal())),XpEnd+padding, YpEnd+padding);
				g.drawString("Fx="+pinForceFormat.format(pAttachment.getTimeAveragedForce().x),XpEnd, YpEnd+forceTextLocator*padding+5);

			}
			
			if (isBPinned()) {
				//g.drawString("F="+pinForceFormat.format(Point.getMag(bPinForceTrack.averagePtVal())),XbEnd+padding, YbEnd+padding);
				g.drawString("Fx="+pinForceFormat.format(bAttachment.getTimeAveragedForce().x),XbEnd, YbEnd+forceTextLocator*padding+5);
				//System.out.println("Ave bPinForce = " + Point.getMag(pinForceTrack.averagePtVal()));
			}
		}
		
		if(showCenterMass){
			g.setColor(Color.MAGENTA);
			xCenter = cm.getPixX();
			yCenter = cm.getPixY();
			pixelRadius = (int) (actinRadius*scale);
			pixelDiameter = (int) (2*actinRadius*scale);
			
			g.fillOval(xCenter-pixelRadius,yCenter-pixelRadius,pixelDiameter,pixelDiameter);
		}
		
		if (showInternalStress  && myStressTracker!=null) {
			myStressTracker.showInternalStress(this,g,scale);
		}
	}
	
	public void setActinLineColor (Graphics g) {
		g.setColor(actinColor);
		//if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) g.setColor(fadedColor);
		if (Sim2D.paintConnected) {
			g.setColor(fadedColor);
			for (int i=0;i<Cluster.bigClusterCt;i++) {
				try {
					if (Cluster.sortedClusters[i].clusterId == clusterId) { g.setColor(Sim2D.clusterColors[i]); }
				} catch (NullPointerException npe) { }
			}
		}
	}
	
	
	public void perimeterRender (Point2D end1, Point2D end2, Graphics g, double scale, double [] offset) {
		Graphics2D g2d=(Graphics2D)g;
		double height=actinRadius;
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
		
	public static void addActin (Actin newActin) {
		newActin.myActinNumber = actinCt;
		theActins[actinCt] = newActin;
		actinCt ++;
	}

	public void die() {
		pEndMonomer = null;
		bEndMonomer = null;
	}
	
		
	
	public void removeMe () {
		super.removeMe();
		int swapId = myActinNumber;
		theActins[swapId] = theActins[actinCt-1];
		theActins[swapId].myActinNumber = swapId;
		actinCt --;
	}
	
	public void printMe() {
		FileOps.reportln(this+"myActinNUmber="+myActinNumber+"  thingNumber="+myThingNumber+"  simTime="+Sim2D.simulationTime+" nMons="+nMonomers+"  pEndMon="+pEndMonomer+"  bEndMon="+bEndMonomer+"  removeMe="+removeMe);
		printMonomerInfo();
	}
	
	public static void printAll () {
		System.out.println("ALL ACTIN PRINT");
		for (int i=0;i<Actin.actinCt;i++) {
			Actin.theActins[i].printMe();
		}
		FileOps.reportln("END ALL ACTIN PRINT");
	}
		
	/**Filaments remain at a fixed length; no treadmilling, no nucleation, no severing. */
	static public String STATIC_FILAMENTS = new String("StaticFilaments");
	
	/** we make fixed number of actins at the outset, whose length before capping
	 * is predetermined.  These grow from barbed ends till they cap, then shrink from
	 * pointed ends to <3 monomers, then they uncap, and jump to a new random position and
	 * orientation and start over.
	 */
	static public String RECYCLED_FILAMENTS = new String("RecycledFilaments");
	
	/** filaments grow and shrink from barbed and pointed ends, exchanging
	 * subunits with the bulk monomer pool.  Nucleation, capping and severing can be turned on
	 * and off separately via the boolean flags Constants.nucleating, Constants.severing,
	 * Constants.capping.  The rates of these reactins are controlled by tunable parameters (below)
	 */
	static public String FULL_DYNAMICS = new String("FullDynamics");

	
	/*********************************************************/
	/**********************  ACTIN PARAMETERS  ***************/
	/*********************************************************/
	/**apply external force as fraction of stall force*/
	static public double fractionStallForceApplied;
	
	/**pulse external force for interval of time*/
	static public double pulseTime;
	
	static public String actinDynamicsMode = Actin.FULL_DYNAMICS;
	
	static public boolean severing;
	
	static public boolean nucleating;
	
	static public boolean capping;
	
	static public boolean useUniformInitialLengthDistribution;
	
	/** Initial number of actin filaments in this simulation. */
	static public int numInitialActins;

	/** Monomer step between myosin binding sites */
	static public int actinBSiteInterval;
	
	/** Total concentration of the actin monomer pool in µM, which is depleted by filament assembly */
	static public double monomerConc;
	
	/** Probability of nucleating a new filament this time step. */
	static public double filamentNucleationProb;

	
	/** Probability of this filament jumping (recycling)this time step. */
	static public double filamentJumpProb;

	/** Probability that a filament will sever per micron per second */
	static public double filamentSeveringProb;
	
	/** Filament capping probability */
	static public double filamentCappingProb;
	
	/** avg length of actin filaments. */
	static public double avgActinLength;
	
	/** standard deviation, in nm, of gaussian actin length distribution,
	 * normalized by average length . */
	static public double stdDevActinLength;
	
	/** Target Actin density in monomers (nm^2) when making actin to a prespecified density */
	static public double targetActinDensity;
	
	/** pointed end polymerization rate (s-µM) */
	static double pointedEndKon;
	
	/** pointed end depolymerization rate (0.1 - 0.25/s) */
	static double pointedEndKoff;

	/** barbed  end polymerization rate (5-50/s-µM) */
	static double barbedEndKon;

	/** pointed end depolymerization rate (/s) */
	static double barbedEndKoff;
	
	/** Second moment of inertia for actin */
	static double actinEI;  // EI in pN-nm  (6e-26 N-m^2)(1e12 pN/N)(1e18 nm^2/m^2)
	
	static double buckleFactor;
	
	/** load on filament toward plus-end... for resistive load benchmarking */
	public static double resistiveLoad;
	
	public static boolean applyConstantFractionStallForce;
	
	static boolean filamentsBuckle;
	
	static boolean monitorStress;

	static boolean showInternalStress;
	
	static boolean trackPinForces;
	
	static String className = new String("main.Actin");
		
	/*********************************************************/
	/**********************  PARAMETERS  *********************/
	/*********************************************************/

	static public boolean setActinDynamicsMode(String s)  {
		if(s.equals(Actin.STATIC_FILAMENTS) || s.equals(Actin.RECYCLED_FILAMENTS) || s.equals(Actin.FULL_DYNAMICS)) {
			actinDynamicsMode = s;
			return true;
		}
		else {
			return false;
		}
	}
	
	static void staticInit() {
		
		Parameters.addParameterSetListener (
			new ParameterSetListener() {
				public void parametersChanged() {
					maxMon = (int) (Sim2D.maxUniqueXDist/monLength)-2;
					numberOfMonomers = monomerConc*Sim2D.xDimension*Sim2D.yDimension*Sim2D.zDimension*Constants.microMolePerSquareNM;
					buckleFactor = (Math.pow(Math.PI,2)*actinEI/4);
				}
			}
		);
			
		Parameters.addParameter(className,"actinDynamicsMode",Actin.FULL_DYNAMICS,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					setActinDynamicsMode(p.getStringValue());
				}
			}
		);
		Parameters.addParameter(className,"fractionStallForceApplied",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						fractionStallForceApplied = p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter(className,"applyConstantFractionStallForce",true,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						applyConstantFractionStallForce = p.getBooleanValue();
					}
				}
			);
		Parameters.addParameter(className,"pulseTime",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						pulseTime = p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter(className,"trackPinForces",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					trackPinForces = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"severing",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					severing = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"filamentsBuckle",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					filamentsBuckle = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"monitorStress",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					monitorStress = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"showInternalStress",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					showInternalStress = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"nucleating",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					nucleating = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"capping",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					capping = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"useUniformInitialLengthDistribution",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					useUniformInitialLengthDistribution = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className,"numInitialActins",2,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					numInitialActins = (int)p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"actinBSiteInterval",1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					actinBSiteInterval = (int)(p.getDoubleValue());
				}
			}
		);
		Parameters.addParameter(className,"monomerConc",0.5,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					monomerConc = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"filamentJumpProb",1.0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					filamentJumpProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"filamentNucleationProb",10,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					filamentNucleationProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"filamentSeveringProb",0.0003,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					filamentSeveringProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"filamentCappingProb",0.01,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					filamentCappingProb = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"avgActinLength",1000,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					avgActinLength = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"stdDevActinLength",1.0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					stdDevActinLength = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"targetActinDensity",0.01,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					targetActinDensity = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"pointedEndKon",0.0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					pointedEndKon = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"pointedEndKoff",0.25,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					pointedEndKoff = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"barbedEndKon",50,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					barbedEndKon = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"barbedEndKoff",0,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					barbedEndKoff = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"actinEI",6.0e-26*1e30,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					actinEI = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className,"resistiveLoad",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						resistiveLoad = p.getDoubleValue();
					}
				}
			);
		
		Sim2D.addProxy(new AgentProxy() {
			public void registerForCollisionDetection() {
				meshActin();
			}
			public void step(double dT) {
				// Actin filaments move and update positions of their monomers.
				stepAllFilaments(dT);
			}
		
			public void reset() {
				for(int i = 0; i < actinCt; i++) {
					theActins[i].die();
					theActins[i] = null;
				}
				actinCt = 0;
			}
		});
		
	}
	
	
}



