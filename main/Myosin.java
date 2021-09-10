package main;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.*;

import util.*;
import collision.*;
import parameters.*;
import gui.*;
import analysis.*;

public class Myosin extends Thing {

	static public String
	GLIDE = new String("Glide"),
	PINNED = new String("Pinned");

	static CollisionDetector myCollisionDetector = new MyoCollisionDetector();

	/** An array containing all the Myosins */
	static public Myosin [] theMyosins = new Myosin [1000000];

	/** Current number of Myosins. */
	static public int myosinCt = 0;

	/** Diameter of a myosin head in nm */
	static double headDiameter;

	/** How close an actin filament must be to the stalk's base for binding to occur. */
	static public double minBindingDistance;

	/** standard orthogonal rest distance for myosin binding */
	static double myoStdOrthDist;

	static Vector cycles = new Vector();

	/** Track duty ratio */
	double timeOn = 0;
	double timeOff = 0;
	double timeOff_tmp = 0;
	double timeOn_tmp = 0;
	double lastTimeOn = 0;
	double lastTimeOff = 0;

	public boolean closeEnoughToBind = false;
	static DecimalFormat dutyRatioOutFormat = new DecimalFormat ("#0.000#; #0.000#");	// a clean format for duty ratio printing


	/**Track number of possible actin binding sites*/
	public int numPossibles;

	/** Track myosin force aves */
	static String hackPath = "/Users/jba/Desktop/";
	static HistogramPlus aveMyoForce = new HistogramPlus (100,0,0.1,hackPath,"aveMonForce",null,false,true,false,false);

	/** biochemical states */
	static final byte
	NONE = 0,
	ATP = 1,
	ADPPi = 2,
	ADP = 3;
	byte atpState = ADPPi;

	/** bound to filament states */
	static final byte
	FREE = 0,
	BOUND = 1,
	BOUND_PS = 2;

	byte boundState = FREE;

	/** Whether this myosin head is at the "A" end or the "B" end of the filament. */
	double orientation;

	/** Position of this myosin in the array of all myosins. */
	int myMyosinNumber;

	/** The MyosinHolder this myosin is attached to. */
	public MyosinHolder myHolder;

	/** True if this myosin head is at "A" end of miniFilament. */
	public boolean iAmAnAMyo;

	/** The actin monomer this myosin head is bound to. */
	public Monomer boundMon;

	/** The point at which the myosin stalk attaches to the myosin holder expressed
	 * with respect to the holder's local coordinate system. */
	Point2D attachPtLocal=new Point2D();

	/** The point at which the myosin stalk attaches in the myosin holder in global coordinates. */
	public Point2D attachPt=new Point2D();

	/** The location of the myosin head.  This is the actin monomer location if bound, the attachPt if free */
	public Point2D bindingSite = new Point2D();

	/** Orthogonal drop to filament bound by this myosin*/
	Point2D orthoPt = new Point2D();
	double orthoOffset = 0; // initial distance from attachPt to orthoPt
	double tangOffset = 0;  // initial signed distance from orthoPt to bindingSite

	Point2D forceSite = new Point2D();
	Point2D headJoint = new Point2D();  // an alternative location for force application
	Point2D headJointUVec = new Point2D(); // actual vector from attachPt to headJoint
	Point2D headJointTargetUVec = new Point2D(); // vector which torque tries to align actual headJointUVec to

	/** Vector reused for rendering across periodic boundaries */
	Point2D renderVec = new Point2D();

	static Color freeColor = new Color (80,80,0);

	/** A unit vector pointing from the miniFilament center to the attachment site for this myosin. */
	Point2D oriUVect = new Point2D();

	/** Unit vector to forceSite.  Used only for angle-of-binding constraints*/
	Point2D bindingUVec = new Point2D();
	/** Copy of unit vector of filament bound to at present */
	Point2D uVecFil = new Point2D();
	/** Force this myosin head exerts along the actin filament. Reused in force calculations. */
	Point2D uVecTang = new Point2D();
	Point2D forceVecTang = new Point2D();
	double forceMagTang = 0;
	double forceMagTangSigned = 0;	// signed +/- if toward +/- end, respectively
	/** Force this myosin exerts orthogonal to the actin filament */
	Point2D forceVecOrtho = new Point2D();
	double forceMagOrtho = 0;

	/** Time-averaged force on this myosin head - used for force-dependent binding kinetics. */
	double forceAv;
	static double maxForce = 0;  // stores max tangential force for all myos... for testing

	/** How far the myosin-actin bind is stretched defined as the distance between the
	 * myosin/miniFilament attachment site and the myosin/Factin binding site minus the stalk length.
	 */
	double bondStretch = 0;

	static Myosin threeBeadMyo;

	/** Color to draw this myosin. */
	Color myosinColor = Color.yellow;

	/** flag used in rendering from files */
	boolean fileMyoOffCortex = false;

	public double usedReleaseRateMult=0;

	public Myosin (MyosinHolder myf) {
		super(0,0);
		myHolder = myf;
		addMyosin(this);
	}



	static public void stepAllMyosins(double dT) {
		for (int i = 0; i < myosinCt; i++) {
			theMyosins[i].step(dT);
		}
	}


	public void setInitialPositionRelativeToMiniFilament(Point2D relMiniFilPosition){
		this.attachPtLocal.set(relMiniFilPosition);
	}

	public Point2D getLocation() {
		return bindingSite;
	}

	public void setBoundState(byte st) {
		boundState = st;
	}

	public void miniStep() {
		updateSites();
		doForces();
	}


	public void step (double dT) {
		//reportRates();
		updateBindingTimes(dT);
		updateSites();
		doForces();
		biochemStep(dT);
		//headCooperativity();
		reset();

	}

	public void biochemStep (double dT){
		switch (atpState) {
		case ATP:
			if (!isFree() && (Math.random() < myoUniformReleaseProb*dT)) { releaseMon(); }
			hydrolize(dT);
			break;
		case ADPPi:
			if (boundState != FREE) { setBoundState(BOUND); }
			dissociatePi(dT);
			break;
		case ADP:
			if (boundState != FREE) { setBoundState(BOUND_PS); }
			dissociateADP(dT);
			break;
		case NONE:
			if (boundState != FREE) { setBoundState(BOUND_PS); }
			atpOnMyo(dT);
			break;
		}
	}

	public void reportRates () {
		System.out.println("atpOnMyo rate is " + atpOnMyoRate*Sim2D.getDeltaT());
		System.out.println("hydrolize (off-fil) rate is " + myoOffFilATP_ADPPiRate*Sim2D.getDeltaT());
		System.out.println("dissociate (on-fil) rate is " + myoOnFilADPPi_ADPRate*Sim2D.getDeltaT());
		System.out.println("adp (on-fil) release rate multiplier is " + adpReleaseRateMult);
	}

	public void atpOnMyo (double dT) {
		if (Math.random() < atpOnMyoRate*dT) {
			setStateATP();
		}
	}

	public void hydrolize (double dT){
		if (boundState != FREE) {
			double hydroProbOnFil = myoOnFilATP_ADPPiRate*dT;
			if (Math.random() < hydroProbOnFil) { setStateADPPi(); }
		} else {
			if (Math.random() < myoOffFilATP_ADPPiRate*dT) { setStateADPPi(); }
			bindingSite.zero();
		}
	}

	public void dissociatePi(double dT) {
		if (boundState != FREE) {
			double disPiOnFil = myoOnFilADPPi_ADPRate*dT;
			if (Math.random() < disPiOnFil) { setStateADP(); }
		} else {
			if (Math.random() < myoOffFilADPPi_ADPRate*dT) { setStateADP(); }
			bindingSite.zero();
		}
	}
	double numADPReleaseEvents=0;
	public void dissociateADP(double dT) {
		double releaseRate = getADPReleaseRate(forceMagTangSigned);
		if (Math.random() < usedReleaseRateMult*releaseRate*dT) {
			setStateNONE();
			lastTimeOff = Sim2D.simulationTime;
			timeOn_tmp = lastTimeOff - lastTimeOn;
			cycles.add(new double [] {timeOn_tmp,timeOff_tmp});
			numADPReleaseEvents++;
		}

	}
	public double getNumADPReleaseEvents(){
		return numADPReleaseEvents;
	}
	public void zeroNumADPReleaseEvents(){
		numADPReleaseEvents=0;
	}

	/** Dissociate ADP with a force-based rate of Ae-kF */
	static public double getADPForceBasedReleaseRate(double F) {

		// the following is based on the curve fit from Veigel et al. 2003 for smooth-muscle myosin
		double k =  1.380662e-2; // Boltzmann's constant in nm-pN-s^-2-K^-1
		double T = 298.15; // temp in K
		double d = 2.7; // distance parameter from Veigel 2003 fit to exp data
		double fudge = 0.6;
		double k1 = myoOnFilADP_NoneRate*Math.exp(-F*fudge*d/(k*T));
		//System.out.println("min kOff = " + myoOnFilADP_NoneRate_local*Math.exp(-Myosin.myoMaxSingleForce*d/(k*T)));
		if (k1 > 2000) k1=2000;
		return k1;
	}

	/** Dissociate ADP with a force-based rate: Catch-Slip from Guo & Guilford 2006*/
	static public double getADPCatchSlipReleaseRate(double F) {

		// the following is based on the curve fit from Guo & Guilford 2006 for smooth-muscle myosin
		double kOffSlip = 15;
		double xSlip = 0.40;
		double kOffCatch = 176;
		double xCatch = -2.5;
		double k =  1.380662e-2; // Boltzmann's constant in nm-pN-s^-2-K^-1
		double T = 298.15; // temp in K

		double bellSlip = kOffSlip*Math.exp(F*xSlip/(k*T));
		double bellCatch = kOffCatch*Math.exp(F*xCatch/(k*T));

		return (bellSlip + bellCatch);
	}

	static public double getMaxForce() {
		return myoStepSize*myoStrokeSpringConstant;
	}

	static public double getADPReleaseRate(double F) {
		if(adpReleaseMode.equals(DISSOCIATE_ADP_CATCH_SLIP))
			return getADPCatchSlipReleaseRate(F);
		else if(adpReleaseMode.equals(DISSOCIATE_ADP_FORCE_BASED))
			return getADPForceBasedReleaseRate(F);
		else
			return myoOnFilADP_NoneRate;
	}


	static public double predictDutyRatio(double F)  {
		double releaseRate = adpReleaseRateMult*getADPReleaseRate(F);
		double t_all = 1./atpOnMyoRate + 1./myoOffFilATP_ADPPiRate + 1./myoUniformBindingProb + 1./myoOnFilADPPi_ADPRate + 1./releaseRate;
		return (1./myoOnFilADPPi_ADPRate + 1./releaseRate)/t_all;
	}

	static public double [] getAverageDutyRatio() {
		if(cycles.size() < 20) return null;
		int cnt = 0;
		double [] cycle;
		double [] dr = new double[3];
		for(Enumeration e = cycles.elements(); e.hasMoreElements(); ) {
			cycle = (double []) e.nextElement();
			dr[0] += cycle[0];
			dr[1] += cycle[1];
			dr[2] += (cycle[0]/(cycle[0]+cycle[1]));
			cnt++;
		}
		cycles = new Vector();
		dr[0]/=cnt;dr[1]/=cnt;dr[2]/=cnt;
		return dr;
	}


	public double getDutyRatio () {
		double totalTime = timeOn+timeOff;
		if (totalTime == 0) {
			return 0;
		} else {
			return timeOn/totalTime;
		}
	}

	public static String getDutyRatioString() {
		return " dutyRatio = " + String.valueOf(dutyRatioOutFormat.format(getAveDutyRatio()));
	}

	public static double getAveDutyRatio () {
		double aveDutyRatio = 0;
		int myosConsidered = 0;
		for (int i=0;i<Myosin.myosinCt;i++) {
			double myoDuty = Myosin.theMyosins[i].getDutyRatio();
			if (myoDuty != 0) {
				aveDutyRatio += myoDuty;
				myosConsidered++;
			}
		}
		if(myosConsidered > 0) {
			return aveDutyRatio/myosConsidered;
		}
		return 0;
	}

	public void setStateNONE () {
		atpState = NONE;
	}

	public void setStateATP () {
		atpState = ATP;
	}

	public void setStateADPPi () {
		atpState = ADPPi;
	}

	public void setStateADP () {
		atpState = ADP;
	}

	private void updateBindingTimes(double dT) {
		if(isFree()) {
			if (closeEnoughToBind) {
				timeOff += dT;
			}
			return;
		}
		timeOn += dT;
	}

	/* Calculate monomer binding sites and force sites with respect to new monomer location. */
	private void updateSites() {

		if(isFree()) return;

		bindingSite.set(boundMon.getLocation());
		if(boundState == BOUND) {
			forceSite.add(bindingSite, tangOffset, boundMon.myFilament.uVect);
			forceSite.wrapPoint();
			//forceVector.getVector(attachPt,forceSite);
			//forceVector.uVect();
			//headJoint.add(attachPt,Constants.myoStalkLength,forceVector);
		}
		else { // state == BOUND_PS
			forceSite.add(bindingSite, tangOffset + myoStepSize,boundMon.myFilament.uVect);
			forceSite.wrapPoint();
			//forceVector.getVector(attachPt,forceSite);
			//forceVector.uVect();
			//headJoint.add(attachPt,Constants.myoStalkLength,forceVector);
		}
	}

	/* If binding angle between myosin and filament is too large then release */
	private boolean tooLargeBindingAngle(Point2D forceUVec) {
		Point2D filVec = boundMon.myFilament.uVect;
		double dotAng = Point2D.dot(forceUVec,filVec);
		double theAng = Math.acos(dotAng);
		if (theAng > myoMaxBindAngle) {
			releaseMon();
			//System.out.println("broke myo bond 'cause angle was too large");
			return true;
		}
		return false;
	}

	/* calculate and apply motor forces. */
	private void doForces() {
		/** This version of doForces calculates two separate spring forces based on the orthogonal drop
		 *  from the myosin attachPt to the filament.  One works to keep that distance constant while
		 *  the other represents the on-filament-axis constraint and powerstroke.
		 */
		if(boundState == FREE) return;
		if(boundMon == null) { releaseMon(); System.out.println("Caught null boundMon in doForces(). MonomerCt = " + Monomer.monomerCt);}
		if(boundState == FREE) return;

		// Check binding angle constraint
		bindingUVec.getVector(attachPt,forceSite);
		bindingUVec.uVect();
		//if (tooLargeBindingAngle(bindingUVec)) { return;} //leave if myosin released due to binding angle

		// Relevant distances and vectors for the two forces
		orthoPt = pointAndLine(attachPt,boundMon.myFilament.pEnd,boundMon.myFilament.bEnd);
		double orthoDist = Point2D.getDistance(orthoPt, attachPt) - myoStdOrthDist;
		if (orthoDist < 0) { orthoDist = 0; }

		double tangDist = Point2D.getDistance(orthoPt, forceSite);

		forceVecOrtho.getVector(attachPt, orthoPt);
		forceVecOrtho.uVect();
		forceVecTang.getVector(orthoPt, forceSite);
		forceVecTang.uVect();

		uVecTang.copy(forceVecTang); // a copy of the tangential force uVec
		uVecFil.copy(boundMon.myFilament.uVect); // copy of current filament's uVec

		// force orthogonal to filament
		forceMagOrtho = myoOrthoSpringConstant*orthoDist;
		forceVecOrtho.scale(forceMagOrtho);

		myHolder.addForce(forceVecOrtho, attachPt);
		forceVecOrtho.scale(-1);
		boundMon.addForce(forceVecOrtho);

		// force on axis of filament
		forceMagTang = myoStrokeSpringConstant*tangDist;
		//if (forceMagTang > Myosin.myoMaxSingleForce) { forceMagTang = Myosin.myoMaxSingleForce; }
		if (!neglectNegativeMyoForces || Point2D.dot(uVecTang, uVecFil) > 0){
			forceVecTang.scale(forceMagTang);
		}
		else if (neglectNegativeMyoForces && Point2D.dot(uVecTang,uVecFil) <= 0){
			forceVecTang.scale(0);
		}

		// store forceMagTangSigned
		if (Point2D.dot(uVecTang, uVecFil) > 0) {
			forceMagTangSigned = forceMagTang;
		} else {
			forceMagTangSigned = -forceMagTang;
		}

		//myHolder.addForce(forceVecTang, attachPt);
		myHolder.addForce(forceVecTang);
		forceVecTang.scale(-1);
		boundMon.addMyoForce(forceVecTang, this);
		//boundMon.addMyoConnect(myHolder);

	}

	private void zeroForces() {
		forceMagTang = 0;
		forceMagOrtho = 0;
	}

	public void usedReleaseRateMult(int whichOne){
		if (whichOne==1){
			usedReleaseRateMult=adpReleaseRateMult;
		}
		else if (whichOne==2){
			usedReleaseRateMult=secondADPReleaseRateMult;
		}
		Random releaseRateMultGen = new Random((long)Math.random());
		usedReleaseRateMult=usedReleaseRateMult+stDevADPReleaseRateMult*2*releaseRateMultGen.nextGaussian()-stDevADPReleaseRateMult;
	}

	public double getTangDist(){
		Point2D tangPt = pointAndLine(attachPt,boundMon.myFilament.pEnd,boundMon.myFilament.bEnd);
		double tangDist = Point2D.getDistance(orthoPt, forceSite);
		return tangDist;
	}

	public double getforceMagTangSigned(){
		return forceMagTangSigned;
	}

	public Point2D getTotalForce() {
		Point2D p = new Point2D();
		p.add(forceVecOrtho,forceVecTang);
		p.scale(-1);
		return p;
	}


	/*private void doForces() {
		if(boundState == FREE) return;
		forceVector.getVector(attachPt,forceSite);
		forceVector.uVect();
		if (tooLargeBindingAngle(forceVector)) { return; } //leave if myosin released due to binding angle
		bondStretch = Point.getDistance(forceSite, attachPt) - minBindingDistance;  // calculate distance between center of particles
		if (bondStretch < 0) { bondStretch = 0; }
		if (!Point.pointOK(forceVector)) {
			System.out.println ("forceVector in Myosin.applyForce() is NaN");
			System.out.println ("bondStretch = " + bondStretch);
			System.out.println ("mFAttachmentSite is " + attachPt.reportVals());
			System.out.println ("forceSite is " + forceSite.reportVals());
			System.out.println ("boundMon is " + boundMon);
			System.out.println ("boundMon.myFilament is " + boundMon.myFilament + " with cm at " + boundMon.myFilament.cm.reportVals());
		}
		forceMag = Constants.myoSpringConstant*bondStretch;
		forceAv = (1-Constants.myoForceAvFactor)*forceAv + Constants.myoForceAvFactor*forceMag;
		forceVector.scale(forceMag);
		myHolder.addForce(forceVector, attachPt);
		forceVector.scale(-1);
		boundMon.addForce(forceVector);
	}*/

	public double getAngle(){
		oriUVect.x = orientation*myHolder.uVect.x;
		oriUVect.y = orientation*myHolder.uVect.y;

		return Point2D.getAngle(boundMon.myFilament.uVect, oriUVect);
	}

	public Point2D getAim(){
		oriUVect.x = orientation*myHolder.uVect.x;
		oriUVect.y = orientation*myHolder.uVect.y;
		return oriUVect;
	}

	public void setConnectedFromHolder () {
		isConnected = myHolder.isConnected;
		connectedLevel = myHolder.connectedLevel;
	}

	public void setToHolderCluster () {
		clusterId = myHolder.clusterId;
		clusterTime = myHolder.clusterTime;
	}

	public void traceTendril (ForceTendril myTendril) {}

	public void setToHolderTendril () {
		clusterId = myHolder.clusterId;
		tendrilTime = myHolder.tendrilTime;
	}

	public void drawYourself (Graphics g, double scale, double [] offset) {
		if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) { return;}

		if (myHolder.offCortex()) { return; }  // don't render if my myosin minifil is off-cortex
		Graphics2D g2d=(Graphics2D)g;
		Composite originalComposite = g2d.getComposite();

		int mfAttachX = attachPt.getPixX();
		int mfAttachY = attachPt.getPixY();
		int pixelDiameter = (int) (0.5*headDiameter*scale);
		int pixelRadius = pixelDiameter/2;

		if (isFree()) {
			g2d.setPaint(freeColor);
		} else {
			g2d.setPaint(Color.yellow);
		}
		g2d.drawOval(mfAttachX-pixelRadius, mfAttachY-pixelRadius, pixelDiameter, pixelDiameter);
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
		g2d.fillOval(mfAttachX-pixelRadius, mfAttachY-pixelRadius, pixelDiameter, pixelDiameter);
		g2d.setComposite(originalComposite);

		if (!isFree()){
			g2d.setPaint(Color.cyan);
			int xBSite = forceSite.getPixX(); //bindingSite.getPixX();
			int yBSite = forceSite.getPixY(); //bindingSite.getPixY();
			int xHeadJoint = headJoint.getPixX();
			int yHeadJoint = headJoint.getPixY();

			g2d.drawOval(xBSite-pixelRadius, yBSite-pixelRadius, pixelDiameter, pixelDiameter);
			g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
			g2d.fillOval(xBSite-pixelRadius, yBSite-pixelRadius, pixelDiameter, pixelDiameter);
			g2d.setComposite(originalComposite);

			/*// Draw the NECK... line rendering that works with crossing boundaries
			if (boundState == BOUND) {
				g2d.setPaint(Color.cyan);
			} else {
				g2d.setPaint(Color.blue);
			}
			g2d.drawLine(mfAttachX, mfAttachY, xHeadJoint, yHeadJoint);
			//g2d.drawLine(xBSite,yBSite,xBSite-neckLineX,yBSite-neckLineY);

			// Draw the HEAD... line rendering that works with crossing boundaries
			g2d.setPaint(Color.cyan);
			g2d.drawLine(xHeadJoint, yHeadJoint, xBSite, yBSite);
			//g2d.drawLine(xBSite,yBSite,xBSite-headLineX,yBSite-headLineY);
			 */
			/** Draw straight line from attachPt to bindingSite... line rendering that works with crossing boundaries */
			if (boundState == BOUND) {
				g2d.setPaint(Color.cyan);
			} else {
				g2d.setPaint(Color.blue);
			}
			renderVec.getVector(attachPt,forceSite); //bindingSite);
			int monX = forceSite.getPixX();
			int monY = forceSite.getPixY();
			g2d.drawLine(mfAttachX, mfAttachY, monX, monY);
			//g2d.drawLine(xBSite, yBSite, -monX, -monY);

			/** Draw orthogonal drop to actin from attachPt... for testing */
			/*g2d.setPaint(Color.orange);
			renderVec.getVector(attachPt,orthoPt);
			int orthX = mfAttachX + renderVec.getPixX();
			int orthY = mfAttachY + renderVec.getPixY();
			g2d.drawLine(mfAttachX, mfAttachY, orthX, orthY);
			//g2d.drawLine(orthoPtX, orthoPtY, -orthX, -orthY);
			 */
		}
	}

	public boolean isFree (){
		//monLinkSanityCk();
		return boundMon == null;
	}

	public boolean canBind () {

		if (myHolder instanceof MyosinMiniFilament) {
			MyosinMiniFilament myMini = (MyosinMiniFilament)myHolder;
			if (iAmAnAMyo) {
				if (myMini.zPosB != 0) { return false; }
			} else {
				if (myMini.zPosA != 0) { return false; }
			}
		}
		return isFree();
	}


	public void bindMon (Monomer mon) {
		if (atpState == ATP) { return; } // cannot bind if head is ATP
		if (!mon.isBindingSite) { return; } // cannot bind if monomer isn't marked as a binding site
		releaseMon();
		mon.bindMyo(this);
		setBoundState(Myosin.BOUND);
		setInitialMonOffset();
		lastTimeOn = Sim2D.simulationTime;
		timeOff_tmp = lastTimeOn - lastTimeOff;

	}

	public void releaseMon () {
		if(boundMon!=null){
			boundMon.releaseMyo();
			zeroForces();
		}
		setBoundState(Myosin.FREE);
		boundMon = null;
	}

	public void setInitialMonOffset () {
		orthoPt = pointAndLine(attachPt,boundMon.myFilament.pEnd,boundMon.myFilament.bEnd);
		bindingSite.set(boundMon.getLocation());
		double pEndToOrtho = Point2D.getDistance(boundMon.myFilament.pEnd, orthoPt);
		double pEndToBSite = Point2D.getDistance(boundMon.myFilament.pEnd, bindingSite);
		tangOffset = pEndToOrtho - pEndToBSite;  // signed appropriately

		orthoOffset = Point2D.getDistance(attachPt, orthoPt);
	}

	public void monLinkSanityCk () {
		if (boundMon != null) {
			if (boundMon.boundMyo != this) { releaseMon(); }
		}
	}

	public static void reportBoundMyos () {
		int boundMyos = 0;
		for (int i=0;i<myosinCt;i++) {
			if (!theMyosins[i].isFree()) { boundMyos++; }
		}
		System.out.println ("Myosins bound to monomers = " + boundMyos);
	}

	public static void addMyosin (Myosin newMyosin) {
		newMyosin.myMyosinNumber = myosinCt;
		theMyosins[myosinCt] = newMyosin;
		myosinCt ++;
	}

	public void die() {
		boundMon = null;
		myHolder = null;
	}

	public void removeMe () {
		super.removeMe();
		int swapId = myMyosinNumber;
		theMyosins[swapId] = theMyosins[myosinCt-1];
		theMyosins[swapId].myMyosinNumber = swapId;
		myosinCt --;
	}

	/** Returns the zero-force expected lifetime (s) of a myosin crossbridge in powerstroke based on transition rates */
	public static double getCrossBridgeLifetime () {
		double expLife = 1/myoOnFilADP_NoneRate + 1/atpOnMyoRate + 1/myoUniformReleaseProb;
		return expLife;
	}

	static public double getCrossbridgeForce() {
		return myoStrokeSpringConstant*myoStepSize;
	}

	/** return point on the line orthogonal to p3 */
	public static Point2D pointAndLine (Point2D p3, Point2D p1, Point2D p2) {
		double x3Mx1 = Point2D.getDiffX(p1, p3);
		double y3My1 = Point2D.getDiffY(p1, p3);
		double x2Mx1 = Point2D.getDiffX(p1, p2);
		double y2My1 = Point2D.getDiffY(p1, p2);
		double p1p2DistSqrd = x2Mx1*x2Mx1 + y2My1*y2My1;

		double u = (x3Mx1*x2Mx1 + y3My1*y2My1)/p1p2DistSqrd;

		double x = p1.x + u*x2Mx1;
		double y = p1.y + u*y2My1;

		Point2D orthoPoint = new Point2D(x,y);
		orthoPoint.wrapPoint();

		return orthoPoint;
	}

	/******* MYOSIN  ********/
	static public final String DISSOCIATE_ADP = "Dissociate_ADP";
	static public final String DISSOCIATE_ADP_CATCH_SLIP = "Dissociate_ADP_Catch_Slip";
	static public final String DISSOCIATE_ADP_FORCE_BASED = "Dissociate_ADP_Force_Based";

	static public String adpReleaseMode = DISSOCIATE_ADP;

	/** motor head biochemical transitions */
	static double atpOnMyoRate;
	static double myoOnFilATP_ADPPiRate;
	static double myoOffFilATP_ADPPiRate; // 10 is 2%, 20 is 4% duty ratio, 100 is 16%, 1000 is 55%, 10000 is 71%,
	static double myoOnFilADPPi_ADPRate;
	static double myoOffFilADPPi_ADPRate;
	static double myoOnFilADP_NoneRate;

	/**multiply ADP release rate by a constant**/
	static public double adpReleaseRateMult;

	/**use a second ADP release rate within the same myosin filament**/
	static public double secondADPReleaseRateMult;

	/** stiffness of the myosin head/F-actin bond.(pN/nm) */
	public static double myoStrokeSpringConstant;

	static double myoOrthoSpringConstant;

	/** maximum binding angle for myosin... release above this */
	static double myoMaxBindAngle;

	/** Diameter of a myosin head (nm) */
	static double myoRadius;

	/** length of the stalk which connects myosin head to the minifilament (nm). */
	static double myoStalkLength;

	/** How far the myosin head steps along the minifilament (nm). */
	static public double myoStepSize;

	/** Whether or not to use force-based release kinetics. */
	static boolean myoUseForceBasedRelease;

	/** Maximun force in a single myosin crossbridge */
	static double myoMaxSingleForce;

	/** Probability of force-independent release. */
	static double myoUniformReleaseProb;

	/** Probability of force-independent binding. */
	static public double myoUniformBindingProb;

	static double myoBindingProbHalfWidth;

	/** myosin force independent rate of power stroking */
	static double myoPowerStrokeProb;

	/** The Cosine of the angle between the actin filament axis and the minifilament axis
	 * must be greater than this number for a head to bind. */
	static public double myoAngleBias;

	/** fraction of current myosin force to add to moving average. Should be between 0
	 * and 1. A value of 1 implies no force memory.  Lower values -> longer memory. */
	static public double myoForceAvFactor;

	/** effective radius (nm) for minifilament-minfilament repulsion. */
	static public double myoCollisionRadius;

	/** The maximaum allowable interval between myosin binding checks.  This
	 * limits the interval set by adaptive mechanism.. */
	static public double myoBindCheckMaxInterval;

	/** For adaptive adjustment of intervals at which we compute myosin bond formation.
	 * Defines the TARGET max probability of forming a new bond in the adaptive
	 * interval calculated over all myosins. */
	static public double myoTargetBindingProbPerTimestep;

	/** For adaptive adjustment of intervals at which we compute myosin bond formation.
	 * Defines the max ALLOWABLE probability of forming a new bond in the adaptive
	 * interval calculated over all myosins. */
	static public double myoMaxBindingProbPerTimestep;

	/** For adaptive adjustment of intervals at which we compute myosin bond formation.
	 * Defines the TARGET move for a filament monomer during
	 * one adaptive timestep interval
	 **/
	static public double targetMonomerMoveFraction;

	/** For adaptive adjustment of intervals at which we compute crosslinker bond formation.
	 * Defines the MAX move allowed for a filament monomer during
	 * one adaptive timestep interval
	 **/
	static public double maxMonomerMoveFraction;

	static public double [] cumBindingProbabilities;
	
	/**remove internal drag*/
	static public boolean neglectNegativeMyoForces;
	
	/**add variation to the release rate across motors*/
	static public double stDevADPReleaseRateMult;

	static String className = new String("main.Myosin");

	static void updateMyoBindingProbabilities()  {
		double maxD = myoBindingProbHalfWidth*Math.sqrt(-2.0*Math.log(0.3));
		int maxI = (int) maxD;
		cumBindingProbabilities = new double[2*maxI+1];
		double sum = 0;
		for(int i = 0; i < maxI; i++) {
			sum += normal(Actin.monLength*(maxI-i),myoBindingProbHalfWidth);
			cumBindingProbabilities[i] = sum;
		}
		sum += normal(0,myoBindingProbHalfWidth);
		cumBindingProbabilities[maxI] = sum;
		for(int i = 1; i <= maxI; i++) {
			sum += normal(Actin.monLength*(i),myoBindingProbHalfWidth);
			cumBindingProbabilities[maxI+i] = sum;
		}
		for(int i = 0; i < 2*maxI+1; i++) {
			cumBindingProbabilities[i] /= sum;
		}
	}

	static double normal(double x,double sigma) {
		return Math.exp(-0.5*x*x/(sigma*sigma))/(sigma*Math.sqrt(2.0*Math.PI));
	}

	static void staticInit() {

		Parameters.addParameterSetListener (
				new ParameterSetListener() {
					public void parametersChanged() {
						headDiameter = 2.0*myoRadius;
						minBindingDistance = myoStalkLength + headDiameter;
						myoStdOrthDist = headDiameter+myoStalkLength;
					}
				}
		);
		Parameters.addParameter(className, "myoTargetBindingProbPerTimestep",0.2,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoTargetBindingProbPerTimestep = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoMaxBindingProbPerTimestep",0.25,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoMaxBindingProbPerTimestep = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "adpReleaseMode",DISSOCIATE_ADP,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				String s = p.getStringValue();
				if(s.equalsIgnoreCase(DISSOCIATE_ADP)) {
					adpReleaseMode = s;
				}
				else if(s.equalsIgnoreCase(DISSOCIATE_ADP_FORCE_BASED)) {
					adpReleaseMode = s;
				}
				else if(s.equalsIgnoreCase(DISSOCIATE_ADP_CATCH_SLIP)) {
					adpReleaseMode = s;
				}
				else {
					System.out.println("Myosin.parameterChanged: got bad value for adpReleaseMode = " + s);
					adpReleaseMode = DISSOCIATE_ADP;
				}
			}
		}
		);
		Parameters.addParameter(className, "myoBindCheckMaxInterval",0.1,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoBindCheckMaxInterval = p.getDoubleValue();
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
		Parameters.addParameter(className, "adpReleaseRateMult",1,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				adpReleaseRateMult= p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "secondADPReleaseRateMult",1,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				secondADPReleaseRateMult= p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "stDevADPReleaseRateMult",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				stDevADPReleaseRateMult= p.getDoubleValue();
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
		Parameters.addParameter(className, "myoUseForceBasedRelease",false,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoUseForceBasedRelease = p.getBooleanValue();
			}
		}
		);
		Parameters.addParameter(className, "atpOnMyoRate",10000,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				atpOnMyoRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoOnFilATP_ADPPiRate",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOnFilATP_ADPPiRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoOffFilATP_ADPPiRate",100,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOffFilATP_ADPPiRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoOnFilADPPi_ADPRate",1000,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOnFilADPPi_ADPRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoOffFilADPPi_ADPRate",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOffFilADPPi_ADPRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoOnFilADP_NoneRate",1000,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOnFilADP_NoneRate = p.getDoubleValue();
			}
		}
		);
		// Spring constant of 0.7pN/nm from Veigel et al Biophys J. 1998
		Parameters.addParameter(className, "myoStrokeSpringConstant",0.7,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoStrokeSpringConstant = p.getDoubleValue();
			}
		}
		);
		// Orthogonal spring of 0.2pN/nm is somewhat random
		Parameters.addParameter(className, "myoOrthoSpringConstant",0.2,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOrthoSpringConstant = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoMaxBindAngle",3*Math.PI/4,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoMaxBindAngle = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoRadius",5,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoRadius = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoStalkLength",15,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoStalkLength = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoStepSize",5.5,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoStepSize = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoMaxSingleForce",5,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoMaxSingleForce = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoUniformReleaseProb",2000,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoUniformReleaseProb = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoUniformBindingProb",30,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoUniformBindingProb = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoBindingProbHalfWidth",10,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoBindingProbHalfWidth = p.getDoubleValue();
				updateMyoBindingProbabilities();
			}
		}
		);
		Parameters.addParameter(className, "myoPowerStrokeProb",1000,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoPowerStrokeProb = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoAngleBias",0.3,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoAngleBias = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "neglectNegativeMyoForces",false,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				neglectNegativeMyoForces = p.getBooleanValue();
			}
		}
		);
		Parameters.addParameter(className, "myoForceAvFactor",0.5,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoForceAvFactor = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoCollisionRadius",50.0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoCollisionRadius = p.getDoubleValue();
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
				stepAllMyosins(dT);
			}
			public void reset() {
				for(int i = 0; i < myosinCt; i++) {
					theMyosins[i].die();
					theMyosins[i] = null;
				}
				myosinCt = 0;
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


