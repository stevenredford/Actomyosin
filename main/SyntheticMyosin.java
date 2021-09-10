package main;

/*
	A rectangular array of myosin filaments. Based loosely on Verkhovsky AB, Svitkina TM, Borisy GG JCB 1995.
 */

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Random;

import collision.*;
import util.*;
import gui.*;
import parameters.*;

public class SyntheticMyosin extends MyosinHolder {

	static CollisionDetector myCollisionDetector = new MiniFilamentCollisionDetector();

	/** array holding all MiniFilaments */
	static public MyosinMiniFilament [] theMiniFilaments = new MyosinMiniFilament [10000];

	/** Current number of minifilaments. */
	static public int miniFilamentCt = 0;	// initialize the LongMiniF count

	/** Color for drawing. */
	Color MiniFilamentColor = Color.LIGHT_GRAY;
	Color fadedMiniFilamentColor = new Color(0.1f,0.1f,0.1f);

	/** color for jumping. */
	Color jumpingColor = Color.DARK_GRAY;

	/** Unique index of this minifilament in the array of all minifilaments. */
	public int myMiniFilamentNumber;

	/** indicates if this minifilament spans one or more boundaries */
	boolean wrapped = false;

	/* First end point of the filament, calculated based on length, centermass, and uVect */
	Point2D end1 = new Point2D();

	/* Second end point of the filament, calculated based on length, centermass, and uVect */
	Point2D end2 = new Point2D();

	Point2D A1 = new Point2D();
	Point2D A2 = new Point2D();
	Point2D B1 = new Point2D();
	Point2D B2 = new Point2D();

	Point2D aLoad = new Point2D();  // resistive load

	double minifilamentStallForce = -1;

	public double [] axialForces;
	public double [] avgInternalStresses;

	public int [] internalStressCts;

	public boolean pinned = false;

	int numBoundMyosins = 0;

	public Point2D initialPos;


	/** Array of myosin heads attached to "A" end of minifilament. */
	Myosin [] myBmyosins;

	/** Array of myosin heads attached to "B" end of minifilament. */
	Myosin [] myAmyosins;

	// CONSTRUCTOR *************************************************************************************************

	public MyosinMiniFilament (double initX, double initY, boolean pinned) {
		super(initX,initY);

		this.pinned = pinned;

		initialPos = new Point2D(cm);
		centroid.set(initialPos);)

		createMyosinHeads (nMyosinHeads);
		addMiniFilament(this);
		evaluateProperties();
		initStressTracking();
	}

	public void restoreInitialPosition() {
		moveTo(initialPos);
	}

	public MyosinMiniFilament (double initX, double initY, double uVectX, double uVectY) {
		this(initX,initY,uVectX,uVectY,false);
	}

	static void setNMyosinHeads(int n) {
		nMyosinHeads = 2*((n+1)/2);  // nMyosinHeads should always be even and there should be at least 2 heads at either end
	}

	public void createMyosinHeads(int n){
		zPosB = 0;//(int)(Math.random()*maxZPos);
		zPosA = 0;//(int)(Math.random()*maxZPos);


		myBmyosins = new Myosin [n];
		myAmyosins = new Myosin [n];
		myMyosins = new Myosin[myBmyosins.length+myAmyosins.length];
		if (useSecondADPReleaseRateMult){
			pickPositionsForSecondADPReleaseRate();
		}
		for(int i=0;i<myBmyosins.length;i++){

			myBmyosins[i]=new Myosin(this);
			if (useSecondADPReleaseRateMult){
				myBmyosins[i].usedReleaseRateMult(bPositionsForSecondADPReleaseRate[i]+1);
			}
			else {
				myBmyosins[i].usedReleaseRateMult(1);
			}
			myMyosins[i]=myBmyosins[i];
			myBmyosins[i].orientation = 1;
			myBmyosins[i].iAmAnAMyo = false;
		}

		for(int i=0;i<myAmyosins.length;i++){
			myAmyosins[i]=new Myosin(this);
			if (useSecondADPReleaseRateMult){
				myAmyosins[i].usedReleaseRateMult(aPositionsForSecondADPReleaseRate[i]+1);
			}
			else {
				myAmyosins[i].usedReleaseRateMult(1);
			}
			myMyosins[i+myBmyosins.length]=myAmyosins[i];
			myAmyosins[i].orientation = -1;
			myAmyosins[i].iAmAnAMyo = true;
		}
		for (int i=0;i<myBmyosins.length;i++) {
			myBmyosins[i].setInitialPositionRelativeToMiniFilament(bOffsets[i]);
		}
		for (int i=0;i<myAmyosins.length;i++) {
			myAmyosins[i].setInitialPositionRelativeToMiniFilament(aOffsets[i]);
		}

		moveMyosins();
	}
	int [] aPositionsForSecondADPReleaseRate =new int[nMyosinHeads];
	int [] bPositionsForSecondADPReleaseRate =new int[nMyosinHeads];
	public void pickPositionsForSecondADPReleaseRate(){

		for (int i=0; i<nHeadsSecondADPReleaseRateMult; i++){
			double randomNumber=Math.random();
			for (int j=0; j<nMyosinHeads; j++){
				if (randomNumber<(j+1.0)/nMyosinHeads){
					if (aPositionsForSecondADPReleaseRate[j]==0){
						aPositionsForSecondADPReleaseRate[j]=1;
						break;
					}
					else if (aPositionsForSecondADPReleaseRate[j]==1){
						i--;
						break;
					}
				}
			}
		}
		for (int i=0; i<nHeadsSecondADPReleaseRateMult; i++){
			double randomNumber=Math.random();
			for (int j=0; j<nMyosinHeads; j++){
				if (randomNumber<(j+1.0)/nMyosinHeads){
					if (bPositionsForSecondADPReleaseRate[j]==0){
						bPositionsForSecondADPReleaseRate[j]=1;
						//System.out.println(j);
						break;
					}
					else if (bPositionsForSecondADPReleaseRate[j]==1){
						i--;
						break;
					}
				}
			}
		}
	
	}

	public MyosinMiniFilament () {
		this(0,0,1,0);
	}

	public MyosinMiniFilament (double initX, double initY) {
		this(initX,initY,1,0);
	}

	static public void fillMFilMesh() {
		MyosinMiniFilament m;
		for (int i = 0; i < miniFilamentCt; i++) {
			m = theMiniFilaments[i];
			Mesh.MFIL_MESH.addPointToMesh(m.myMiniFilamentNumber,m.cm,Myosin.myoCollisionRadius);
		}
	}

	static public void stepAllMinifilaments(double dT) {
		for (int i = 0; i < miniFilamentCt; i++) {
			if (!theMiniFilaments[i].removeMe) { theMiniFilaments[i].step(dT); }
		}
	}

	static public void turnResistiveLoadOn() {
		resistiveLoadOn = true;
	}

	static public void turnResistiveLoadOff() {
		resistiveLoadOn = false;
	}

	/** Resistive load... simulates loads to reproduce Howard Fig 16.5 */
	public void resistiveLoad (double load) {
		aLoad.copy(uVect);
		aLoad.scale(-load);
		addForce(aLoad);
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
	static private void doMinifilamentTurnover() {

		double turnoverRate = Sim2D.getDeltaT()*miniFilamentTurnoverRate*miniFilamentCt;
		int n = (int)Math.floor(turnoverRate);
		for(int i = 0; i < n; i++) {
			minifilamentTurnover();
		}

		double residual = turnoverRate - n;
		if(Math.random() < residual)  {
			minifilamentTurnover();
		}
	}

	static private void minifilamentTurnover() {
		int choice = (int)Math.floor(Math.random()*miniFilamentCt);
		MyosinMiniFilament mf = theMiniFilaments[choice];
		mf.releaseAll();
		mf.die();
		mf.removeMe();
		System.out.println("minifilament dies");
	}

	static private void doMinifilamentRecruitment() {
		double recruitmentRate = Sim2D.getDeltaT()*miniFilamentRecruitmentRate*Sim2D.xDimension*Sim2D.yDimension;
		int n = (int)Math.floor(recruitmentRate);
		for(int i = 0; i < n; i++) {
			minifilamentRecruitment();
		}

		double residual = recruitmentRate - n;
		if(Math.random() < residual)  {
			minifilamentRecruitment();
		}
	}

	static private void minifilamentRecruitment() {

		double initX = Math.random()*Sim2D.xDimension;
		double initY = Math.random()*Sim2D.yDimension;
		double randomAng = 2*Math.PI*Math.random();
		MyosinMiniFilament.makeMiniFilament(initX, initY, Math.cos(randomAng), Math.sin(randomAng));
		System.out.println("minifilament recruited");

	}

	public void step (double dT) {

		stressCalculations();

		//System.out.println ("zPosA,zPosB = " + zPosA + "," + zPosB);

		//***********************************************************************************************
		// PARALLEL DIFFUSION FORCES - moving in the direction of uVect
		// you are going to diffuse with this step size in the direction of uVect
		diffStepPar = Math.sqrt(2*curDPar*dT)* generator.nextGaussian();

		// Convert the step in parallel direction to x,y components of forces in the universal coordinate system
		brownianForceParallel.set((myosinfilTransGamma*diffStepPar) / dT,uVect);

		//***********************************************************************************************
		//PERPENDICULAR DIFFUSION FORCES - moving in the direction of upVect
		diffStepPerp = Math.sqrt(2*curDPerp*dT)* generator.nextGaussian();

		brownianForcePerpendicular.set((myosinfilPerpGamma*diffStepPerp) / dT,upVect);

		//***********************************************************************************************
		//TOTAL TRANSLATIONAL FORCES
		//calculate the total force acting on the center of mass


		if(resistiveLoadOn) {
			double l = sizeOfResistiveLoad();
			if(getNumBoundMyosins() >0) resistiveLoad(l); // add a constant resistive load
		}

		if (Sim2D.brownianOn) {
			totalTranslationalForce.inc(brownianForceParallel);
			totalTranslationalForce.inc(brownianForcePerpendicular);
		}

		velocity.set(1/myosinfilTransGamma,totalTranslationalForce);	//WHAT DRAG COEFFICIENT TO USE?

		// Move the center of mass.
		if (!pinned) cm.inc(dT,velocity);

		cm.wrapPoint();


		//************************************************************************************************
		// ROTATIONAL DIFFUSION FORCES
		diffStepRot = Math.sqrt(2*curDRot*dT)* generator.nextGaussian();

		brownianTorque = myosinRotGam * diffStepRot / dT;

		if (Sim2D.brownianOn) { torque += brownianTorque; }

		angularVel =  torque / myosinRotGam;

		//************************************************************************************************
		// TOTAL ANGULAR FORCES
		deltaTheta = angularVel * dT;		//contribution to change in theta from diffusion
		if (pinned) deltaTheta = 0;

		thetaNOW = Math.atan2 (uVect.y, uVect.x);
		thetaNEXT = thetaNOW +deltaTheta;

		uVect.x = Math.cos (thetaNEXT);
		uVect.y = Math.sin (thetaNEXT);
		uVect.uVect();
		upVect.orthogonalVector(uVect);

		resetEnds();
		// Move the Myosins
		moveMyosins();

		reset();
	}


	private void resetEnds() {
		// calculate the position of the ends
		end1.add(cm,-halfLength,uVect);
		end2.add(cm,halfLength,uVect);
		boolean end1Wrapped = end1.wrapPoint();
		boolean end2Wrapped = end2.wrapPoint();
		if (end1Wrapped | end2Wrapped) {
			wrapped = true;
		} else { wrapped = false; }
	}

	public void rotateTo(double theta) {
		uVect.x = Math.cos (theta);
		uVect.y = Math.sin (theta);
		uVect.uVect();
		upVect.orthogonalVector(uVect);
		resetEnds();
		moveMyosins();
	}


	public void moveTo(Point2D p) {
		cm.x = p.x; cm.y = p.y;
		end1.add(cm,-halfLength,uVect);
		end2.add(cm,halfLength,uVect);

		cm.wrapPoint();
		resetEnds();
		moveMyosins();
	}

	public void move(Point2D inc) {
		cm.x += inc.x; cm.y += inc.y;
		end1.add(cm,-halfLength,uVect);
		end2.add(cm,halfLength,uVect);

		cm.wrapPoint();
		resetEnds();
		moveMyosins();
	}

	public void setDiffMultipliers (double transScale, double rotScale) {
		curDPar = transScale*diffPar;
		curDPerp = transScale*diffPerp;
		curDRot = rotScale*diffRot;
	}

	// Determination of buckling of internal segments
	public void stressCalculations () {
		if(Actin.monitorStress) {
			findAxialForces();
			findInternalForces();
		}
	}

	public boolean allMyosFree () {
		for (int i=0;i<myBmyosins.length;i++) {
			if (!myBmyosins[i].isFree()) return false;
		}
		for (int i=0;i<myAmyosins.length;i++) {
			if (!myAmyosins[i].isFree()) return false;
		}
		return true;
	}

	public boolean myosFreeA () {
		for (int i=0;i<myAmyosins.length;i++) {
			if (!myAmyosins[i].isFree()) return false;
		}
		return true;
	}

	public boolean myosFreeB () {
		for (int i=0;i<myBmyosins.length;i++) {
			if (!myBmyosins[i].isFree()) return false;
		}
		return true;
	}

	public boolean offCortex () {
		if (zPosB != 0 & zPosA != 0) { return true; }
		return false;
	}


	public void moveMyosins () {
		numBoundMyosins = 0;
		Point2D temp=new Point2D();
		for (int i=0;i<myBmyosins.length;i++) {
			temp.x=cm.x
			+myBmyosins[i].attachPtLocal.x*uVect.x
			+myBmyosins[i].attachPtLocal.y*upVect.x;
			temp.y=cm.y
			+myBmyosins[i].attachPtLocal.x*uVect.y
			+myBmyosins[i].attachPtLocal.y*upVect.y;

			//wrap
			temp.wrapPoint();
			myBmyosins[i].attachPt.set(temp);
			if (myBmyosins[i].isFree()){
				myBmyosins[i].bindingSite.set(temp);
			}
			else numBoundMyosins++;
		}
		for (int i=0;i<myAmyosins.length;i++) {
			temp.x=cm.x
			+myAmyosins[i].attachPtLocal.x*uVect.x
			+myAmyosins[i].attachPtLocal.y*upVect.x;
			temp.y=cm.y
			+myAmyosins[i].attachPtLocal.x*uVect.y
			+myAmyosins[i].attachPtLocal.y*upVect.y;

			//wrap
			temp.wrapPoint();
			myAmyosins[i].attachPt.set(temp);
			if (myAmyosins[i].isFree()){
				myAmyosins[i].bindingSite.set(temp);
			}
			else numBoundMyosins++;
		}
	}


	public void reset () {
		// reset forces and velocities
		super.reset();
	}

	private void initStressTracking() {
		if(Actin.monitorStress) {
			axialForces = new double[nMyosinHeads];
			avgInternalStresses = new double[nMyosinHeads-1];
			internalStressCts = new int[nMyosinHeads-1];
		}
	}

	int AStressIndex(int i) {
		return (nMyosinHeads-1-i)/2;
	}

	int BStressIndex(int i) {
		return (nMyosinHeads+i)/2;
	}

	public int getStressIndex(double arcL) {
		int n = nMyosinHeads/2 - 1;
		int ind = 0;
		for(int i = 0; i < n; i++) {
			arcL -= myoSpacing;
			if(arcL <= 0) {
				internalStressCts[ind]++;
				return ind;
			}
			ind++;
		}
		arcL -= length;
		if(arcL <= 0) {
			internalStressCts[ind]++;
			return ind;
		}
		ind++;
		for(int i = 0; i < n; i++) {
			arcL -= myoSpacing;
			if(arcL <= 0) {
				internalStressCts[ind]++;
				return ind;
			}
			ind++;
		}
		return -1; // someting bad happened.
	}

	/** return the coordinates of the first end of exended minifilament (i.e extended to attachment point of distal-most heads)
	 * NOTE:  This ignores periodic boundary conditions and is thus useful for intersection tests. */
	public Point2D getEndPoint1() {
		Point2D p = new Point2D();
		p.add(end1,-headLength,uVect);
		return p;
	}

	/** return the coordinates of the second end of exended minifilament (i.e extended to attachment point of distal-most heads)
	 * NOTE:  This ignores periodic boundary conditions and is thus useful for intersection tests. */
	public Point2D getEndPoint2() {
		Point2D p = new Point2D();
		p.add(end1,length+headLength,uVect);
		return p;
	}


	// calculates centerpoints of stress swegemnts along this filament
	public Point2D [] getStressCenters () {
		double spacing = myoSpacing;
		int centerIndex = nMyosinHeads/2;
		Point2D [] centers = new Point2D[nMyosinHeads-1];
		Point2D p = new Point2D(cm);

		centers[centerIndex] = Point2D.Copy(cm);

		p.add(-halfLength-0.5*spacing,uVect);
		for(int i = centerIndex-1; i >= 0; i--) {
			centers[i] = new Point2D(p);
			p.add(-spacing,uVect);
		}

		p.set(cm);
		p.add(halfLength+0.5*spacing,uVect);
		for(int i = centerIndex+1; i < nMyosinHeads-1; i++) {
			centers[i] = new Point2D(p);
			p.add(spacing,uVect);
		}
		return centers;
	}

	// Load axial force components and points of application
	public void findAxialForces () {

		for(int i = 0; i < axialForces.length; i++) {
			axialForces[i] = 0;
		}
		for (int i=0;i < myAmyosins.length;i++) {
			if (!myAmyosins[i].isFree()) {
				axialForces[AStressIndex(i)] = Point2D.dot(myAmyosins[i].getTotalForce(),uVect);
			}
		}
		for (int i=0;i < myBmyosins.length;i++) {
			if (!myBmyosins[i].isFree()) {
				axialForces[BStressIndex(i)] = Point2D.dot(myBmyosins[i].getTotalForce(),uVect);
			}
		}
	}

	public void findInternalForces () {
		double APointForceSum = 0;	// force on A end
		double BPointForceSum = 0;	// force on B end

		for (int i=0;i<nMyosinHeads;i++) { BPointForceSum += axialForces[i]; }  // sum all axial forces
		for (int i=0;i<nMyosinHeads-1;i++) {

			BPointForceSum -= axialForces[i];
			APointForceSum += axialForces[i];

			if (APointForceSum > 0 && BPointForceSum < 0) {  // this is the only condition for internal compression
				if (APointForceSum > Math.abs(BPointForceSum)) {
					registerInternalForce(i,BPointForceSum);
				} else {										// register smaller of two oppposing forces
					registerInternalForce(i,-APointForceSum);
				}
			} else if (APointForceSum < 0 && BPointForceSum > 0) { // condition for internal tension
				if (Math.abs(APointForceSum) > BPointForceSum) {
					registerInternalForce(i,BPointForceSum);
				} else {										// register smaller of two oppposing forces
					registerInternalForce(i,-APointForceSum);
				}
			} else {
				registerInternalForce(i,0);  // shouldn't ever set value here... just in case though
			}
		}
	}

	public double getInternalStress(int ind) {
		return avgInternalStresses[ind];
	}

	private void registerInternalForce(int i, double f) {
		avgInternalStresses[i] = (1-Sim2D.forceAverageFraction)*avgInternalStresses[i] + Sim2D.forceAverageFraction*f;
	}

	public void evaluateProperties () {
		double naturalLogTerm = (double) Math.log(0.5*totalLength/(Myosin.myoRadius));
		double PiViscosityLength= (double) Math.PI*Constants.viscosity*totalLength;

		// double pValue = (2*radius)/length;
		// the Sigmas would normally be looked up from a table given the pValue.
		// Jacked from page 107 of Jonathan Howard's Mechanics of Motor Proteins and the Cytoskeleton
		double rotSigma = -0.66;
		double parSigma = -0.20;
		double perpSigma = 0.86;

		double pValue = Myosin.myoRadius/halfLength;

		//System.out.println();
		//System.out.println("pValue="+pValue+"  radius="+radius+"  length="+length);
		boolean found=false;
		for(int i=0;i<Drag.cylinderDragValues.length&&!found;i++){
			if(pValue>=Drag.cylinderDragValues[i][Drag.pValueIndex]){
				parSigma=Drag.cylinderDragValues[i][Drag.parIndex];
				perpSigma=Drag.cylinderDragValues[i][Drag.perpIndex];
				rotSigma=Drag.cylinderDragValues[i][Drag.rotIndex];
				found=true;
				//System.out.println("using pValue="+pValue+"  index="+i+"  parSigma="+parSigma+"  perpSigma="+perpSigma+"  rotSigma="+rotSigma);
			}
		}

		if(!found){
			int defaultIndex=Drag.cylinderDragValues.length-1;
			parSigma=Drag.cylinderDragValues[defaultIndex][Drag.parIndex];
			perpSigma=Drag.cylinderDragValues[defaultIndex][Drag.perpIndex];
			rotSigma=Drag.cylinderDragValues[defaultIndex][Drag.rotIndex];
			//System.out.println("pValue outside of range so using default");
		}


		untweakedMyosinRotationalGamma = (PiViscosityLength*totalLength*totalLength)/(3*(naturalLogTerm+rotSigma));
		untweakedMyosinfilTransGamma = (2*PiViscosityLength)/(naturalLogTerm+parSigma);
		untweakedMyosinfilPerpGamma = (4*PiViscosityLength)/(naturalLogTerm+perpSigma);

		reTweakViscosities();

		setDiffMultipliers(1,1);	// default to no multiplier for diffusion terms

		if (diffPerp > diffPar) {
			outOfPlaneDiff = outOfPlaneDiffFactor * diffPerp;
		} else {
			outOfPlaneDiff = outOfPlaneDiffFactor * diffPar;
		}
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
		myosinRotGam = untweakedMyosinRotationalGamma*tweakDragFactor;
		myosinfilTransGamma = untweakedMyosinfilTransGamma*tweakDragFactor;
		myosinfilPerpGamma = untweakedMyosinfilPerpGamma*tweakDragFactor;

		diffPar = Constants.kT / myosinfilTransGamma;
		diffPerp = Constants.kT / myosinfilPerpGamma;
		diffRot =  Constants.kT / myosinRotGam;
	}

	public void drawYourself (Graphics g, double scale, double [] offset) {
		//if (Sim2D.paintOnlyConnected & !isConnected) { return;}

		if (offCortex()) { System.out.println("off!"); return; } // no rendering if off cortex
		Graphics2D g2d=(Graphics2D)g;
		g2d.setPaint(MiniFilamentColor);
		if (Sim2D.paintConnected & connectedLevel > Sim2D.paintConnectedLevel) g2d.setPaint(fadedMiniFilamentColor);

		// draw from end1
		A1.add(end1,-halfWidth,upVect);
		A2.add(end1,length,uVect,-halfWidth,upVect);
		B1.add(end1,halfWidth,upVect);
		B2.add(end1,length,uVect,halfWidth,upVect);
		drawPerimeter(g2d);

		if (wrapped) {		// if spanning a boundary we need to draw from both ends
			// draw from end2
			A1.add(end2,-length,uVect,-halfWidth,upVect);
			A2.add(end2,-halfWidth,upVect);
			B1.add(end2,-length,uVect,halfWidth,upVect);
			B2.add(end2,halfWidth,upVect);
			drawPerimeter(g2d);
		}

		if(showInfo){
			float END1X=(float)((end1.x-offset[0])*scale);
			float END1Y=(float)((end1.y-offset[0])*scale);

			float END2X=(float)((end2.x-offset[0])*scale);
			float END2Y=(float)((end2.y-offset[0])*scale);

			float CMX=(float)((cm.x-offset[0])*scale);
			float CMY=(float)((cm.y-offset[0])*scale);

			//g2d.drawString("A1",(int)A1X,(int)A1Y);
			//g2d.drawString("A2",(int)A2X,(int)A2Y);
			//g2d.drawString("B1",(int)B1X,(int)B1Y);
			//g2d.drawString("B2",(int)B2X,(int)B2Y);

			g2d.drawString("end1",(int)END1X,(int)END1Y);
			g2d.drawString("end2",(int)END2X,(int)END2Y);

			int EndPixelRadius = (int) (10*scale);
			int EndPixelDiameter = (int) (2*10*scale);

			if(showEndsInfo){
				int A1X=A1.getPixX();
				int A1Y=A1.getPixY();

				int A2X=A2.getPixX();
				int A2Y=A2.getPixY();

				int B1X=B1.getPixX();
				int B1Y=B1.getPixY();

				int B2X=B2.getPixX();
				int B2Y=B2.getPixY();
				g2d.setPaint(Color.CYAN);
				g.fillOval((int)(A1X-EndPixelRadius), (int)(A1Y-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);
				g2d.setPaint(Color.WHITE);
				g.fillOval((int)(A2X-EndPixelRadius), (int)(A2Y-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);
				g2d.setPaint(Color.RED);
				g.fillOval((int)(B1X-EndPixelRadius), (int)(B1Y-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);
				g2d.setPaint(Color.YELLOW);
				g.fillOval((int)(B2X-EndPixelRadius), (int)(B2Y-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);

				g2d.setPaint(Color.MAGENTA);
				g.fillOval((int)(CMX-EndPixelRadius), (int)(CMY-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);

				g2d.setPaint(Color.YELLOW);
				g.fillOval((int)(END1X-EndPixelRadius), (int)(END1Y-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);

				g2d.setPaint(Color.CYAN);
				g.fillOval((int)(END2X-EndPixelRadius), (int)(END2Y-EndPixelRadius), EndPixelDiameter, EndPixelDiameter);
			}

			if(showUVect){
				int U1X=(int)((cm.x-offset[0])*scale);
				int U1Y=(int)((cm.y-offset[0])*scale);

				double lengthOfUVectDrawing=20;
				double newx=cm.x+lengthOfUVectDrawing*uVect.x;
				double newy=cm.y+lengthOfUVectDrawing*uVect.y;
				int U2X=(int)((newx-offset[0])*scale);
				int U2Y=(int)((newy-offset[0])*scale);

				g2d.setPaint(Color.cyan);
				g2d.drawLine(U1X, U1Y, U2X, U2Y);


				int UP1X=cm.getPixX();
				int UP1Y=cm.getPixY();

				newx=cm.x+lengthOfUVectDrawing*upVect.x;
				newy=cm.y+lengthOfUVectDrawing*upVect.y;
				int UP2X=(int)((newx-offset[0])*scale);
				int UP2Y=(int)((newy-offset[0])*scale);

				g2d.setPaint(Color.magenta);
				g2d.drawLine(UP1X, UP1Y, UP2X, UP2Y);
			}
		}
		/* int U1X=(int)((centermass.x-offset[0])*scale);
	    	 int U1Y=(int)((centermass.y-offset[0])*scale);

	    	 g2d.drawString(this+" myID= "+Integer.toString(myMiniFilamentNumber)+"  "+Double.toString(brownianTorque),U1X,50+myMiniFilamentNumber*20);*/
	}

	public void drawPerimeter (Graphics2D g2d) {
		int A1X=A1.getPixX();
		int A1Y=A1.getPixY();

		int A2X=A2.getPixX();
		int A2Y=A2.getPixY();

		int B1X=B1.getPixX();
		int B1Y=B1.getPixY();

		int B2X=B2.getPixX();
		int B2Y=B2.getPixY();

		perimeter.reset();
		perimeter.moveTo(A1X,A1Y);
		perimeter.lineTo(A2X,A2Y);
		perimeter.lineTo(B2X,B2Y);
		perimeter.lineTo(B1X,B1Y);
		perimeter.closePath();

		g2d.draw(perimeter);
		Composite originalComposite = g2d.getComposite();
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
		g2d.fill(perimeter);
		g2d.setComposite(originalComposite);
	}

	static public void makeMiniFilament(double initX, double initY, double uVectX, double uVectY) {
		MyosinMiniFilament mf=new MyosinMiniFilament (initX, initY, uVectX, uVectY,false);
	}

	public static void addMiniFilament (MyosinMiniFilament newMiniFilament) {
		addMyosinHolder(newMiniFilament);
		newMiniFilament.myMiniFilamentNumber = miniFilamentCt;
		theMiniFilaments[miniFilamentCt] = newMiniFilament;
		miniFilamentCt ++;
	}

	public void die()  {
		myBmyosins = null;
		myAmyosins = null;
		myMyosins = null;
	}


	public void removeMe () {
		super.removeMe();
		int swapId = myMiniFilamentNumber;
		theMiniFilaments[swapId] = theMiniFilaments[miniFilamentCt-1];
		theMiniFilaments[swapId].myMiniFilamentNumber = swapId;
		miniFilamentCt --;
	}

	public int getNumBoundMyosins() {
		return numBoundMyosins;
	}

	public int getNumAHeadsBound() {
		int cnt = 0;
		for(int i = 0; i < myAmyosins.length; i++) {
			if(!myAmyosins[i].isFree()) {
				cnt++;
			}
		}
		return cnt;
	}

	public int getNumBHeadsBound() {
		int cnt = 0;
		for(int i = 0; i < myBmyosins.length; i++) {
			if(!myBmyosins[i].isFree()) {
				cnt++;
			}
		}
		return cnt;
	}


	public Myosin getFreeMyosin() {
		if(getNumBoundMyosins() == nMyosinHeads) return null;
		Myosin m = null;
		int choice;
		while(m == null) {
			choice = (int) (nMyosinHeads*Math.random());
			if (myMyosins[choice].isFree()){
				m = myMyosins[choice];
			}
		}
		return m;
	}


	/** Predict the stall force as (# heads)*(~dutyrato at stall)*(max force per head) */
	static public double predictStallForce()  {
		return nMyosinHeads*Myosin.predictDutyRatio(Myosin.getMaxForce())*Myosin.getMaxForce();
	}

	static public int predictNumHeadsAtStall() {
		return (int) (nMyosinHeads*Myosin.predictDutyRatio(Myosin.getMaxForce()));
	}

	public double sizeOfResistiveLoad() {
		if(minifilamentStallForce > 0) {
			return myoResistiveLoad*minifilamentStallForce;
		}
		else {
			return myoResistiveLoad*predictStallForce();
		}
	}

	public double getStallForce() {
		return minifilamentStallForce;
	}

	public void setStallForce(double f) {
		minifilamentStallForce = f;
	}

	/** Number of minifilaments in this simulation. */
	static public int numInitialMinifilaments;

	/** Number of individual heads at each end of a myosin minifilament. */
	static public int nMyosinHeads;

	/** Magnitude of a resistive load for calculating force-velocity curves, expressed as a
	 * fraction of the estimated stall force, which will vary with other parameters.. */
	static public double myoResistiveLoad;

	/** Spacing of myosin heads along minifilament long-axis */
	static public double myoSpacing;

	/**Offset in placement of top myosin heads from bottom heads along filament */
	static double myoOffset;

	/** Length of a minifilament rod from end to end. */
	static public double length;

	/** Half the length of a minifilament rod. */
	static double halfLength;

	/** length that heads extend beyond the rod. */
	static double headLength;

	/** total length, rod plus heads.*/
	static public double totalLength;

	/** Width of the minifilament central rod. */
	static double width;

	/** Half the width of the central rod. */
	static double halfWidth;

	/** square of collision radius used for minifilamernt collision calculations. */
	static public double collisionRadSqrd;

	static public double miniFilamentTurnoverRate;

	static public double miniFilamentRecruitmentRate;

	static public boolean miniFilamentsTurnover;

	static boolean useSecondADPReleaseRateMult;
	static double nHeadsSecondADPReleaseRateMult=0;

	static String className = new String("main.MyosinMinifilament");

	static Point2D [] aOffsets;

	static Point2D [] bOffsets;

	static private void computeHeadOffsets() {

		aOffsets = new Point2D [nMyosinHeads];
		bOffsets = new Point2D [nMyosinHeads];

		//A's are the top row, and B's are the bottom row
		double spacing=0.5*myoSpacing;
		double offset=myoOffset;
		Point2D temp=new Point2D();
		for (int i=0;i<bOffsets.length;i+=2) {
			temp.x=halfLength +i*spacing;
			temp.y=halfWidth;
			//temp.y=0;
			bOffsets[i] = new Point2D(temp);
			if (nMyosinHeads!=1) {
				temp.x=halfLength + i*spacing + offset;
				temp.y=-halfWidth;
				//temp.y=-0;
				bOffsets[i+1] = new Point2D(temp);
			}
		}
		for (int i=0;i<aOffsets.length;i+=2) {
			temp.x=-halfLength - i*spacing-offset;
			temp.y=halfWidth;
			//temp.y=0;
			aOffsets[i] = new Point2D(temp);
			if (nMyosinHeads!=1) {
				temp.x=-halfLength - i*spacing;
				//temp.y=-0;
				temp.y=-halfWidth;
				aOffsets[i+1] = new Point2D(temp);
			}
		}
	}


	static public void staticInit() {
		Parameters.addParameterSetListener (
				new ParameterSetListener() {
					public void parametersChanged() {
						halfLength = 0.5*length;
						halfWidth = 0.5*width;
						collisionRadSqrd = Myosin.myoCollisionRadius*Myosin.myoCollisionRadius;
						headLength = (nMyosinHeads/2 -1)*myoSpacing;
						totalLength = length+2*headLength;
						computeHeadOffsets();
					}
				}
		);

		Parameters.addParameter(className, "numInitialMinifilaments",1,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				numInitialMinifilaments = (int)p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "nMyosinHeads",2,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				setNMyosinHeads((int) p.getDoubleValue());
			}
		}
		);
		Parameters.addParameter(className,"MyoResistiveLoad",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoResistiveLoad = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "miniFilLength",200,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				length = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "miniFilWidth",10,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				width = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoOffset",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoOffset = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "myoSpacing",5,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				myoSpacing = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "miniFilamentTurnoverRate",0.1,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				miniFilamentTurnoverRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "miniFilamentRecruitmentRate",0.01,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				miniFilamentRecruitmentRate = p.getDoubleValue();
			}
		}
		);
		Parameters.addParameter(className, "miniFilamentsTurnover",false,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				miniFilamentsTurnover = p.getBooleanValue();
			}
		}
		);

		Parameters.addParameter(className, "useSecondADPReleaseRateMult",false,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				useSecondADPReleaseRateMult = p.getBooleanValue();
			}
		}
		);
		Parameters.addParameter(className, "nHeadsSecondADPReleaseRateMult",0,
				new ParameterListener() {
			public void parameterChanged(Parameter p) {
				nHeadsSecondADPReleaseRateMult = p.getDoubleValue();
			}
		}
		);

		Sim2D.addProxy(new AgentProxy() {
			public void registerForCollisionDetection() {
				fillMFilMesh();
			}
			public void checkCollisions(double simTime) {
				myCollisionDetector.checkCollisions(simTime);
			}
			public void init() {
				myCollisionDetector.init();
			}
			public void step(double dT) {
				//  Minifilaments move and update the positions of their heads
				stepAllMinifilaments(dT);
				if(miniFilamentsTurnover) {
					doMinifilamentTurnover();
					doMinifilamentRecruitment();
				}
			}
			public void reset() {
				for(int i = 0; i < miniFilamentCt; i++) {
					theMiniFilaments[i].die();
					theMiniFilaments[i] = null;
				}
				miniFilamentCt = 0;
			}
		});
	}


}

