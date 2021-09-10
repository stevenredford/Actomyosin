package main;

/*
	A filamentous particle with Brownian diffusion and elastic collisions with other particles
*/

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Random;

/** A fascin bundle - just a renamed actin filament.  */
public class Bundle extends Actin {
	
	/** The list of all fascin bundles. */
	static Bundle [] theBundles = new Bundle [50000];
	
	/** Current number of fascin bundles. */
	static int bundleCt = 0;
			
	/** The position of this fascin bundle in the array of all fascin bundles. */
	int myBundleNumber;
	
	// CONSTRUCTORS ***** CONSTRUCTORS ***** CONSTRUCTORS ***** CONSTRUCTORS***** CONSTRUCTORS

	/** Make a new filament at (0,0), angle = 0, length = defaultMonomerLength, not static, not Pinned */
	public Bundle () {
		super(0,0,0,defaultMonomerLength,defaultMonomerLength,false,null,null);
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
	public Bundle (double initX, double initY) {
		super(initX,initY,0, defaultMonomerLength,defaultMonomerLength,false,null,null);
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
	public Bundle (double initX, double initY, double initAng, int numMons, int finalMons) {
		super(initX,initY,initAng, numMons,finalMons,false,null,null);
	}

	/**
	 * Constructor
	 *
	 * Make a new filament with the specified monomer as its barbed end.
	 *
	 * @param    newBEndMonomer      a  Monomer
	 *
	 */
	public Bundle (Monomer newB, Actin parent) {
		super(0,0);
		indexOffset = parent.indexOffset;
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
		countMonomers();
		lengthFromMonomerCount(nMonomers);
		cm.add(pEnd,length/2,uVect);
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
		evaluateProperties();
	
				
		capped = false;
		
		//System.out.println(" length="+length+"    SEVER NEW ACTIN"+this+" pEndMonomer="+pEndMonomer+"  bEndMonomer="+bEndMonomer+"   newBEndMonomer="+newBEndMonomer);
		addBundle(this);
		evaluateProperties();
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
	public Bundle (double initX, double initY, double initAng, int numMons, int finalMons, boolean isStatic, FilamentAttachment minusAttachment, FilamentAttachment plusAttachment) {
		super (initX,initY,initAng,numMons,finalMons,isStatic,minusAttachment,plusAttachment);
		uVect.set(Math.cos(initAng), Math.sin(initAng));
		uVect.uVect();
		upVect.orthogonalVector(uVect);
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
		
		for (int i=0;i<3;i++) { segLength[i] = 0; }
		createMyMonomers (numMons-2);
		length = monLength*(numMons-1);
		physicalLength = monLength*numMons;
		useCenterMassAndUVectToDetermineEnds();
		segment();
		wrapCoordinates();
		evaluateProperties();
	
		updateAllMonomerPositions();
		
		capped = false;
		
		addBundle(this);
		evaluateProperties();
	}
	
	/** Create and static variables before first simulation run. */
	static public void initDiagnostics() {
	}
	
	static public void sampleDiagnostics(double tm){
	}

	/** Reset everything in before starting next simulation run. */
	static public void resetAll() {
		for(int i = 0; i < bundleCt; i++) {
			theBundles[i].die();
			theBundles[i] = null;
		}
		bundleCt = 0;
	}
	
	public static void addBundle (Bundle newBundle) {
		newBundle.myBundleNumber = bundleCt;
		theBundles[bundleCt] = newBundle;
		bundleCt ++;
	}

	public void die() {
	}
	
		
	
	public void removeMe () {
		super.removeMe();
		int swapId = myBundleNumber;
		theBundles[swapId] = theBundles[bundleCt-1];
		theBundles[swapId].myBundleNumber = swapId;
		bundleCt --;
	}

	static public void writeHistograms() {
	}
	
	public static void printAll () {
	}

	public void setActinColor(){
		actinColor = Color.yellow;
	}

	
}
