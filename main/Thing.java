package main;

/*
/ Thing.... the superclass for all moving objects in this demonstration
*/

import java.awt.*;

import util.*;

public class Thing extends Object {
	
	/** List of all things. */
	static Thing [] theThings = new Thing [1000000];
	
	/** Current number of all things. */
	static int thingCt = 0;
	
	/** Connectivity tracking */
	public boolean isConnected = false;
	public int connectedLevel = 1000000;	// crazy-large number by default
	public double clusterTime = 0;
	public int clusterId = 0;
	
	/** Tendril tracking */
	public double tendrilTime = 0;
	public int [] tendrilIDs = new int [100];
	public int tendrilCt = 0;
	
	/** List of dead things scheduled for removal. */
	static Thing [] deadThings = new Thing [1000000];
	
	/** Number of dead things scheduled for removal. */
	static int deadCt = 0;
	
	/** Position of this thing in the array of all things. */
	int myThingNumber;
	
	/** if true this Thing will be eliminated */
	public boolean removeMe = false;
	
	/** Center of mass of this thing. */
	public Point2D cm = new Point2D();
	
	/** maximum x position this Thing can occupy */
	static double maxX = Sim2D.xDimension;
	
	/** maximum y position this Thing can occupy */
	static double maxY = Sim2D.yDimension;
		
	/** An object passed to and from line-line and line-point intersect tests */
	public class RetObj {
		
		Point2D conPt1, conPt2, ray1, ray2, ray3, ray4;
		double conDist = 0;
		boolean collision = false;
		
		public RetObj () {
			conPt1 = new Point2D();
			conPt2 = new Point2D();
			ray1 = new Point2D();
			ray2 = new Point2D();
			ray3 = new Point2D();
		}
		
		public void reset () {
			collision = false;
		}
		
	}
	
	public void reset () {
		isConnected = false;
		connectedLevel = 1000000;
	}
	
	/** Reset everything in before startiung next simulation run. */
	static public void resetAll() {
		for(int i = 0;i < thingCt; i++) {
			theThings[i] = null;
		}
		thingCt = 0;
	}
		
			
	
	/** Make a thing with coordinates (initX,intiY). */
	public Thing (double initX, double initY) {
		cm.set(initX,initY);
		addThing(this);
	}
	
	/** Make a thing at initPos */
	public Thing (Point2D initPos) {
		cm.set(initPos);
		addThing(this);
	}
	
	/** Override this to move this object each time-step */
	public void step (double dT) {
		
	}
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		// put the code here to draw the object on "g"
	}
	
	/** Override this to add a force to this thing. */
	public void addForce(Point2D f){
	}
	
	/** Add this thing to the list of all things. */
	public static void addThing (Thing newThing) {
		newThing.myThingNumber = thingCt;
		theThings[thingCt] = newThing;
		thingCt++;
	}
	/** Add a thing to the list of things scheduled for removal (= death). */
	public static void addDeadThing (Thing deadThing) {
		deadThings[deadCt] = deadThing;
		deadCt++;
	}
	
	/** Schedule this thing for removal. */
	public void setForRemove () {
		removeMe = true;
		addDeadThing(this);
	}
	
	/** Remove this thing */
	public void removeMe () {
		//System.out.println("In removeMe: " + Thread.currentThread());
		int swapId = myThingNumber;
		theThings[swapId] = theThings[thingCt-1];
		theThings[swapId].myThingNumber = swapId;
		thingCt --;
	}
		
	/** Remove all things in the list of dead things. */
	public static void removeArrayedDeadThings () {
		for (int i=0;i<deadCt;i++) {
			if(deadThings[i] != null) {
				deadThings[i].removeMe();
			}
		}
		deadCt = 0;
	}
	
	/** Return a random point in the simulation arena. */
	public static Point2D randomPointInArena () {
		double initX = Math.random()*Sim2D.xDimension;
		double initY = Math.random()*Sim2D.yDimension;
		return new Point2D(initX,initY);
	}
	
	
	/*8 Override this to return this things position. */
	public Point2D getLocation(){
		return null;
	}
	
}

