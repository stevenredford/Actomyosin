package main;

public class ForceTendril {
	static int maxTendrils = 10000;
	public static ForceTendril [] theTendrils = new ForceTendril[maxTendrils];
	static ForceTendril [] sortedTendrils = new ForceTendril[maxTendrils];
	public static int tendrilCt = 0;
	
	static boolean followMyos = true;
	static boolean followXLs = true;
	
	static int maxTendrilElements = 10000;
	public Thing [] tendrilElements = new Thing[maxTendrilElements];
	public int elementCt = 0;
	public int tendrilID = 0;
	public boolean sorted = false;
	
	public ForceTendril() {
		
	}
	
	public void addElement (Thing newE) {
		tendrilElements[elementCt] = newE;
		newE.clusterId = tendrilID;
		newE.clusterTime = Sim2D.simulationTime;
		elementCt++;
	}
	
	public void zero () {
		elementCt = 0;
	}
	
	public static void startTendril (Myosin startMyo) {
		if (tendrilCt >= maxTendrils) { System.out.println("** Max Tendrils Exceeded **"); return; }
		if (theTendrils[tendrilCt] == null) { 
			theTendrils[tendrilCt] = new ForceTendril(); 
			theTendrils[tendrilCt].tendrilID = tendrilCt;
		}
		theTendrils[tendrilCt].zero();
		startMyo.traceTendril(theTendrils[tendrilCt]);
		tendrilCt++;
	}
	
	public static void resetClusters () {
		tendrilCt = 0;
	}
	
	public static void sortClustersBySize () {
		
	}
	
	public static void setNLargestClusters (int n) {
		if (tendrilCt == 0) return;
		for (int i=0;i<tendrilCt;i++) { theTendrils[i].sorted = false; }  // reset all flags
		
		int numSorted = 0;
		while (numSorted < n) {
			sortedTendrils[numSorted] = getNextLargestNotSorted (numSorted);
			sortedTendrils[numSorted].sorted = true;
			numSorted++;
		}			
	}
	
	public static ForceTendril getNextLargestNotSorted (int numSorted) {
		int largestSize = 0;
		ForceTendril largestTendril = null;
		ForceTendril curTendril;
		for (int i=0;i<tendrilCt;i++) {
			curTendril = theTendrils[i];
			if (!curTendril.sorted && curTendril.elementCt >= largestSize) { 
				largestSize = curTendril.elementCt; 
				largestTendril = curTendril;
			}
		}
		return largestTendril;
	} 
	
	public static int getLargestClusterSize() {
		int largestSize = 0;
		for (int i=0;i<tendrilCt;i++) {
			if (theTendrils[i].elementCt > largestSize) { largestSize = theTendrils[i].elementCt; }
		}
		return largestSize;
	}
	
	public static void reportClusterStats () {
		System.out.println("Number of clusters: " + tendrilCt);
		if (tendrilCt > 0) {
			int n=10;
			setNLargestClusters(n);
			for (int i=0;i<n;i++) {
				System.out.println("cluster#" + i + " is " + sortedTendrils[i] + " and has " + sortedTendrils[i].elementCt + " elements");
			}
		}
	}
}
