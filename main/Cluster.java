package main;

public class Cluster {
	static int maxClusters = 10000;
	public static Cluster [] theClusters = new Cluster[maxClusters];
	static Cluster [] sortedClusters = new Cluster[maxClusters];
	static int bigClusterCt = 0;
	public static int clusterCt = 0;
	
	static boolean followMyos = true;
	static boolean followXLs = true;
	
	static int maxClusterElements = 10000;
	public Thing [] clusterElements = new Thing[maxClusterElements];
	public int elementCt = 0;
	public int clusterId = 0;
	public boolean sorted = false;
	
	public Cluster () {
		
	}
	
	public void addElement (Thing newE, boolean incrementCt) {
		clusterElements[elementCt] = newE;
		newE.clusterId = clusterId;
		newE.clusterTime = Sim2D.simulationTime;
		if(incrementCt) elementCt++;
	}
	
	public void zero () {
		elementCt = 0;
	}
	
	public static void startCluster (Actin newE) {
		if (clusterCt >= maxClusters) { System.out.println("** Max Clusters Exceeded **"); return; }
		if (theClusters[clusterCt] == null) {
			theClusters[clusterCt] = new Cluster();
			theClusters[clusterCt].clusterId = clusterCt;
		}
		theClusters[clusterCt].zero();
		newE.setConnected(theClusters[clusterCt]);
		clusterCt++;
	}
	
	public static void resetClusters () {
		clusterCt = 0;
	}
	
	public static void sortClustersBySize () {
		
	}
	
	public static void setNLargestClusters (int n) {
		bigClusterCt = 0;
		if (clusterCt == 0) return;
		for (int i=0;i<clusterCt;i++) { theClusters[i].sorted = false; }  // reset all flags
		
		int numSorted = 0;
		while (numSorted < n) {
			sortedClusters[numSorted] = getNextLargestNotSorted (numSorted);
			if(sortedClusters[numSorted] == null) return;
			sortedClusters[numSorted].sorted = true;
			numSorted++;
			bigClusterCt++;
		}
	}
	
	public static Cluster getNextLargestNotSorted (int numSorted) {
		int largestSize = 0;
		Cluster largestCluster = null;
		Cluster curCluster;
		for (int i=0;i<clusterCt;i++) {
			curCluster = theClusters[i];
			if (!curCluster.sorted && curCluster.elementCt >= largestSize) {
				largestSize = curCluster.elementCt;
				largestCluster = curCluster;
			}
		}
		return largestCluster;
	}
	
	public static int getLargestClusterSize() {
		int largestSize = 0;
		for (int i=0;i<clusterCt;i++) {
			if (theClusters[i].elementCt > largestSize) { largestSize = theClusters[i].elementCt; }
		}
		return largestSize;
	}
	
	public static void reportClusterStats () {
		System.out.println("Number of clusters: " + clusterCt);
		if (clusterCt > 0) {
			int n=10;
			setNLargestClusters(n);
			for (int i=0;i<bigClusterCt;i++) {
				System.out.println("cluster#" + i + " is " + sortedClusters[i] + " and has " + sortedClusters[i].elementCt + " elements");
			}
		}
	}
}
