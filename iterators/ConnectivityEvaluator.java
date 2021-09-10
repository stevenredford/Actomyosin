package iterators;


import java.util.*;
import java.awt.Color;
import java.io.*;

import collision.CollisionDetector;

import main.*;
import util.*;
import io.*;
import parameters.*;
import gui.*;

public class ConnectivityEvaluator extends Evaluator {
	
	double X1,Y1,X2,Y2;
	Point2D p1 = new Point2D();
	Point2D p2 = new Point2D();
	Point2D segUVect = new Point2D();
	Point2D forceSum = new Point2D();
	double forceOrtho = 0;
	int maxAs = 1000;
	Actin [] crossingAs = new Actin[maxAs];
	double [] crossingAsArc = new double[maxAs];
	int crossingAsCt = 0;
	Color forceLineColor = Color.orange;
	double maxFraction;
	double whereMax;
	
;	RemoteImageWriter remoteImageWriter;
	WorldFrame theFrame;
	
	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus clusterSizeHisto;
	
	public ConnectivityEvaluator() {}

	public void init(String path, String name) throws Exception  {
		super.init(path,name);
	}
	
	/** Create ClusterSize histo at this point when we know how many filaments so we can set appropriate binCt. */
	public void doFirstEval() {
		maxFraction = -1e20;
		whereMax = -1;
		clusterSizeHisto = new HistogramPlus(100,0,Actin.actinCt,dataPath,"clusterSizes",null,true,false,false,false);
	}
	
	public void evaluate (double tm) {
		updateClusterSizeHisto();
		double val = clusterSizeHisto.getMaxValue()/Actin.actinCt;
		if(val > maxFraction) {
			maxFraction = val;
			whereMax = tm;
		}
	}
	
	public void writeDataFileHeader () {
		clusterSizeHisto.printHeader(dataPW);
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "ConnectivityEval.dat";
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("ConnectivityEvaluator.mkDataFile(): error creating File"); }
	}

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if (false) {
			
		}
		else super.loadParameter(tag,in);
	}

	public void resetGraphics(boolean remote) {
		
		try
        {
			if(!remote) {
				theFrame = new WorldFrame(false);
				showGraphics();
			}
		} catch (Exception e)
        {
            System.err.println(e);
            e.printStackTrace();
        }
	}

	public void showGraphics() {
		boolean state = Sim2D.paintConnected;
		Parameters.setParameter("paintConnected",true);
		Actin.traceAllClusters();
		Cluster.setNLargestClusters(Sim2D.numClustersToPaint);
		theFrame.showAll();
		Parameters.setParameter("paintConnected",state);
	}
	
	public void pngFromBufferedImage(String path, String name) {
		PNGer.pngFromBufferedImage(theFrame.getImage(),new String(path + File.separator + name));
	}
		
	public void remoteWriteMovieFrameToPng(String path, String name) {
		if(remoteImageWriter == null) {
			remoteImageWriter = new RemoteImageWriter(WorldFrame.width,WorldFrame.height);
			remoteImageWriter.initialize();
		}
		boolean state = Sim2D.paintConnected;
		Parameters.setParameter("paintConnected",true);
		Actin.traceAllClusters();
		Cluster.setNLargestClusters(Sim2D.numClustersToPaint);
		remoteImageWriter.draw();
 		PNGer.pngFromBufferedImage(remoteImageWriter.getImage(),new String(path + File.separator + name));
		Parameters.setParameter("paintConnected",state);
	}
		
		
	void updateClusterSizeHisto() {
		clusterSizeHisto.clearBins();
		Actin.traceAllClusters();
		for (int i=0;i<Cluster.clusterCt;i++) {
			clusterSizeHisto.addValue(Cluster.theClusters[i].elementCt);
			//System.out.println ("adding cluster with " + Cluster.theClusters[i].elementCt + " elements to histogram");
		}
		clusterSizeHisto.writeToFile(dataPW);
	}
	
	public double getFracElementsInLargestCluster () {
		double elementsInLargestCluster = Cluster.getLargestClusterSize();
		double totalElements = Actin.actinCt + Crosslinker.crosslinkerCt + MyosinMiniFilament.miniFilamentCt;
		return elementsInLargestCluster/totalElements;
	}

	
	public boolean stop(double tm)  {
		if(super.stop(tm)) return true;
		
		return false;
	}
	
	public String getHeaderString() {
		return new String("\t" + name + "\t" + "MaxClusterFraction" + "\t" + "WhenMax");
	}
	
	public String getDataString() {
		return new String("\t" + name + "\t" + maxFraction + "\t" + whereMax);
	}
	
	
	public boolean hasData() {
		return true;
	}
	
	public boolean hasMovie() {
		return true;
	}

	private void storeData(int i){
		
	}
	
}
