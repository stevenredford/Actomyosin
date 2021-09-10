package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import analysis.*;
import gui.*;

public class MyoDensEvaluator extends Evaluator {
	

	

	public static int myoDensGridSize = 250;
	DataFrame theFrame;
	DataMesh theMesh;
;	RemoteImageWriter remoteImageWriter;
	MyosinHeadGrabber theGrabber;

	
	Point2D cenOfDisc;
	public int radialMyoDensArraySize = 100;
	public double [] radialMyoDens = new double [radialMyoDensArraySize];
	public double radialMyoDensDistInc = Sim2D.xDimension/radialMyoDensArraySize;
	

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
	}
	
	public void reset()  {
		super.reset();
		if(theMesh == null) {
			theGrabber = new MyosinHeadGrabber();
			theMesh = new DataMesh (myoDensGridSize, DataMesh.SCALAR, "MyoDensity", theGrabber, 15);
		}
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "MonDens.dat";
	}
	
	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("MyoDensEvaluator.mkDataFile(): error creating File"); }
	}
	
	public void resetGraphics(boolean remote) {
		
		try
        {
			if(!remote) {
				theFrame = new DataFrame(theMesh,"Myosin Density");
				theFrame.updateGraphics();
			}
		} catch (Exception e)
        {
            System.err.println(e);
            e.printStackTrace();
        }
	}

	public void showGraphics() {
		theFrame.updateGraphics();
	}
	
	public void pngFromBufferedImage(String path, String name) {
		PNGer.pngFromBufferedImage(theMesh.getPanel().getImage(),new String(path + File.separator + name));
	}
		
	public void remoteWriteMovieFrameToPng(String path, String name) {
		if(remoteImageWriter == null) {
			remoteImageWriter = new RemoteImageWriter(WorldFrame.width,WorldFrame.height,theMesh);
			remoteImageWriter.initialize();
		}
		remoteImageWriter.draw();
 		PNGer.pngFromBufferedImage(remoteImageWriter.getImage(),new String(path + File.separator + name));
	}
		

	/** This Evaluator writes data files. */
	public boolean hasData()  {
		return true;
	}

	/** This Evaluator writes movie frames. */
	public boolean hasMovie()  {
		return true;
	}
	
	public void evaluate(double tm)  {
		writeMyoDensData();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}
	
	public double [] mapOfRadialMyoDens() {
		clearRadialMyoDens();
		for(int i=0;i<Myosin.myosinCt; i++) {
			mapToRadialMyoDens(Myosin.theMyosins[i]);
		}
		return radialMyoDens;
	}
	
	public void clearRadialMyoDens () {
		for (int i=0;i<radialMyoDensArraySize;i++) {
			radialMyoDens[i] = 0;
		}
		cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
		radialMyoDensDistInc = Sim2D.xDimension/radialMyoDensArraySize;
	}
		
	public void mapToRadialMyoDens (Myosin m) {
		Point2D cm = m.getLocation();
		double dist = Point2D.getDistance(cm,cenOfDisc);
		int rLoc = (int)(dist/radialMyoDensDistInc);
		radialMyoDens[rLoc]++;
	}
	

	
	
	public void writeMyoDensData () {
		double [] data = mapOfRadialMyoDens();
		for(int i = 0; i < data.length; i++) {
			dataPW.print(data[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.flush();
	}
				
	
	public void cleanUp()  {}
				

}
