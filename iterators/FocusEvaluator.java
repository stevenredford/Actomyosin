package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import analysis.*;
import gui.*;

public class FocusEvaluator extends Evaluator {
	

	

	public double avgFocusArea = 0;
	public int numFoci = 0;
	public double focusSizeThreshold = 0;
	public static int myoDensGridSize = 250;
	DataFrame theFrame;
	DataMesh theMesh;
;	RemoteImageWriter remoteImageWriter;
	BlobFinder theFinder;

	

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equals("FocusSizeThreshold")) {
			focusSizeThreshold = in.nextDouble();
		}
		else super.loadParameter(tag, in);
	}
	
	public void reset()  {
		super.reset();
		if(theMesh == null) {
			theFinder = new BlobFinder();
			theFinder.setFocusSizeThreshold(focusSizeThreshold);
			theMesh = new DataMesh (myoDensGridSize, DataMesh.SCALAR, "MyoDensity", theFinder, 15);
		}
	}
	
	private void calculateFocusStatistics() {
		theFinder.grabData();
		numFoci = theFinder.getNumBlobs();
		avgFocusArea = theFinder.getAverageBlobArea();
	}
		
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "MyosinFoci.dat";
	}
	
	
	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("FocusEvaluator.mkDataFile(): error creating File"); }
	}

	public void resetGraphics(boolean remote) {
		
		try
        {
			if(!remote) {
				theFrame = new DataFrame(theMesh,"Myosin Blobs");
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

	public String getHeaderString() {
		return new String("\t" + "NumberOfFoci" + "\t" + "AverageFocusArea");
	}
	
	public String getDataString() {
		calculateFocusStatistics();
		return new String("\t" + numFoci + "\t" + avgFocusArea);
	}

	/** This Evaluator writes movie frames. */
	public boolean hasMovie()  {
		return true;
	}
	
	public void evaluate(double tm)  {
		writeFocusData();
	}
	
	public void writeDataFileHeader () {
		dataPW.println(new String("Simulation time" + "\t" + "numFoci" + "\t" + "avgFocusArea"));
		dataPW.flush();
	}

	public void writeFocusData () {
		calculateFocusStatistics();
		dataPW.println(Sim2D.simulationTime + "\t" + numFoci + "\t" + avgFocusArea);
		dataPW.flush();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}
				
	
	public void cleanUp()  {}
				

}
