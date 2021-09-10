package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import heatMap.*;

public class MonDensEvaluator extends Evaluator {
	
	
	public static int monDensGridSize = 200;
	public double [][] monDens = new double [monDensGridSize][monDensGridSize];
	int [][] monDensCt = new int [monDensGridSize][monDensGridSize];
	
	public int radialMonDensArraySize = 100;
	public double [] radialMonDens = new double [radialMonDensArraySize];
	public double radialMonDensDistInc = Sim2D.xDimension/radialMonDensArraySize;
	public double monDensDistInc = 0.5*Sim2D.xDimension/(radialMonDensArraySize*Math.sqrt(2.));
	Point2D cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
	

	HeatMap remotePanel;

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
	}
	
	public void reset()  {
		super.reset();
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "MonDens.dat";
	}
	
		
	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("MonDensEvaluator.mkDataFile(): error creating File"); }
	}

	public void remoteWriteMovieFrameToPng(String path, String name) {
		mapAllToMonDensHeatMap();
		String monDensFileName = path + File.separator + name;
    	remotePanel.updateData(monDens, true,true);
		PNGer.pngFromBufferedImage(remotePanel.bufferedImage,monDensFileName);
	}

	public void showGraphics() {
		mapAllToMonDensHeatMap();
		MonDensMap.updateMonDensMap(monDens,false);
	}
	
	public void pngFromBufferedImage(String path, String name) {
		String monDensFileName = path + File.separator + name;
		PNGer.pngFromBufferedImage(DotAngMap.panel.bufferedImage,monDensFileName);
	}


	public void resetGraphics(boolean remote) {
		try
        {
			mapAllToMonDensHeatMap();
			if(remote) {
				remotePanel = new HeatMap(monDens, true, Gradient.GRADIENT_BLACK_TO_WHITE, true);
			}
			else {
				MonDensMap.createAndShowGUI(monDens,remote);
			}
        }
        catch (Exception e)
        {
            System.err.println(e);
            e.printStackTrace();
        }
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
		writeMonDensData();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}

	public double [] mapOfRadialMonDens() {
		clearRadialMonDens();
		for (int i=0;i<Monomer.monomerCt;i++) {
			mapToRadialMonDens(Monomer.theMonomers[i]);
		}
		return radialMonDens;
	}
	
	public void clearRadialMonDens () {
		for (int i=0;i<radialMonDensArraySize;i++) {
			radialMonDens[i] = 0;
		}
		monDensDistInc = 0.5*Sim2D.xDimension/(radialMonDensArraySize*Math.sqrt(2.));
		cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
	}
		
	public void mapToRadialMonDens (Monomer m) {
		if (m.getFilament() == null) return;
		Point2D cm = m.getLocation();
		double dist = Point2D.getDistance(cm,cenOfDisc);
		int rLoc = (int)(dist/radialMonDensDistInc);
		radialMonDens[rLoc]++;
	}
	
	public void mapAllToMonDensHeatMap() {
		clearMonDens();
		for (int i=0;i<Monomer.monomerCt;i++) {
			mapToMonDensHeatMap(Monomer.theMonomers[i]);
		}
	}

	public void clearMonDens () {
		for (int i=0;i<monDensGridSize;i++) {
			for (int j=0;j<monDensGridSize;j++) {
				monDensCt[i][j] = 0;
				monDens[i][j] = 0;
			}
		}
		cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
		monDensDistInc = Sim2D.xDimension/monDensGridSize;
	}
		
	public void mapToMonDensHeatMap (Monomer m) {
		if (m.getFilament() == null) return;
		Point2D cm = m.getLocation();
		int xLoc = (int)(cm.x/monDensDistInc);
		int yLoc = (int)(cm.y/monDensDistInc);
		monDensCt[xLoc][yLoc]++;
		monDens[xLoc][yLoc]++;
	}
	
	
	public void writeMonDensData () {
		double [] data = mapOfRadialMonDens();
		for(int i = 0; i < data.length; i++) {
			dataPW.print(data[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.flush();
	}
			
	
	public void cleanUp()  {}
				

}
