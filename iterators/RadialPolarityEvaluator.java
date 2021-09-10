package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import heatMap.*;

public class RadialPolarityEvaluator extends Evaluator {
	

	Point2D rVec = new Point2D();
	Point2D cm = new Point2D();
	double dist;
	
	public int radialArraySize = 25;
	public double radialInc = 0.5*Sim2D.xDimension/(radialArraySize*Math.sqrt(2.));
	Point2D cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
	
	public double [] radialPolarityArray = new double [radialArraySize];
	public int [] radialArrayCts = new int[radialArraySize];
	
	public static int radialPolarityGridSize = 200;
	public double [][] radialPolarityGrid = new double [radialPolarityGridSize][radialPolarityGridSize];
	int [][] gridCt = new int [radialPolarityGridSize][radialPolarityGridSize];
	double angleDistInc;
			
	HeatMap remotePanel;
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
	}
	

	public void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "Radialpolarities.dat";
	}
	
	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("RadialPolarityEvaluator.mkDataFile(): error creating File"); }
	}


	public void clearRadialPolarities () {
		for (int i=0;i<radialArraySize;i++) {
			radialPolarityArray[i] = 0;
			radialArrayCts[i] = 0;
		}
	}
		
	public void mapToRadialArray (Monomer m) {
		Actin f = m.getFilament();
		if (f == null) return;
		Point2D cm = m.getLocation();
		dist = Point2D.getDistance(cm,cenOfDisc);
		int rLoc = (int)(dist/radialInc);
		if(rLoc >= radialArraySize) rLoc = radialArraySize-1;
		rVec.getVector(cm,cenOfDisc);
		rVec.uVect();
		double rDot = Point2D.dot(f.uVect, rVec);
		radialPolarityArray[rLoc]+=rDot;
		radialArrayCts[rLoc]++;
	}
	
	public void mapAllToRadialArray() {
		clearRadialPolarities();
		for (int i=0;i<Monomer.monomerCt;i++) {
			mapToRadialArray(Monomer.theMonomers[i]);
		}
		normalizeRadialPolarities();
	}
	
	public void averageRadialPolarities () {
		for (int i=0;i<radialArraySize;i++) {
			if (radialArrayCts[i]==0) {
				radialPolarityArray[i] = 0;
			} else {
				radialPolarityArray[i] /= radialArrayCts[i];
			}
		}
	}
	
	/** scale radial polarity values by the max bin count */
	public void normalizeRadialPolarities () {
		int maxCt = 0;
		for (int i=0;i<radialArraySize;i++) {
			if (radialArrayCts[i] > maxCt) maxCt = radialArrayCts[i];
		}
		for (int i=0;i<radialArraySize;i++) {
			radialPolarityArray[i] /= maxCt;
		}
	}
	
	public void mapAllToGrid () {
		clearRadialPolarityGrid();
		for (int i=0;i<Monomer.monomerCt;i++) {
			mapToGrid(Monomer.theMonomers[i]);
		}
		averageGrid();
	}

	public void clearRadialPolarityGrid () {
		for (int i=0;i<radialPolarityGridSize;i++) {
			for (int j=0;j<radialPolarityGridSize;j++) {
				gridCt[i][j] = 0;
				radialPolarityGrid[i][j] = 0;
			}
		}
		cenOfDisc = new Point2D(Sim2D.xDimension/2,Sim2D.yDimension/2);
		angleDistInc = Sim2D.xDimension/radialPolarityGridSize;
	}
	
	public void averageGrid () {
		for (int i=0;i<radialPolarityGridSize;i++) {
			for (int j=0;j<radialPolarityGridSize;j++) {
				if (gridCt[i][j]==0) {
					radialPolarityGrid[i][j] = 0;
				} else {
					radialPolarityGrid[i][j] /= gridCt[i][j];
				}
			}
		}
	}
		
	public void mapToGrid (Monomer m) {
		Actin f = m.getFilament();
		if (f == null) return;
		Point2D cm = m.getLocation();
		rVec.getVector(cm,cenOfDisc);
		rVec.uVect();
		double rDot = Point2D.dot(f.uVect, rVec);
		int xLoc = (int)(cm.x/angleDistInc);
		int yLoc = (int)(cm.y/angleDistInc);
		if(xLoc == radialArraySize) xLoc = radialArraySize-1;
		if(yLoc == radialArraySize) yLoc = radialArraySize-1;
		gridCt[xLoc][yLoc]++;
		radialPolarityGrid[xLoc][yLoc]+=rDot;
	}

	/** This Evaluator writes data files. */
	public boolean hasData()  {
		return true;
	}

	/** This Evaluator writes movie frames? */
	public boolean hasMovie()  {
		return true;
	}
	
	public void evaluate(double tm)  {
		writeRadialPolarityData();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}
	
	public String getHeaderString() {
		return new String("\t" + name + "\t" + "Max" + "\t" + "WhereMax" + "\t" + "Min" + "\t" + "WhereMin"  + "\t" + "Mean");
	}
	
	public String getDataString() {
		double max = -1e20, min = 1e20, mean = 0;
		double whereMin = -1, whereMax = -1;
		mapAllToRadialArray();
		for(int i = 0; i < radialArraySize; i++) {
			if(radialPolarityArray[i] > max) {
				max = radialPolarityArray[i];
				whereMax = i;
			}
			if(radialPolarityArray[i] < min) {
				min = radialPolarityArray[i];
				whereMin = i;
			}
			mean += radialPolarityArray[i];
		}
		mean /= radialArraySize;
		whereMax/=radialArraySize;
		whereMin/=radialArraySize;
		return new String("\t" + name + "\t" + max + "\t" + whereMax + "\t" + min  + "\t" + whereMin + "\t" + mean);
	}

	public void remoteWriteMovieFrameToPng(String path, String name) {
		mapAllToGrid();
		String radialGridFileName = path + File.separator + name;
    	remotePanel.updateData(radialPolarityGrid, true,true);
		PNGer.pngFromBufferedImage(remotePanel.bufferedImage,radialGridFileName);
	}

	public void showGraphics() {
		mapAllToGrid();
		DotAngMap.updateRadialPolarityMap(radialPolarityGrid,false);
	}
	
	public void pngFromBufferedImage(String path, String name) {
		String radialGridFileName = path + File.separator + name;
		PNGer.pngFromBufferedImage(DotAngMap.panel.bufferedImage,radialGridFileName);
	}

	public void resetGraphics(boolean remote) {
		try
        {
			mapAllToGrid();
			if(remote) {
				remotePanel = new HeatMap(radialPolarityGrid, true, Gradient.GRADIENT_BLACK_TO_WHITE, true);
			}
			else {
				DotAngMap.createAndShowGUI(radialPolarityGrid,remote);
			}
        }
        catch (Exception e)
        {
            System.err.println(e);
            e.printStackTrace();
        }
	}
	public void writeRadialPolarityData () {
		mapAllToRadialArray();
		dataPW.print("radialPolarity\t");
		for(int i = 0; i < radialPolarityArray.length; i++) {
			dataPW.print(radialPolarityArray[i] + "\t");
		}
		dataPW.print("\n");
		dataPW.flush();
	}
	
	public void cleanUp()  {}
				

}
