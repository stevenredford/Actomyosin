package iterators;

import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;
import heatMap.*;
import parameters.*;

public class GlobalStressEvaluator extends Evaluator {
	

	static String
		
		RADIAL = new String("Radial"),
		CIRCUMFERENTIAL = new String("Circumferential"),
		HORIZONTAL = new String("Horizontal"),
		VERTICAL = new String("Vertical");
	
	String arrayOrientation;

	double totalStress = 0,totalMyoStress = 0,totalLinkerStress = 0,totalActinStress = 0;
	double posMyoStress = 0, posLinkerStress = 0, posActinStress = 0;
	double negMyoStress = 0, negActinStress = 0, negLinkerStress = 0;
	double maxMyoStress, maxActinStress, maxLinkerStress;
	
	public static int stressGridSize = 200;
	public double [][] stresses = new double [stressGridSize][stressGridSize];
	int [][] stressCt = new int [stressGridSize][stressGridSize];
	double stressDistInc;
	
	HeatMap remotePanel;
		
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equalsIgnoreCase("orientation")) {
			arrayOrientation = in.nextString();
			if(!(arrayOrientation.equals(VERTICAL) || arrayOrientation.equals(HORIZONTAL) || arrayOrientation.equals(RADIAL) || arrayOrientation.equals(CIRCUMFERENTIAL))) {
				throw new Exception("GlobalStressEvaluator.loadParameter(): got bad value for orientation = " + arrayOrientation);
			}
		}
		else {
			super.loadParameter(tag, in);
		}
	}
	
	public void init(String path, String name) throws Exception  {
		super.init(path,name);
		Parameters.setParameter("monitorStress",true);
	}
	
	public void reset()  {
		super.reset();
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "Stress.dat";
	}
			
	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("GlobalStressEvaluator.mkDataFile(): error creating File"); }
	}

	public void remoteWriteMovieFrameToPng(String path, String name) {
		mapAllToStress();
		String stressFileName = path + File.separator + name;
    	remotePanel.updateData(stresses, true,true);
		PNGer.pngFromBufferedImage(remotePanel.bufferedImage,stressFileName);
	}

	public void showGraphics() {
		mapAllToStress();
		DotAngMap.updateRadialPolarityMap(stresses,false);
	}
	
	public void pngFromBufferedImage(String path, String name) {
		PNGer.pngFromBufferedImage(DotAngMap.panel.bufferedImage,new String(path + File.separator + name));
	}

	public void resetGraphics(boolean remote) {
		try
        {
			mapAllToStress();
			if(remote) {
				remotePanel = new HeatMap(stresses, true, Gradient.GRADIENT_BLACK_TO_WHITE, true);
			}
			else {
				DotAngMap.createAndShowGUI(stresses,remote);
			}
        }
        catch (Exception e)
        {
            System.err.println(e);
            e.printStackTrace();
        }
	}

	
	public void mapAllToStress () {
		clearStresses();
		for (int i=0;i<Monomer.monomerCt;i++) {
			mapToStress(Monomer.theMonomers[i]);
		}
		averageStressValues();
	}

	public void clearStresses () {
		for (int i=0;i<stressGridSize;i++) {
			for (int j=0;j<stressGridSize;j++) {
				stressCt[i][j] = 0;
				stresses[i][j] = 0;
			}
		}
		stressDistInc = Sim2D.xDimension/stressGridSize;
	}
	
	public void averageStressValues () {
		for (int i=0;i<stressGridSize;i++) {
			for (int j=0;j<stressGridSize;j++) {
				if (stressCt[i][j]==0) {
					stresses[i][j] = 0;
				} else {
					stresses[i][j] /= stressCt[i][j];
				}
			}
		}
	}
		
	public void mapToStress (Monomer m) {
		Point2D cm = m.getLocation();
		int xLoc = (int)(cm.x/stressDistInc);
		int yLoc = (int)(cm.y/stressDistInc);
		stressCt[xLoc][yLoc]++;
		stresses[xLoc][yLoc]+=m.avgMonStress;
	}
	
	private void measureGlobalStress() {
		totalMyoStress = totalLinkerStress = totalActinStress = 0;
		posMyoStress = posLinkerStress = posActinStress = 0;
		negMyoStress = negActinStress = negLinkerStress = 0;
		maxMyoStress = maxActinStress = maxLinkerStress = 0;
		
		getForcesForAllActins();
		getForcesForAllMyos();
		getForcesForAllLinkers();
		
		double L = Sim2D.yDimension;
		posMyoStress /= L;
		posLinkerStress /= L;
		posActinStress /= L;
		negMyoStress /= L;
		negLinkerStress /= L;
		negActinStress /= L;
		totalMyoStress /= L;
		totalLinkerStress /= L;
		totalActinStress /= L;
		totalStress = totalMyoStress + totalLinkerStress + totalActinStress;
	}

	private void getForcesForAllActins() {
		Monomer curM;
		double stress;
		double orthoComponent;
		for(int i = 0; i < Monomer.monomerCt; i++) {
			curM = Monomer.theMonomers[i];
			orthoComponent = Math.abs(curM.myFilament.uVect.x);
			stress = curM.avgMonStress*orthoComponent*Actin.monLength;
			if(stress > 0)
				posActinStress += stress;
			else negActinStress += stress;
			if(stress > maxActinStress) maxActinStress = stress;
			totalActinStress += stress;
		}
	}

	private void getForcesForAllMyos() {
		MyosinMiniFilament curM;
		double  stress;
		double orthoComponent;
		int cInd = MyosinMiniFilament.nMyosinHeads/2 - 1;
		for(int i = 0; i < MyosinMiniFilament.miniFilamentCt; i++) {
			curM = MyosinMiniFilament.theMiniFilaments[i];
			orthoComponent = Math.abs(curM.uVect.x);
			for(int j = 0; j < MyosinMiniFilament.nMyosinHeads-1; j++) {
				stress = curM.getInternalStress(j)*orthoComponent;
				if(j == cInd) {
					stress *= MyosinMiniFilament.length;
				} else stress *= MyosinMiniFilament.myoSpacing;
				if(stress > 0)
					posMyoStress += stress;
				else negMyoStress += stress;
				if(curM.getInternalStress(j) > maxMyoStress) maxMyoStress = curM.getInternalStress(j);
				totalMyoStress += stress;
			}
		}
		
	}
	
	private void getForcesForAllLinkers() {
		Crosslinker curL;
		double  stress;
		double orthoComponent;
		for(int i = 0; i < Crosslinker.crosslinkerCt; i++) {
			curL = Crosslinker.theCrosslinkers[i];
			if(curL.isBound()) {
				orthoComponent = Math.abs(curL.uVect.x);
				stress = curL.forceAv*orthoComponent*curL.getLength();
				if(stress > 0)
					posLinkerStress += stress;
				else negLinkerStress += stress;
				if(curL.forceAv > maxLinkerStress) maxLinkerStress = curL.forceAv;
				totalLinkerStress += stress;
			}
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
		measureGlobalStress();
		writeStressData();
	}
	
	public boolean decide(double tm) {
		return super.decide(tm);
	}
	
	public void writeDataFileHeader()  {
		dataPW.print("posMyoStress" + "\t" + "negMyoStress" + "\t" + "totalMyoStress" + "\t");
		dataPW.print("posActinStress" + "\t" + "negActinStress" + "\t" + "totalActinStress" + "\t");
		dataPW.print("posLinkerStress" + "\t" + "negLinkerStress" + "\t" + "totalLinkerStress" + "\t");
		dataPW.print("totalStress" + "\n");
	}
		
	private void writeStressData()  {

		dataPW.print(posMyoStress + "\t" + negMyoStress + "\t" + totalMyoStress + "\t");
		dataPW.print(posActinStress + "\t" + negActinStress + "\t" + totalActinStress + "\t");
		dataPW.print(posLinkerStress + "\t" + negLinkerStress + "\t" + totalLinkerStress + "\t");
		dataPW.print(totalStress + "\n");
		dataPW.flush();
	}

	
	
	public void cleanUp()  {}
				

}
