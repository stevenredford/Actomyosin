package iterators;


import java.util.*;
import java.io.*;

import util.*;
import main.*;
import io.*;

public class ThreeBeadEvaluator extends Evaluator  {

	/** Histograms, etc for tracking filament position in Three-Bead Assay */
	Actin threeBeadActin;
	Myosin threeBeadMyo;
	HistogramPlus unboundPos;
	HistogramPlus boundPos;
	HistogramPlus tempBound;
	HistogramPlus localAvePos;
	HistogramPlus pt05Histo;
	HistogramPlus pt1Histo;
	HistogramPlus pt2Histo;
	HistogramPlus pt3Histo;
	HistogramPlus pt4Histo;
	boolean trackingBound = false;
	double bindingTimeStart = 0;
	double bindingTimeEnd = 0;
	final double bindingTimeThreshold = 0.2; // record binding events greater than this time
	final double threeBeadAssayTime = 30.0;	// write histos and exit run after this simulated time
	int histosAdded = 0;
	
	boolean firstEval = true;
	
	double nameOfInput;
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
	}
	
	public void reset()  {
		super.reset();
		firstEval = true;
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		unboundPos = new HistogramPlus (50,-25,25,dataPath,"unboundPos",null,false,true,false,false);
		boundPos = new HistogramPlus (50,-25,25,dataPath,"boundPos",null,false,true,false,false);
		tempBound = new HistogramPlus (50,-25,25,dataPath,"tempBound",null,false,true,false,false);
		localAvePos = new HistogramPlus (50,-25,25,dataPath,"localAvePos",null,false,true,false,false);
		pt05Histo = new HistogramPlus (50,-25,25,dataPath,"pt05Histo",null,false,true,false,false);
		pt1Histo = new HistogramPlus (50,-25,25,dataPath,"pt1Histo",null,false,true,false,false);
		pt2Histo = new HistogramPlus (50,-25,25,dataPath,"pt2Histo",null,false,true,false,false);
		pt3Histo = new HistogramPlus (50,-25,25,dataPath,"pt3Histo",null,false,true,false,false);
		pt4Histo = new HistogramPlus (50,-25,25,dataPath,"pt4Histo",null,false,true,false,false);
	}

	
	/** Override this method to specify whether this Evaluator writes data files.  These
	 * files MUST be written to the path defined by the variable dataPata to be compiled
	 * by the iterator.
	 *
	 * */
	public boolean hasData()  {
		return true;
	}

	public boolean stop(double tm)  {
		//if (Sim2D.simulationTime > threeBeadAssayTime) {
		if(tm > nextEval) {
			if (firstEval){
				new MyosinSurface(new Point2D(Sim2D.xDimension/2.,Sim2D.yDimension/2.));
				threeBeadMyo = Myosin.theMyosins[0];
				threeBeadActin = Actin.theActins[0];
				firstEval = false;
			} else {
		
		if (threeBeadMyo.isFree()) {
			if (trackingBound) {
				bindingTimeEnd = Sim2D.simulationTime;
				double boundTime = bindingTimeEnd - bindingTimeStart;
				if (boundTime > 0.05 & boundTime < 0.1) HistogramPlus.sumHistograms(tempBound, pt05Histo);
				if (boundTime > 0.1 & boundTime < 0.2) HistogramPlus.sumHistograms(tempBound, pt1Histo);
				if (boundTime > 0.2 & boundTime < 0.3) HistogramPlus.sumHistograms(tempBound, pt2Histo);
				if (boundTime > 0.3 & boundTime < 0.4) HistogramPlus.sumHistograms(tempBound, pt3Histo);
				if (boundTime > 0.4) HistogramPlus.sumHistograms(tempBound, pt4Histo);

				if (boundTime > bindingTimeThreshold) {
					HistogramPlus.sumHistograms(tempBound, boundPos);
					double localMean = tempBound.getMean();
					localAvePos.addValue(localMean);
					histosAdded++;
					System.out.println ("boundTime was " + boundTime + " ending at " + bindingTimeEnd + " with mean position of " + localMean);
				}
				trackingBound = false;
				tempBound.clearBins();
			} else {
				unboundPos.addValue(1500-threeBeadActin.cm.x);
			}
		} else {
			if (trackingBound) {
				tempBound.addValue(1500-threeBeadActin.cm.x);
			} else {
				trackingBound = true;
				bindingTimeStart = Sim2D.simulationTime;
				tempBound.addValue(1500-threeBeadActin.cm.x);
			}
		}
			}
			nextEval += interval;
		}
		if(super.stop(tm)) {
			boundPos.writeToFile();
			localAvePos.writeToFile();
			pt05Histo.writeToFile();
			pt1Histo.writeToFile();
			pt2Histo.writeToFile();
			pt3Histo.writeToFile();
			pt4Histo.writeToFile();
			return true;
		}
		return false;
	}
	
	
	
	public void cleanUp()  {}
			
}
