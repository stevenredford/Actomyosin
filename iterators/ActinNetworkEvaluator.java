package iterators;



import java.util.*;
import java.io.*;

import main.*;
import util.*;
import io.*;

public class ActinNetworkEvaluator extends Evaluator implements ActinListener  {

	double meanLength;
	double stdvLength;
	double meanLifetime;
	double stdvLifetime;
	double iterNum = 0;
	
	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus lifetimeHistogram;
	
	/** A histogram that keeps running track of distribution of filament lengths.. */
	HistogramPlus lengthHistogram;
	
				
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		super.loadParameter(tag, in);
	}
	
	public void init(String path, String name) throws Exception {
		super.init(path,name);
		Actin.setMonomerListener(this);
	}

		public void reset()  {
		super.reset();
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		lifetimeHistogram = new HistogramPlus(20,0,100,dataPath,"actinLifetimes",null,false,true,false,false);
		lengthHistogram = new HistogramPlus(50,0,0.5*Sim2D.xDimension,dataPath,"actinLengths",null,false,true,false,false);
	}
	
	public void monomerReleased(Monomer m) {
		double lifetime = Sim2D.simulationTime - m.birthTime;
		if(m.birthTime > 0) {
			lifetimeHistogram.addValue(lifetime);
		}
	}
		
	void updateLengthHistogram() {
		lengthHistogram.clearBins();
		for (int i=0;i<Actin.actinCt;i++) {
			lengthHistogram.addValue(Actin.theActins[i].physicalLength);
		}
	}

	
	void writeHistograms() {
		lengthHistogram.writeToFile();
		lifetimeHistogram.writeToFile();
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
		if(tm > nextEval) {
			updateLengthHistogram();
			iterNum++;
			nextEval += interval;
		}
		if(super.stop(tm)) {
			meanLifetime = lifetimeHistogram.getMean();
			stdvLifetime = lifetimeHistogram.getStdev();
			meanLength = lengthHistogram.getMean();
			stdvLength = lengthHistogram.getStdev();
			iterNum = 0;
			writeHistograms();
			return true;
		}
		return false;
	}
	
	
	public String getHeaderString() {
		return new String("MeanLifetime" + "\t" + "StdevLifetime" + "\t" + "MeanLength" + "\t"+ "StdevLength");
	}
	
	public String getDataString() {
		return new String(meanLifetime + "\t" + stdvLifetime + "\t" + meanLength + "\t" + stdvLength);
	}
	
	public void cleanUp()  {}
			
}
