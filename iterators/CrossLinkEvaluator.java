package iterators;


import java.util.*;
import java.io.*;

import main.*;
import util.*;
import io.*;

public class CrossLinkEvaluator extends Evaluator implements CrosslinkListener  {

	XLVariant  myXLVariant;
	String xlVariantType = null;
	double avgNumCrosslinks = 1;
	double avgNumLinkers = 2;
	double meanBondLifetime = 3;
	double meanXLLifetime = 3;
	double meanXlPerMonomer = 4;
	double stdevXlPerMonomer = 5;
	static double eqlTime;
	double iterNum = 0;
	
	/** A histogram that keeps running track of distribution of filament subunit lifetimes.. */
	HistogramPlus crosslinkLifetimes;
	
	/** A histogram that keeps running track of linker/Factin bond lifetimes. */
	HistogramPlus bondLifetimes;
			
	/** A histogram that keeps running track of singly-bound linkers per unit length of filament. */
	HistogramPlus linkersPerMonomer;
	
	/** A histogram that keeps running track of crosslinks per unit length of filament. */
	HistogramPlus crosslinksPerMonomer;

	String bundleFileName;
	
	File bundleFile;
		
	PrintWriter bundlePW;

	double nameOfInput;
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equals("NameOfInput"))  {
			nameOfInput = in.nextDouble();
		}
		else if(tag.equals("BundleInterval"))  {
			interval = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("XLVariant")) {
			xlVariantType = in.nextString();
		}
		else super.loadParameter(tag, in);
	}
	
	/** Called at the beginning of a single iteration. */
	public void init(String path, String name) throws Exception  {
		super.init(path,name);
		Crosslinker.addBondListener(this);
		myXLVariant = Crosslinker.getXLVariant(xlVariantType);
		eqlTime = 3/(myXLVariant.getRecruitmentProb()*Actin.getTotalFilamentLength() + myXLVariant.getReleaseProb());

	}

	public void reset()  {
		super.reset();
	}
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		crosslinkLifetimes = new HistogramPlus(50,0,3./myXLVariant.getReleaseProb(),dataPath,"xlLifetimes",null,false,true,false,false);
		bondLifetimes = new HistogramPlus(50,0,3./myXLVariant.getReleaseProb(),dataPath,"bondLifetimes",null,false,true,false,false);
		linkersPerMonomer = new HistogramPlus(200,0,0.5,dataPath,"linkersPerFilament",null,false,true,false,false);
		crosslinksPerMonomer = new HistogramPlus(200,0,0.1,dataPath,"crosslinksPerFilament",null,false,true,false,false);
		bundleFileName = dataPath + File.separator + "bundledata.dat";
		makeDataFile();
	}

	private void makeDataFile() throws Exception {
		bundleFile = new File(bundleFileName);
		FileWriter bundleFW = new FileWriter(bundleFile);
		bundlePW = new PrintWriter(bundleFW,true);
	}
	
	private void writeBundleData() {
		bundlePW.println("I am some output");
		bundlePW.flush();
	}

	
	/** Override this method to specify whether this Evaluator writes data files.  These
	 * files MUST be written to the path defined by the variable dataPata to be compiled
	 * by the iterator.
	 *
	 * */
	public boolean hasData()  {
		return bundleFile != null;
	}

	public boolean stop(double tm)  {
		if((tm > nextEval)&&(tm > eqlTime)) {
			updateLinkerHistograms();
			if (iterNum==0){
				avgNumCrosslinks = getNumberCrosslinks();
				avgNumLinkers = getNumberSinglyBoundLinkers();
			} else {
				avgNumCrosslinks += getNumberCrosslinks();
				avgNumLinkers += getNumberSinglyBoundLinkers();
			}
			iterNum++;
			writeBundleData();
			nextEval += interval;
		}
		if(super.stop(tm)) {
			avgNumCrosslinks = avgNumCrosslinks/iterNum;
			avgNumLinkers = avgNumLinkers/iterNum;
			iterNum = 0;
			meanBondLifetime = bondLifetimes.getMean();
			meanXLLifetime = crosslinkLifetimes.getMean();
			meanXlPerMonomer = crosslinksPerMonomer.getMean();
			stdevXlPerMonomer = crosslinksPerMonomer.getStdev();
			bundlePW.close();
			return true;
		}
		return false;
	}
	
	public void singlyBoundLinkerReleased(double lifeTime) {
		bondLifetimes.addValue(lifeTime);
	}
	
	public void doublyBoundLinkerReleased(double lifeTime) {
		crosslinkLifetimes.addValue(lifeTime);
	}
	
	private void updateLinkerHistograms() {
		int crosslinkCt,linkerCt;
		for(int i = 0; i < Actin.actinCt; i++ ) {
			crosslinkCt = linkerCt = 0;
			for(Monomer m = Actin.theActins[i].pEndMonomer; m != null; m = m.next) { // run through the linked list of monomers
				if(m.boundLinker!=null) {
					if(m.boundLinker.isBound()) {
						crosslinkCt++;
					}
					else { // must be singly bound...
						linkerCt++;
					}
				}
			}
			crosslinksPerMonomer.addValue(((double)crosslinkCt)/Actin.theActins[i].nMonomers);
			linkersPerMonomer.addValue(((double)linkerCt)/Actin.theActins[i].nMonomers);
		}
	}

	private double getNumberCrosslinks() {
		double ct = 0;
		for(int i = 0; i < Crosslinker.crosslinkerCt; i++) {
			if(Crosslinker.theCrosslinkers[i].isBound())
				ct++;
		}
		return ct;
	}
		
	private double getNumberSinglyBoundLinkers() {
		double ct = 0;
		for(int i = 0; i < Crosslinker.crosslinkerCt; i++) {
			if(Crosslinker.theCrosslinkers[i].isFree())
				ct++;
		}
		return ct;
	}
	
	public String getHeaderString() {
		return new String("AvgNumCrosslinks" + "\t" + "AvgNumLinkers" + "\t" + "MeanBondLifetime" + "\t"+ "MeanXLLifetime" + "\t" + "meanXlPerActin" + "\t" + "stdevXlPerActin");
	}
	
	public String getDataString() {
		return new String(avgNumCrosslinks + "\t" + avgNumLinkers + "\t" + meanBondLifetime + "\t" + meanXLLifetime + "\t" + meanXlPerMonomer + "\t" + stdevXlPerMonomer);
	}
	
	public void cleanUp()  {}
			
}
