package iterators;


import java.util.*;
import java.io.*;

import main.*;
import io.*;
import util.*;

public class NetworkMorphologyEvaluator extends Evaluator {
	
	XLVariant  myXLVariant;
	String xlVariantType = null;
	double [] avgNumCrosslinks;
	double [] avgNumLinkers;
	double [] avgNumCrosslinksPerFil;
	double [] avgNumLinkersPerFil;
	double [] meanXlPerMonomer;
	double [] stdevXlPerMonomer;
	double [] bundledRatio;
	double eqlTime;
	int evalNum = 0;
	double nextBundleEval = 0;
	double bundleInterval = 10;
	double[] angleCorrFun1;
	double[] angleCorrFun2;
	double[] angleCorrFun3;
	double[] angleCorrFun4;
	double[] timeList;
	
	boolean [] checked;
	
	int []	numLinkers;
	
	int[][] linkerList;
	
	int [] bundleNumber;

	double nameOfInput;
	
	String bundleFileName;
	
	File bundleFile;
		
	PrintWriter bundlePW;
	
	public void init(String path, String name) throws Exception  {
		super.init(path,name);
		
		myXLVariant = Crosslinker.getXLVariant(xlVariantType);
		eqlTime = 3/(myXLVariant.getRecruitmentProb()*Actin.getTotalFilamentLength() + myXLVariant.getReleaseProb());
		
		int numEvals = (int)(Sim2D.runTime/interval) + 1;
		//angleCorrFun1 = new double[numEvals];
		//angleCorrFun2 = new double[numEvals];
		//angleCorrFun3 = new double[numEvals];
		//angleCorrFun4 = new double[numEvals];
		timeList = new double[numEvals];
		avgNumCrosslinks = new double[numEvals];
		avgNumLinkers = new double[numEvals];
		avgNumCrosslinksPerFil = new double[numEvals];
		avgNumLinkersPerFil = new double[numEvals];
		meanXlPerMonomer = new double[numEvals];
		stdevXlPerMonomer = new double[numEvals];
		bundledRatio = new double[numEvals];
		resetBundleData();
;
	}
		
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		bundleFileName = dataPath + File.separator + "bundledata.dat";
		makeDataFile();
	}

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equals("NameOfInput"))  {
			nameOfInput = in.nextDouble();
		}
		else if(tag.equals("BundleInterval"))  {
			bundleInterval = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("XLVariant")) {
			xlVariantType = in.nextString();
		}
		else super.loadParameter(tag, in);
	}
		
	public void reset()  {
		super.reset();
		nextBundleEval =0;
		resetBundleData();
		evalNum = 0;
	}
	
	public void resetBundleData() {
		checked = new boolean[Actin.actinCt];
		numLinkers = new int[Actin.actinCt];
		linkerList = new int[Actin.actinCt][200];
		bundleNumber = new int[Actin.actinCt];
	}
		
	
	public boolean stop(double tm)  {
		if((tm >= nextEval)) {
			avgNumCrosslinks[evalNum] = getNumberCrosslinks();
			avgNumLinkers[evalNum] = getNumberSinglyBoundLinkers();
			avgNumCrosslinksPerFil[evalNum] = avgNumCrosslinks[evalNum]/Actin.actinCt;
			avgNumLinkersPerFil[evalNum] = avgNumLinkers[evalNum]/Actin.actinCt;
			bundledRatio[evalNum]=detectBundles();
			evalNum++;
			nextEval += interval;
		}
		if(tm > nextBundleEval){
			detectBundles();
			nextBundleEval+=bundleInterval;
		}
		if(super.stop(tm)) {
			bundlePW.close();
			return true;
		}
		return false;
	}
	
	public void calculateBundleCorr(int evalNum) {
		Actin [] actins = Actin.theActins;
		int actinCt = Actin.actinCt;
		double cosTheta = 0;
		double rij = 0;
		double nFactor = 0;
		double cosAngle = 0;
		double cosAngler = 0;
		double cosAngle2 = 0;
		double cosAngle2r = 0;
		for (int i=0;i<actinCt;i++) {
			for (int j=0;j<actinCt;j++) {
				if(i!=j) {
				cosTheta=Point2D.dot(actins[i].uVect, actins[j].uVect);
				rij=Point2D.getDistance(actins[i].cm,actins[j].cm);
				//corrFun+=cosTheta*cosTheta/rij;
				cosAngle+=cosTheta;
				cosAngler+=cosTheta/rij;
				cosAngle2+=cosTheta*cosTheta;
				cosAngle2r+=cosTheta*cosTheta/rij;
				}
			}
		}
		nFactor = actinCt*(actinCt-1)*Math.sqrt(actinCt/(Sim2D.xDimension*Sim2D.yDimension));
		angleCorrFun1[evalNum] = cosAngle/(actinCt*(actinCt-1));
		angleCorrFun2[evalNum] = cosAngler/nFactor;
		angleCorrFun3[evalNum] = cosAngle2/(actinCt*(actinCt-1));
		angleCorrFun4[evalNum] = cosAngle2r/nFactor;
		timeList[evalNum] = Sim2D.simulationTime;
	}
	
	
	public double detectBundles() {
		Actin [] actins = Actin.theActins;
		int actinCt = Actin.actinCt;
		resetBundleData();
		scanActinLinkers();
		int bundleCounter = 1;
		//detect sequences of doubly linked filaments and associate sequence them with a bundle number
		for(int i=0;i<actinCt;i++){
			if (!checked[i]){
				detectBundledActin(i,bundleCounter);
				bundleCounter++;
			}
		}
		
		// bundle analysis
		double[] bundleFilNum = new double[bundleCounter];
		double[] bundleOrientation = new double[bundleCounter];
		for(int i=1;i<bundleCounter;i++){
			for(int j=0;j<actinCt;j++){
				if (bundleNumber[j]==i){
					bundleFilNum[i]++;
					bundleOrientation[i]+=Math.atan(actins[j].upVect.y/actins[j].upVect.x);
				}
			}
			bundleOrientation[i]/=bundleFilNum[i];
			//System.out.println("bundle number " + i + " has " + bundleFilNum[i] + " filaments at an average orientation: " + bundleOrientation[i]);
		}
		//System.out.println("number of bundles is " + (bundleCounter-1));
		
		//calculate the ratio of bundled actin filaments
		int numBundledFilaments = 0;
		for (int i=0;i<bundleCounter;i++){
			if(bundleFilNum[i]>1){
				numBundledFilaments+=bundleFilNum[i];
			}
		}
		
		writeBundleData();
		
			
		return (numBundledFilaments/Actin.actinCt);
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
				
	public String getHeaderString() {
		return new String("bundleratio" + "\t" + "#singlelinks" + "\t" + "#doublelinks" + "\t" + "#linkerpool" + "\t" + "#filaments" + "\t" + "filamentlength" + "\t" + "#monomers" + "\t" + "overlab");
	}
	
	public String getDataString() {
		int totalMonomers = 0;
		for (int i=0;i<Actin.actinCt;i++){
			totalMonomers += Actin.theActins[i].nMonomers;
		}
		double overlap = 2*Sim2D.xDimension*Sim2D.yDimension/(Actin.actinCt*Actin.avgActinLength);
		
		return new String(detectBundles() + "\t" + getNumberSinglyBoundLinkers() + "\t" + getNumberCrosslinks() + "\t" + Crosslinker.totalLinkerPool + "\t" + Actin.actinCt + "\t" + Actin.avgActinLength + "\t" + overlap);
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
		
	void detectBundledActin(int filID,int bundleCounter){
		if (checked[filID]) return;
		checked[filID] = true;
		bundleNumber[filID] = bundleCounter;
		int linker1FilId=0;
		int linker2FilId=0;
		for (int i=0;i<numLinkers[filID];i++){
			//System.out.println("thisActin.myNumOfLinkers=" + thisActin.myNumOfLinkers + " i=" + i);
			if (Crosslinker.theCrosslinkers[linkerList[filID][i]].isBound()){
				linker1FilId = getOtherActinId(filID,Crosslinker.theCrosslinkers[linkerList[filID][i]]);
				for (int j=i+1;j<numLinkers[filID];j++){
					//System.out.println("thisActin.myNumOfLinkers=" + thisActin.myNumOfLinkers + " i=" + i + " j=" + j);
					if (Crosslinker.theCrosslinkers[linkerList[filID][j]].isBound()){
						linker2FilId = getOtherActinId(filID,Crosslinker.theCrosslinkers[linkerList[filID][j]]);
						if(linker1FilId==linker2FilId){
							//System.out.println("actin number " + thisActin.myActinNumber + " is doubly linked to filament number " + linker1FilId + " linkers numbers= " + thisActin.myLinkersList[i] + " " + thisActin.myLinkersList[j]);
							detectBundledActin(linker1FilId,bundleCounter);
						}
					}
				}
			}
		}
	}
	
	
	static int getOtherActinId(int thisActinId,Crosslinker linker){
		if(thisActinId==linker.bound1.myFilament.myActinNumber){
			return linker.bound2.myFilament.myActinNumber;
		}
		else{
			return linker.bound1.myFilament.myActinNumber;
		}
	}
	
	
	public void detectNetwork() {
		resetBundleData();
		scanActinLinkers();
		int networkCounter = 1;
		//detect sequences of linked filaments and identify them with a sequence number
		for(int i=0;i<Actin.actinCt;i++){
			if (!checked[i]){
				detectLinkedActin(i,networkCounter);
				networkCounter++;
			}
		}
		//bundle analysis
		double[] networkFilNum = new double[networkCounter];
		double[] networkOrientation = new double[networkCounter];
		for(int i=1;i<networkCounter;i++){
			for(int j=0;j<Actin.actinCt;j++){
				if (bundleNumber[j]==i){
					networkFilNum[i]++;
					networkOrientation[i]+=Math.atan(Actin.theActins[j].upVect.y/Actin.theActins[j].upVect.x);
				}
			}
			networkOrientation[i]/=networkFilNum[i];
			//System.out.println("bundle number " + i + " has " + bundleFilNum[i] + " filaments at an average orientation: " + bundleOrientation[i]);
		}
		System.out.println("number of network is " + (networkCounter-1));
		
		// write the correlation function to a file in the current run directory
		String networkPath = basePath + File.separator + "networkdata." + String.valueOf(Sim2D.simulationTime);
		int jj=1;
		File networkFile = new File(networkPath + "." + String.valueOf(jj) + ".dat");
		
		while (networkFile.isDirectory()) {
			System.out.println (networkFile.getName() + " exists.... keeping file name BUT changing directory name");
			networkFile = new File(networkPath + "." + String.valueOf(jj) + ".dat");
			jj++;
		}
		try {
			FileWriter networkFW = new FileWriter(networkFile);
			PrintWriter networkPW = new PrintWriter(networkFW,true);
			networkPW.printf("%16s %16s %16s","bundle number","Num of Filaments","Orientation");
			networkPW.print("\n");
			for(int j=1;j<networkCounter;j++) {
				networkPW.printf("%16d %16.6f %16.6f",j,networkFilNum[j],networkOrientation[j]);
				networkPW.print("\n");  // print a new line character.
			}
			networkPW.flush();
			networkPW.close();
		} catch (IOException ioe) { System.out.println("An error creating dataFile"); }
	}
	
	
	public void scanActinLinkers(){
		int filID;
		for (int i=0;i<Crosslinker.crosslinkerCt;i++) {
			if(Crosslinker.theCrosslinkers[i].bound1!=null){
				filID = Crosslinker.theCrosslinkers[i].bound1.myFilament.myActinNumber;
				linkerList[filID][numLinkers[filID]]=i;
				numLinkers[filID]++;
			}
			if(Crosslinker.theCrosslinkers[i].bound2!=null){
				filID = Crosslinker.theCrosslinkers[i].bound2.myFilament.myActinNumber;
				linkerList[filID][numLinkers[filID]]=i;
				numLinkers[filID]++;
			}
		}
	}
	
	void detectLinkedActin(int filID,int bundleCounter){
		if (checked[filID]) return;
		checked[filID] = true;
		bundleNumber[filID] = bundleCounter;
		for (int i=0;i<numLinkers[filID];i++){
			if (Crosslinker.theCrosslinkers[linkerList[filID][i]].isBound()){
				int otherActin = getOtherActin(filID,Crosslinker.theCrosslinkers[linkerList[filID][i]]);
				detectLinkedActin(otherActin,bundleCounter);
			}
		}
	}
	
	int getOtherActin(int filID,Crosslinker linker){
		if(filID==linker.bound1.myFilament.myActinNumber){
			return linker.bound2.myFilament.myActinNumber;
		}
		else{
			return linker.bound1.myFilament.myActinNumber;
		}
	}
	
	public void cleanUp()  {}
	
	
}

