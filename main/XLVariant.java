/**
 * XLInfoPacket.java
 *
 * @author Created by Omnicore CodeGuide
 */

package main;

import java.awt.Color;

import io.*;
import iterators.*;

public class XLVariant implements CrosslinkListener
{
	String name = null;
	
	/** Size of total pool for this variant. */
	double totalPool = -1;
	
	/** Rest length for force calculations. */
	double restLength = -1;
	
	/** Stiffness of this linker. */
	double springConstant = -1;
	
	/** Probability this linker will bind a filament from bulk solution. */
	double recruitmentProb = -1;
	
	/** Probability a singly-bound linker will bind a second filament if it is within reach. */
	double bindingProb = -1;
	
	/** multiplier for force-based linker release. */
	double forceBasedReleaseFac = -1;
	
	/** linkers only bind filament if contact angle exceeds this angle. */
	double minBindingAngle = -1;
	
	/** Maximum stretch before a linker releases. */
	double maxStretch = -1;
	
	/** Probabality that a bound linker will release a bond. */
	double releaseProb = -1;
	
	double crosslinkerCt;
	
	/** an array of probabilities of binding to a monomer in a filament based on the distance from this linker */
	double[][] distanceProbTable;
	
	double distInterval = 0.1;
	
	Color ColorRod = Color.black;
	Color ColorFree = Color.black;
	Color ColorBound = Color.black;

	
	public XLVariant(String n) {
		name = n;
		Crosslinker.addBondListener(this);
	}
	
	public void init(double totalLinkerPool,
					 double xlRestLength,
					 double xlSpringConstant,
					 double xlBindingProb,
					 double xlRecruitmentProb) {
		
		if(totalPool < 0) {
			totalPool = totalLinkerPool;
		}
		if(bindingProb < 0) {
			bindingProb = xlBindingProb;
		}
		if(restLength < 0) {
			restLength = xlRestLength;
		}
		if(springConstant < 0) {
			springConstant = xlSpringConstant;
		}
		if(recruitmentProb < 0) {
			recruitmentProb = xlRecruitmentProb;
		}
		
		makeProbDistanceTable();
		crosslinkerCt = 0;
	}
	
	
	public String getName() {
		return name;
	}
		
	public double getSpringConstant() {
		return springConstant;
	}
	
	public boolean hasSpringConstant() {
		return springConstant >= 0;
	}
	
	public double getReleaseProb() {
		return releaseProb;
	}
	
	public boolean hasReleaseProb() {
		return releaseProb >= 0;
	}
	
	public double getMaxStretch() {
		return maxStretch;
	}
	
	public boolean hasMaxStretch() {
		return maxStretch >= 0;
	}
	
	public double getForceBasedReleaseFac() {
		return forceBasedReleaseFac;
	}
	
	public boolean hasForceBasedReleaseFac() {
		return forceBasedReleaseFac >= 0;
	}
	
	public double getMinBindingAngle() {
		return minBindingAngle;
	}
	
	public boolean hasMinBindingAngle() {
		return minBindingAngle >= 0;
	}
	
	
	public double getRestLength() {
		return restLength;
	}
	
	public boolean hasRestLength() {
		return restLength >= 0;
	}

	public double getRecruitmentProb() {
		return recruitmentProb;
	}
	
	public boolean hasRecruitmentProb() {
		return recruitmentProb >= 0;
	}
	
	public double getBindingProb() {
		return bindingProb;
	}
	
	public boolean hasBindingProb() {
		return bindingProb >= 0;
	}
	
	public double [][] getDistanceProbabilityTable() {
		return distanceProbTable;
	}
	
	public double getDistInterval() {
		return distInterval;
	}
	
	//color of XL
	public Color getColorRod() {
		return ColorRod;
	}
	public boolean hasColorRod() {
		return ColorRod != Color.black;
	}
	public Color getColorFree() {
		return ColorFree;
	}
	public boolean hasColorFree() {
		return ColorFree != Color.black;
	}
	public Color getColorBound() {
		return ColorBound;
	}
	public boolean hasColorBound() {
		return ColorBound != Color.black;
	}
	
	private void makeProbDistanceTable(){
		int nBins = (int)((Crosslinker.xlMaxBindingDistance)/distInterval)+2;
		int maxMonNumber = 2*(int)Math.floor(Crosslinker.xlMaxBindingDistance/Actin.monLength)+1;
		distanceProbTable = new double[nBins][maxMonNumber+3];
		double[][] realProbTable = new double[nBins][maxMonNumber];
		double Ai = 0;
		double Rij = 0;
		int nMonomers = 0;
		double probTot;
		for (int i=0; i<nBins; i++){
			Ai = Math.sqrt(Math.pow(Crosslinker.xlMaxBindingDistance, 2)-Math.pow(i*distInterval, 2));
			nMonomers = (int)Math.floor(Ai/Actin.monLength);
			distanceProbTable[i][0] = nMonomers;
			distanceProbTable[i][1] = i*distInterval;
			probTot=0;
			for (int j=0; j<nMonomers; j++){
				Rij = Math.sqrt(Math.pow(i*distInterval,2)+Math.pow(j*Actin.monLength,2));
				realProbTable[i][j] = bindingProb*Math.exp(-Math.pow((Rij-restLength),2)*springConstant/Constants.kT);
				probTot+=2*realProbTable[i][j];
				distanceProbTable[i][j+3]=probTot;
			}
			distanceProbTable[i][2] += distanceProbTable[i][nMonomers+2];
		}

	}
	
	/** returns a random choice from the wieghted distribution d, assuming that entries of d are distributed in (0,1]
	 * in monotnically increasing order.
	 */
	public int selectFromWeightedDistribution(double [] d) {
		double choice = Math.random();
		for(int i = 0; i < d.length; i++) {
			if(choice < d[i]) {
				return i;
			}
		}
		return -1;
	}
			

	public void doCrossLinkerBinding(double [] monomerDist,double dT) {
		Monomer m;
		Crosslinker c;
		double bindingProbability = dT*recruitmentProb*Actin.getTotalFilamentLength();
		if (bindingProbability>1){
			bindingProbability=1;
		}
		double totalBound = bindingProbability*(totalPool - crosslinkerCt);

		int n = (int)Math.floor(totalBound);
		for(int i = 0; i < n; i++) {
			int choice = selectFromWeightedDistribution(monomerDist);
			m = Actin.theActins[choice].bindCrossLinker();
			if(m != null) {
				addCrossLinker(m);
			}
		}
		
		double residual = totalBound - n;
		if(Math.random() < residual) {
			int choice = selectFromWeightedDistribution(monomerDist);
			m = Actin.theActins[choice].bindCrossLinker();
			if(m != null) {
				addCrossLinker(m);
			}
		}
	}
	
	public void addCrossLinker(Monomer m) {
		Crosslinker c = new Crosslinker(m);
		c.setXLVariant(this);
		crosslinkerCt++;
	}
		

	public void singlyBoundLinkerReleased(double lifeTime) {
		crosslinkerCt--;
	}
	
	public void doublyBoundLinkerReleased(double lifeTime) {}
	
	
	final public void load(AMInputStream in)  throws Exception {
		
		String tag = in.nextTag();
		while(!tag.equals("endXLVariant"))  {
			loadParameter(tag, in);
			tag = in.nextTag();
		}
	}
	/** Override this method to load in specific parameters for XLVariant.  The basic format
	 * to use is:
	 *
		if(tag.equals("Parameter1"))  {
			variable1 = in.nextDouble();
		}
		else if(tag.equals("Parameter2"))  {
			variable2 = in.nextDouble();
		}
	 * .
	 * .
	 * .
	 * 	else if(tag.equals("ParameterN"))  {
			variableN = in.nextDouble();
		}

	 *
	 *
	 *  */
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equalsIgnoreCase("name"))  {
			name = in.nextString();
		}
		else if(tag.equalsIgnoreCase("SpringConstant"))  {
			springConstant = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("Restlength"))  {
			restLength = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("xlRecruitmentProb"))  {
			recruitmentProb = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("totalPool"))  {
			totalPool = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("ColorRod"))  {
			ColorRod = in.nextColor();
		}
		else if(tag.equalsIgnoreCase("ColorFree"))  {
			ColorFree = in.nextColor();
		}
		else if(tag.equalsIgnoreCase("ColorBound"))  {
			ColorBound = in.nextColor();
		}

		else throw new Exception("XLVariant.loadParameter(): got bad tag = " + tag);
	}
	
}

