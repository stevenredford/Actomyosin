package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;


public class PolarOrMixedBundle  extends Initializer
{

	Vector templates = new Vector();


	public void init() throws Exception {
		for(Enumeration e = templates.elements(); e.hasMoreElements(); ) {
			((Template)e.nextElement()).init();
		}
	}

	public void loadParameter(String tag, AMInputStream in)  throws Exception {

		if(tag.equals("Filament"))  {
			Template t = new Template();
			t.loadTemplate(in);
			templates.add(t);
		}
		else throw new Exception("MakeSpecifiedActinFilaments.loadParameter(): got bad tag = " + tag);
	}

	class Template {

		double xPosition;
		double yPosition;
		int numInitialActins;
		boolean setAvgFilsInXsection;
		boolean constantActinSpacing;
		/*
		 * polarityBias=1 if all actin filaments in one direction, 0 if equal numbers in both directions
		 */
		double fractionBarbedEndsRight;
		double crosslinkDensity;
		double viscosityTweak=1;
		double textLoc = -1;
		FilamentAttachment minusAttachmentTemplate = null;
		FilamentAttachment plusAttachmentTemplate = null;
		

		public void loadTemplate(AMInputStream in)  throws Exception {

			String tag = in.nextTag();

			while(!tag.equals("endFilament")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}
				else if(tag.equals("InitYPosition"))  {
					yPosition = in.nextDouble();
				}

				else if(tag.equalsIgnoreCase("viscosityTweak"))  {
					viscosityTweak = in.nextDouble();
				}
				else if(tag.equalsIgnoreCase("crosslinkDensity"))  {
					crosslinkDensity = in.nextDouble();
				}
				else if(tag.equalsIgnoreCase("setAvgFilsInXsection"))  {
					setAvgFilsInXsection = in.nextBoolean();
				}
				
				else if(tag.equalsIgnoreCase("constantActinSpacing"))  {
					constantActinSpacing = in.nextBoolean();
				}
			
				else if(tag.equalsIgnoreCase("fractionBarbedEndsRight"))  {
					fractionBarbedEndsRight = in.nextDouble();
				}

				else if(tag.equals("TextAbove"))  {
					textLoc = 1;
				}

				else throw new Exception("MakeSpecifiedActinFilaments.Template.loadParameter(): got bad tag = " + tag);
				tag = in.nextTag();
			}
		}

		public void init() throws Exception {
			double bundleLeftX=Sim2D.xDimension/2+xPosition-MiscellaneousParameters.bundleLength/2;
			double bundleWidth;
			int monomers = (int)(MiscellaneousParameters.actinLength/Actin.monLength);
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			if (setAvgFilsInXsection){
				numInitialActins=(int)(MiscellaneousParameters.bundleLength/MiscellaneousParameters.actinLength*MiscellaneousParameters.avgFilsXsection);
				if (constantActinSpacing){
					bundleWidth=MiscellaneousParameters.bundleSpacing*(MiscellaneousParameters.avgFilsXsection-1);
				}
				else{
					bundleWidth=MiscellaneousParameters.bundleWidth;
				}
				
			}
			else{
				numInitialActins=Actin.numInitialActins;
				bundleWidth=MiscellaneousParameters.bundleWidth;
			}
			double bundleBottomY = Sim2D.yDimension/2+yPosition+bundleWidth/2;
			for(int i = 0; i < numInitialActins; i++){
				double initX = bundleLeftX+Math.random()*MiscellaneousParameters.bundleLength;
				double initY;
				if (setAvgFilsInXsection && constantActinSpacing){
					int rowNumber=(int) (MiscellaneousParameters.avgFilsXsection*Math.random());
					initY=bundleBottomY-rowNumber*MiscellaneousParameters.bundleSpacing;
					//System.out.println(rowNumber);
				}
				else {
					initY=bundleBottomY-Math.random()*MiscellaneousParameters.bundleWidth;
				}
				//System.out.println(initY);
				//double initY = bundleBottomY-Math.random()*MiscellaneousParameters.bundleWidth;
				double orientation;
				if (Math.random()<MiscellaneousParameters.fractionBarbedEndsRight){
					orientation=0;
				}
				else{
					orientation=Math.PI;
				}
				Actin nuActin = new Actin(initX, initY,orientation, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),null,null);
				nuActin.initPolarity=orientation;
				for(Monomer m=nuActin.pEndMonomer; m!=null; m=m.next){
					if (Math.random()<crosslinkDensity){
						new Crosslinker(m);
					}
				}
				nuActin.forceTextLocator = (int)textLoc;
				nuActin.setTweakDragFactor(viscosityTweak);
			}
		}

	}

}

