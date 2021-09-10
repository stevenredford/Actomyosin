package initializers;
/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class SixFilamentMixedAlignedBundleSpecifiedOrientations extends Initializer{
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
		


		/*
		 * polarityBias=1 if all actin filaments in one direction, 0 if equal numbers in both directions
		 */
		
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

				else if(tag.equals("TextAbove"))  {
					textLoc = 1;
				}
				else if(tag.equals("MinusEndAttachment"))  {
					String attach_name = "no_name";
					try {
						attach_name = in.nextString();
						Class c = Class.forName("main." + attach_name);
						minusAttachmentTemplate = (FilamentAttachment)c.newInstance();
						minusAttachmentTemplate.load(in);
					} catch(Exception e) {
						System.out.println("MakeSpecifiedActinFilaments.loadParameter: Problem loading MinusAttachment " + attach_name + ":" + e.toString());
						throw(e);
					}
				}
				else if(tag.equals("PlusEndAttachment"))  {
					String attach_name = "no_name";
					try {
						attach_name = in.nextString();
						Class c = Class.forName("main." + attach_name);
						plusAttachmentTemplate = (FilamentAttachment)c.newInstance();
						plusAttachmentTemplate.load(in);
					} catch(Exception e) {
					System.out.println("MakeSpecifiedActinFilaments.loadParameter: Problem loading PlusAttachment " + attach_name + ":" + e.toString());
						throw(e);
					}
				}

				else throw new Exception("MakeSpecifiedActinFilaments.Template.loadParameter(): got bad tag = " + tag);
				tag = in.nextTag();
			}
		}

		public void init() throws Exception {
			
			double [] yLocations={-2.5*MiscellaneousParameters.bundleSpacing,-1.5*MiscellaneousParameters.bundleSpacing,-0.5*MiscellaneousParameters.bundleSpacing,0.5*MiscellaneousParameters.bundleSpacing,1.5*MiscellaneousParameters.bundleSpacing,2.5*MiscellaneousParameters.bundleSpacing};
			int [] usedBarbedEndLocations=new int[6];
			
			if (MiscellaneousParameters.barbedEndLocationIndicator==1){
				int[] barbedEndLocations={1,1,1,0,0,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==2){
				int[] barbedEndLocations={1,1,0,1,0,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==3){
				int[] barbedEndLocations={1,1,0,0,1,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==4){
				int[] barbedEndLocations={1,1,0,0,0,1};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==5){
				int[] barbedEndLocations={1,0,1,1,0,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==6){
				int[] barbedEndLocations={1,0,1,0,1,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==7){
				int[] barbedEndLocations={1,0,1,0,0,1};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==8){
				int[] barbedEndLocations={0,1,1,1,0,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==9){
				int[] barbedEndLocations={0,1,1,0,1,0};
				usedBarbedEndLocations=barbedEndLocations;
			}
			else if(MiscellaneousParameters.barbedEndLocationIndicator==10){
				int[] barbedEndLocations={0,1,1,0,0,1};
				usedBarbedEndLocations=barbedEndLocations;
			}
			
			
		
			
			
			FilamentAttachment mAtt, pAtt;
			if(minusAttachmentTemplate != null) {
				mAtt = minusAttachmentTemplate.makeNewInstance();
			}
			else mAtt = null;
			
			if(plusAttachmentTemplate != null) {
				pAtt = plusAttachmentTemplate.makeNewInstance();
			}
			else pAtt = null;
			double initX = Sim2D.xDimension/2+xPosition;
			int monomers = (int)(MiscellaneousParameters.actinLength/Actin.monLength);
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			
			for (int i=0; i<=5; i++){
				double orientation=0;
				double initY=Sim2D.yDimension/2+yLocations[i];
				//System.out.println(initY);
				if (usedBarbedEndLocations[i]==0){
					orientation=0;
				}
				else if (usedBarbedEndLocations[i]==1){
					orientation=Math.PI;
				}
				Actin nuActin = new Actin(initX, initY,orientation, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
				nuActin.initPolarity=orientation;
				System.out.println(nuActin.cm.y);
				nuActin.forceTextLocator = (int)textLoc;
				nuActin.setTweakDragFactor(viscosityTweak);
			}
			
	
			

		}

	}

}
