package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class MakeAnchoredActinBundle  extends Initializer
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
		else throw new Exception("MakeSpecifiedActinBundle.loadParameter(): got bad tag = " + tag);
	}

	class Template {

		double xPosition;
		double bunchWidth;
		int numFilaments;
		double bundleLength;
		double anchorLength;
		boolean crosslinkers;
		FilamentAttachment rightAnchorAttachmentTemplate = null;
		FilamentAttachment leftAnchorAttachmentTemplate = null;


		public void loadTemplate(AMInputStream in)  throws Exception {

			String tag = in.nextTag();

			while(!tag.equals("endFilament")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}
				else if(tag.equals("BundleLength"))  {
					bundleLength = in.nextDouble();
				}
				else if(tag.equals("AnchorLength"))  {
					anchorLength = in.nextDouble();
				}

				else if(tag.equals("bunchWidth"))  {
					bunchWidth = in.nextDouble();
				}

				else if(tag.equals("RightAnchorAttachment"))  {
					String attach_name = "no_name";
					try {
						attach_name = in.nextString();
						Class c = Class.forName("main." + attach_name);
						rightAnchorAttachmentTemplate= (FilamentAttachment)c.newInstance();
						rightAnchorAttachmentTemplate.load(in);
					} catch(Exception e) {
						System.out.println("MakeSpecifiedActinFilaments.loadParameter: Problem loading RightAnchorAttachment " + attach_name + ":" + e.toString());
						throw(e);
					}
				}
				else if(tag.equals("LeftAnchorAttachment"))  {
					String attach_name = "no_name";
					try {
						attach_name = in.nextString();
						Class c = Class.forName("main." + attach_name);
						leftAnchorAttachmentTemplate = (FilamentAttachment)c.newInstance();
						leftAnchorAttachmentTemplate.load(in);
					} catch(Exception e) {
						System.out.println("MakeSpecifiedActinFilaments.loadParameter: Problem loading LeftAnchorAttachment " + attach_name + ":" + e.toString());
						throw(e);
					}
				}
				else throw new Exception("MakeSpecifiedActinFilaments.Template.loadParameter(): got bad tag = " + tag);
				tag = in.nextTag();
			}
		}

		public void init() throws Exception {
			FilamentAttachment rAtt, lAtt;
			double initX = Sim2D.xDimension/2+xPosition;
			double initY = Sim2D.yDimension/2;
			int monomers = (int)(bundleLength/Actin.monLength);
			int anchorMonomers = (int)(anchorLength/Actin.monLength);
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
		
			if (rightAnchorAttachmentTemplate != null){
				rAtt = rightAnchorAttachmentTemplate.makeNewInstance();
			}
			else rAtt=null;
			
			if (leftAnchorAttachmentTemplate != null){
				lAtt = leftAnchorAttachmentTemplate.makeNewInstance();
			}
			else lAtt=null;
			Actin rightAnchor = new Actin(initX+bundleLength/2-anchorLength/2,initY,0, anchorMonomers,anchorMonomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),null,rAtt);
			Actin leftAnchor = new Actin(initX-bundleLength/2+anchorLength/2,initY,Math.PI, anchorMonomers,anchorMonomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),null,lAtt);
			Monomer rightAnchorMonomer= rightAnchor.bEndMonomer;
			Monomer leftAnchorMonomer=leftAnchor.bEndMonomer;
			for (int i=0; i<Actin.numInitialActins;i++){
			

				if (Actin.numInitialActins>1){

					initY=Sim2D.yDimension/2+bunchWidth/2-i*bunchWidth/(Actin.numInitialActins);
				}

				Actin newActin = new Actin(initX,initY,0,monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),null,null);
				Monomer crosslinkedMonomer=newActin.bEndMonomer;
				if (i>0){
					rightAnchorMonomer=rightAnchorMonomer.prev;
				}
				for (int j=0; j<i; j++){
					crosslinkedMonomer=crosslinkedMonomer.prev;
				}
				new Crosslinker(rightAnchorMonomer,crosslinkedMonomer);
			
				
			}
			for (int i=0; i<Actin.numInitialActins;i++){
		
				if (Actin.numInitialActins>1){

					initY=Sim2D.yDimension/2+bunchWidth/2-i*bunchWidth/(Actin.numInitialActins);
				}

				Actin newActin = new Actin(initX,initY,Math.PI,monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),null,null);
				Monomer crosslinkedMonomer=newActin.bEndMonomer;
				if (i>0){
					leftAnchorMonomer=leftAnchorMonomer.prev;
				}
				for (int j=0; j<i; j++){
					crosslinkedMonomer=crosslinkedMonomer.prev;
				}
				new Crosslinker(leftAnchorMonomer,crosslinkedMonomer);
			}
		}

	}

}

