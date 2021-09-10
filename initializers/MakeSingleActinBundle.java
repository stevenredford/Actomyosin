package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class MakeSingleActinBundle  extends Initializer
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
		double length;
		boolean crosslinkers;
		FilamentAttachment minusAttachmentTemplate = null;
		FilamentAttachment plusAttachmentTemplate = null;


		public void loadTemplate(AMInputStream in)  throws Exception {

			String tag = in.nextTag();

			while(!tag.equals("endFilament")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}

				else if(tag.equals("Length"))  {
					length = in.nextDouble();
				}

				else if(tag.equals("bunchWidth"))  {
					bunchWidth = in.nextDouble();
				}
				else if(tag.equals("Crosslinkers"))  {
					crosslinkers = in.nextBoolean();
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
			FilamentAttachment mAtt, pAtt;
			double initX = Sim2D.xDimension/2+xPosition;
			double initY = Sim2D.yDimension/2;
			int monomers = (int)(length/Actin.monLength);
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			Monomer m1180 = new Monomer();
			Monomer m2180 = new Monomer();
			Monomer m10 = new Monomer();
			Monomer m20= new Monomer();
			int whichMonomer180 =1;
			int whichMonomer0=1;
			for (int i=0; i<Actin.numInitialActins;i++){
				if(minusAttachmentTemplate != null) {
					mAtt = minusAttachmentTemplate.makeNewInstance();
				}
				else mAtt = null;

				if(plusAttachmentTemplate != null) {
					pAtt = plusAttachmentTemplate.makeNewInstance();
				}
				else pAtt = null;

				if (Actin.numInitialActins>1){

					initY=Sim2D.yDimension/2+bunchWidth/2-i*bunchWidth/(Actin.numInitialActins);
				}

				Actin newActin = new Actin(initX, initY,Math.PI, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
				if (crosslinkers){
					if (whichMonomer180==-1){
					m1180= newActin.bEndMonomer;
					}
					else {
						m1180=newActin.bEndMonomer.prev;
					}
					if (i>0){
						new Crosslinker (m1180, m2180);
					}
					if (whichMonomer180==1){
						m2180=newActin.bEndMonomer;
					}
					else {
						m2180=newActin.bEndMonomer.prev;
					}
					whichMonomer180*=-1;
				}
				
			}
			for (int i=0; i<Actin.numInitialActins;i++){
				if(minusAttachmentTemplate != null) {
					mAtt = minusAttachmentTemplate.makeNewInstance();
				}
				else mAtt = null;

				if(plusAttachmentTemplate != null) {
					pAtt = plusAttachmentTemplate.makeNewInstance();
				}
				else pAtt = null;
				if (Actin.numInitialActins>1){

					initY=Sim2D.yDimension/2+bunchWidth/2-i*bunchWidth/(Actin.numInitialActins);
				}

				Actin newActin = new Actin(initX, initY,0, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
				if (crosslinkers){
					if (whichMonomer0==-1){
					m10= newActin.bEndMonomer;
					}
					else {
						m10=newActin.bEndMonomer.prev;
					}
					if (i>0){
						new Crosslinker (m10, m20);
					}
					if (whichMonomer0==1){
						m20=newActin.bEndMonomer;
					}
					else {
						m20=newActin.bEndMonomer.prev;
					}
					whichMonomer0*=-1;
				}
			}
		}

	}

}



