package initializers;

/**
 * Initializer.java
 *
 * @author Created by Omnicore CodeGuide
 */

import java.util.*;

import io.*;
import main.*;

public class MakeSpecifiedActinBundles  extends Initializer
{
	
	Vector templates = new Vector();
	
	public void init() throws Exception {
		FilamentAttachment mAtt, pAtt;
		for(Enumeration e = templates.elements(); e.hasMoreElements(); ) {
			((Template)e.nextElement()).init();
		}
	}
	
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equals("Bundle"))  {
			Template tm = new Template();
			tm.loadTemplate(in);
			templates.add(tm);
		}
		else throw new Exception("MakeSpecifiedActinBundles.loadParameter(): got bad tag = " + tag);
	}
	
	class Template {
		
		double xPosition;
		double yPosition;
		double orientation;
		double length;
		FilamentAttachment minusAttachmentTemplate = null;
		FilamentAttachment plusAttachmentTemplate = null;
	
		public void loadTemplate(AMInputStream in)  throws Exception {
			
			String tag = in.nextTag();
			
			while(!tag.equals("endBundle")) {
				if(tag.equals("InitXPosition"))  {
					xPosition = in.nextDouble();
				}
				else if(tag.equals("InitYPosition"))  {
					yPosition = in.nextDouble();
				}
				else if(tag.equals("Length"))  {
					length = in.nextDouble();
				}
				else if(tag.equals("Orientation"))  {
					orientation = 2.0*Math.PI*in.nextDouble()/360.0;
				}
			else if(tag.equals("MinusEndAttachment"))  {
				String attach_name = "no_name";
				try {
					attach_name = in.nextString();
					Class c = Class.forName("main." + attach_name);
					minusAttachmentTemplate = (FilamentAttachment)c.newInstance();
					minusAttachmentTemplate.load(in);
				} catch(Exception e) {
					System.out.println("MakeSpecifiedActinBundles.loadParameter: Problem loading MinusAttachment " + attach_name + ":" + e.toString());
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
				System.out.println("MakeSpecifiedActinBundles.loadParameter: Problem loading PlusAttachment " + attach_name + ":" + e.toString());
					throw(e);
				}
			}
				else throw new Exception("MakeSpecifiedActinBundles.Template.loadParameter(): got bad tag = " + tag);
				tag = in.nextTag();
			}
		}
		
		public void init() throws Exception {
			FilamentAttachment mAtt, pAtt;
			double initX = Sim2D.xDimension/2+xPosition;
			double initY = Sim2D.yDimension/2+yPosition;
			int monomers = (int)(length/Actin.monLength);
			if (monomers > Actin.maxMon) { monomers = Actin.maxMon; }
			
			if(minusAttachmentTemplate != null) {
				mAtt = minusAttachmentTemplate.makeNewInstance();
			}
			else mAtt = null;
			
			if(plusAttachmentTemplate != null) {
				pAtt = plusAttachmentTemplate.makeNewInstance();
			}
			else pAtt = null;

			new Bundle(initX, initY,orientation, monomers,monomers,Actin.actinDynamicsMode.equals(Actin.STATIC_FILAMENTS),mAtt,pAtt);
		}

	}

}

