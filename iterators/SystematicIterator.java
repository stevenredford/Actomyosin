
package iterators;


import java.util.*;
import java.io.*;
import java.net.*;

import io.*;
import parameters.*;

public class SystematicIterator extends Iterator {

	int numSamplesPerPar;
			
	public void loadParameter(String tag, AMInputStream in) throws Exception {
		if(tag.equals("NumSamplesPerPar")) {
			numSamplesPerPar = in.nextInt();
		}
		else {
			super.loadParameter(tag, in);
		}
	}
		
	/*  One step in a recursive iteration that samples a regular grid in n-dimensional parameter space, writing each
	set of sampled parameter values to an output file for some client to read.
	 */
	private void oneStep(int depth, int [] numLevels, Parameter [] pars, RecursionRider r) throws Exception {
		for(int i = 0; i < numLevels[depth]; i++)  {
			pars[depth].setDiscreteLevel(i,numLevels[depth]);
			if(depth == pars.length-1) {
				if(r != null) {
					issueNextJob(r);
					r.cnt++;
				}
			}
			else
				oneStep(depth+1,numLevels,pars, r);
		}
	}
	
	public void serveFiles() throws Exception {
		
		File [] list = (new File(basePath)).listFiles(new FileFilter() {
			public boolean accept(File f) {
				return f.getName().contains(baseName) && f.isDirectory();
			}}
		);
		File f;
		RecursionRider r = new RecursionRider(0,list);
	
		Parameter [] pars = getParameterArray();
		int [] numSamples = getNumSamplesPerPar(pars);
		
		for(int i = 0; i < numIteratesToMake; i++) {
			oneStep(0,numSamples,pars,r);
		}
	
		boolean allDone = false;
		while(!allDone)  {
			allDone = true;
			try {
				Thread.sleep(100);
			}
			catch(Exception e) {}
			for(int i = 0; i < r.dirs.length; i++) {
				f = new File(r.dirs[i] + File.separator + "CLIENT_READY");
				if(f.exists()) {
					f = new File(r.dirs[i] + File.separator + "ALL_DONE");
					if(!f.exists()) {
						f.createNewFile();
					}

				}
				else {
					allDone = false;
				}
			}
		}
			try {
				Thread.sleep(1000);
			}
			catch(Exception e) {}
		collectData(pars,r.dirs);
	}
	
	int [] getNumSamplesPerPar(Parameter [] pars) {
		int [] nSamples = new int[pars.length];
		for(int i = 0; i < pars.length; i++)  {
			if(pars[i].getNumSamples() > 0) {
				nSamples[i] = pars[i].getNumSamples();
			}
			else {
				nSamples[i] = numSamplesPerPar;
			}
		}
		return nSamples;
	}
		
	private File getNextAvailableClient(RecursionRider r) {
		File ready;
		for(int i = 0; i < r.dirs.length; i++) {
			ready = new File(r.dirs[i].getAbsolutePath() + File.separator + "CLIENT_READY");
			if(ready.exists()) {
				return r.dirs[i];
			}
		}
		return null;
	}
	
	
	private void issueNextJob(RecursionRider r) throws Exception{
		File client = getNextAvailableClient(r);
		while(client == null)  {
			try {
				Thread.sleep(100);
			}
			catch(Exception e) {}
			client = getNextAvailableClient(r);
		}
		File f = new File(client.getPath() + File.separator + "CLIENT_READY");
		f.delete();
		f = new File(client.getPath() + File.separator + "tmp");
		AMOutputStream out = new AMOutputStream(new FileOutputStream(f));
		writeParamsToVary(out);
		out.close();
		f.renameTo((new File(client.getPath() + File.separator + "PARFILE.par")));
	}
	
	
	public void setParameter(String name, double val) {
		Parameter par;
		for(Enumeration e = paramsToVary.elements(); e.hasMoreElements();) {
			par = (Parameter)e.nextElement();
			if(par.nameMatches(name)) {
				par.setValue(val);
				Parameters.setParameter(name,val);
				return;
			}
		}
	}

		
	public void writeParamsToVary(AMOutputStream out) {
		Parameter p;
		out.printTag("ParameterValues");
		out.increaseIndent();
		try {
			for(Enumeration e = paramsToVary.elements(); e.hasMoreElements(); ) {
				p = (Parameter)e.nextElement();
				out.printTag(p.getName(),p.getDoubleValue());
			}
		} catch(Exception e) {}
		out.decreaseIndent();
		out.printTag("endParameterValues");
	}
	
	

	public class RecursionRider  {
		
		int cnt = 0;
		File [] dirs;

		public RecursionRider(int cnt, File [] dirs) {
			this.cnt = 0;
			this.dirs = dirs;
		}
		
	}
}
