
package iterators;


import java.util.*;
import java.io.*;
import java.net.*;

import io.*;
import parameters.*;

public class ParStreamIterator extends Iterator {

	String 	parFileToSample = null;
	
	AMInputStream parameterStream;
	
		
	public void loadParameter(String tag, AMInputStream in) throws Exception {
		if(tag.equalsIgnoreCase("ParSampleFile")) {
			parFileToSample = in.nextString();
			sampleMode = FROM_PAR_FILE;
		}
		else if(tag.equals("ParametersToVary")) {
			String parName = in.nextTag();
			while(!parName.equals("endParametersToVary")) {
				if(!Parameters.isValidParameterName(parName)) {
					throw new Exception("Iterator.loadParameter(): got bad name for parameter to vary = " + parName);
				}
				Parameter p = new Parameter("Iterator",parName,Parameter.DOUBLE);
				addParameter(p);
				parName = in.nextTag();
			}
		}
		else {
			super.loadParameter(tag, in);
		}
	}
	
	public void serveFiles() throws Exception {
		File [] dirs = (new File(basePath)).listFiles(new FileFilter() {
			public boolean accept(File f) {
				return f.getName().contains(baseName) && f.isDirectory();
			}}
		);
		File f;
		
		Parameter [] pars = getParameterArray();
		issueJobs(pars,dirs);

		boolean allDone = false;
		while(!allDone)  {
			allDone = true;
			try {
				Thread.sleep(100);
			}
			catch(Exception e) {}
			for(int i = 0; i < dirs.length; i++) {
				f = new File(dirs[i] + File.separator + "CLIENT_READY");
				if(f.exists()) {
					f = new File(dirs[i] + File.separator + "ALL_DONE");
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
		collectData(pars,dirs);
	}
	
	private void issueJobs(Parameter [] pars, File [] dirs) throws Exception {
	
		if(parFileToSample == null)  {
			throw new Exception("ParStreamIterator:  must specify a ParFileToSample to use this iterator!");
		}

		parameterStream = new AMInputStream(new File(basePath + File.separator + parFileToSample));
		
		try {
			while(true) {
				getNextParameters(pars);
				issueNextJob(dirs);
			}
		} catch(EOFException e) {
			System.out.println(e.toString());
		}
	}
	
	private void getNextParameters(Parameter [] pars) throws Exception {
		for(int i = 0; i < pars.length; i++) {
			pars[i].setValue(parameterStream.nextDouble());
		}
		parameterStream.skipLine();
	}

		
	
	private File getNextAvailableClient(File [] dirs) {
		File ready;
		for(int i = 0; i < dirs.length; i++) {
			ready = new File(dirs[i].getAbsolutePath() + File.separator + "CLIENT_READY");
			if(ready.exists()) {
				return dirs[i];
			}
		}
		return null;
	}
	
	
	private void issueNextJob(File [] dirs) throws Exception{
		File client = getNextAvailableClient(dirs);
		while(client == null)  {
			try {
				Thread.sleep(100);
			}
			catch(Exception e) {}
			client = getNextAvailableClient(dirs);
		}
		File f = new File(client.getPath() + File.separator + "CLIENT_READY");
		f.delete();
		f = new File(client.getPath() + File.separator + "tmp");
		AMOutputStream out = new AMOutputStream(new FileOutputStream(f));
		writeParamsToVary(out);
		out.close();
		f.renameTo((new File(client.getPath() + File.separator + "PARFILE.par")));
	}
}
