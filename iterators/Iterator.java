package iterators;

import java.util.*;
import java.io.*;
import java.net.*;

import main.*;
import io.*;
import parameters.*;
;

public class Iterator  {
	
	static String
		RANDOM = new String("Random"),
		FROM_PAR_FILE = new String("FromParFile");

	/** The Sim2D object that runs a simulation. */
	Sim2D theRunner;
	
	/** A list of the evaluators that from the gauntlet each parameter set has to run. */
	Vector theEvaluators = new Vector();
		
	/** A list of the parameters to vary. */
	Vector paramsToVary = new Vector();
	
	/** Directory from which the search was launched. */
	String basePath = "";
	
	/** base name of input file. */
	String baseName = "";
	
	String errorPath = null;
	
	PrintStream errorOut = null;
	
	File errorFile;
	
	int numIteratesToMake = 1;
	
	String sampleMode = RANDOM;
		
	AMOutputStream out;
	
	public Iterator()  {}
	
	public void init(Sim2D runner, String runMode) throws Exception  {
		setBasePath(runner.getBasePath());
		setBaseName(runner.getBaseName());
		setRunner(runner);
		if(!(runMode.equals(Sim2D.SERVE) || runMode.equals(Sim2D.GATHER))) {
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator) en.nextElement();
				eval.init(basePath,baseName);
			}
		}
		else  {
			theRunner.initialize();
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator) en.nextElement();
				eval.getName();
			}
		}
				
	}
	
	public void setOutputStream(AMOutputStream out) {
		this.out = out;
	}
	
	public void makeErrorPath(int cnt)  throws Exception  {
		File errorDir = new File(basePath + File.separator + "error" + cnt);
		if(!errorDir.mkdir()) {
			throw new Exception("Iterator.setErrorPath: Can't make new directory for path = " + basePath + File.separator + "error" + cnt);
		}
		errorPath = errorDir.getAbsolutePath();
		errorFile = new File(errorPath+File.separator + baseName + cnt + ".err");
		errorOut = new PrintStream(errorFile);
		System.setOut(errorOut);
	}
		
	public void reset(int cnt) {
		theRunner.reset();
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator) en.nextElement();
			eval.reset();
		}
	}
	
	public void iterate() throws Exception {
		boolean endOfParamStream = false;
		int cnt = 0;
		while(!endOfParamStream && cnt < numIteratesToMake)  {
			try  {
				getNextParameters();
				doOneIterate(cnt++);
			} catch(EOFException e)  {
				endOfParamStream = true;
				System.out.println("iterator stream reached eof");
			}
		}
	}

	public void doOneIterate(int cnt)  {
		try {
			theRunner.setMoviePath(new String("movie" + cnt));
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator) en.nextElement();
				eval.setOutputPaths(cnt);
			}
			makeErrorPath(cnt);
			theRunner.run();
			errorOut.close();
		} catch(Exception e) {
			recordError(cnt,e);
			reset(cnt);
			errorOut.close();
			return;
		}
		if(gotHit()) {
			recordHit(cnt);
		}
		reset(cnt);
	}
	
	private boolean gotHit() {
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator) en.nextElement();
			if(!eval.gotHit()) {
				return false;
			}
		}
		return true;
	}

	
	public void recordError(int cnt, Exception e)  {
		out.printTag("Error",cnt);
		out.increaseIndent();
		out.printTag("ParameterValues",paramsToString());
		out.println(e.toString());
		out.printStackTrace(e);
		out.decreaseIndent();
		out.printTag("endError");
	}

	public void recordHit(int cnt)  {
		out.printTag("Hit",cnt);
		out.increaseIndent();
		out.printTag("ParameterValues",paramsToString());
		
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator) en.nextElement();
			eval.recordHit(out,cnt);
		}
		if(Sim2D.writeMovieFrames) {
			out.printTag("MoviePath",new String("\"" + theRunner.getMoviePath() + "\""));
			String path = new String(basePath + File.separator + baseName + "_ff_" + cnt + ".png");
			out.printTag("FinalFrame",new String("\"" + path + "\""));
			theRunner.writeRemoteFrame(path);
		}
		out.printTag("ErrorPath",new String("\"" + errorPath + "\""));
		out.decreaseIndent();
		out.printTag("endHit");
	}
	
	public void gatherFiles() throws Exception {
		File [] list = (new File(basePath)).listFiles(new FileFilter() {
			public boolean accept(File f) {
				return f.getName().contains(baseName) && f.isDirectory();
			}}
		);
		for(int i = 0; i < list.length; i++) {
			System.out.println(list[i]);
		}
		Parameter [] pars = getParameterArray();
		for(int i = 0; i < pars.length; i++) {
			System.out.println(pars[i]);
		}
		collectData(pars,list);
	}
		
	
	void collectData(Parameter [] pars, File [] dirs) throws Exception {
				
		String pathPrefix = new String("");
		String targetPath = new String(basePath + File.separator + baseName);
		File targetDir = new File(targetPath);
		targetDir.mkdir();
		
		AMInputStream in;
		AMOutputStream spreadsheet = new AMOutputStream(new FileOutputStream(targetPath + File.separator + "ss.txt"));
		AMOutputStream summary = new AMOutputStream(new FileOutputStream(targetPath + File.separator + "summary.txt"));
		File oldPath, newPath,dataDir, ffDir, movieDir, errorDir, evalDir, parFile;
		int cnt = 0;
		
		
		String ss_data = new String("cnt\t");
		for(int i = 0; i < pars.length; i++) {
			ss_data = ss_data + pars[i].getName() + "\t";
		}
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator)en.nextElement();
			ss_data = ss_data + eval.getHeaderString() + "\t";
		}
		spreadsheet.println(ss_data);
		
		
		for(int i = 0; i < dirs.length; i++) {
			try {
				in = new AMInputStream(new File(dirs[i] + File.separator + "iteratorOut"));
				String tag;
				while(true) {
					while(!(tag = in.nextTag()).equals("Hit")) {}
					summary.printTag("Hit",cnt);
					summary.increaseIndent();
					ss_data = new String(cnt + "\t");

					while(!(tag = in.nextTag()).equals("endHit")) {
						if(tag.equals("ParameterValues")) {
							double val;
							String parString = new String("");
							for(int j = 0; j < pars.length; j++) {
								val = in.nextDouble();
								parString = parString + val + "\t";
							}
							summary.printTag("ParameterValues",parString);
							ss_data = ss_data + parString;
						}
						else if(tag.equals("Data")) {
							String dataString = in.nextString();
							ss_data = ss_data + dataString + "\t";
							summary.printTag(pathPrefix + "_Data",new String("\"" + dataString + "\""));
						}
						else if(tag.equals("DataPath")) {
							String s = in.nextString();
							oldPath = new File(s.substring(s.lastIndexOf(baseName),s.length()));
							dataDir = new File(targetPath + File.separator + "DataFiles");
							if(!dataDir.exists()) {
								dataDir.mkdir();
							}
							newPath = new File(dataDir + File.separator + baseName + "_data_" + cnt);
							oldPath.renameTo(newPath);
							summary.printTag(pathPrefix + "_DataPath",new String("\"" + newPath + "\""));
								
						}
						else if(tag.equals("FinalFrame")) {
							String s = in.nextString();
							oldPath = new File(s.substring(s.indexOf(baseName),s.length()));
							ffDir = new File(targetPath + File.separator + "final_frames");
							if(!ffDir.exists()) {
								ffDir.mkdir();
							}
							newPath = new File(ffDir + File.separator + baseName + "_ff_" + cnt + ".png");
							oldPath.renameTo(newPath);
							summary.printTag(pathPrefix + "_FinalFrame",new String("\"" + newPath + "\""));
						}
						else if(tag.equals("MoviePath")) {
							String s = in.nextString();
							oldPath = new File(s.substring(s.lastIndexOf(baseName),s.length()));
							movieDir = new File(targetPath + File.separator + "movies");
							if(!movieDir.exists()) {
								movieDir.mkdir();
							}
							newPath = new File(movieDir + File.separator + baseName + "_movie_" + cnt);
							oldPath.renameTo(newPath);
							summary.printTag(pathPrefix + "_MoviePath",new String("\"" + newPath + "\""));
						}
						else if(tag.equals("ErrorPath")) {
							String s = in.nextString();
							oldPath = new File(s.substring(s.lastIndexOf(baseName),s.length()));
							errorDir = new File(targetPath + File.separator + "errors");
							if(!errorDir.exists()) {
								errorDir.mkdir();
							}
							newPath = new File(errorDir + File.separator + baseName + "_error_" + cnt);
							oldPath.renameTo(newPath);
							summary.printTag(pathPrefix + "_ErrorPath",new String("\"" + newPath + "\""));
						}
						else if(tag.equals("Evaluator")) {
							String evalName = in.nextString();
							pathPrefix = new String(evalName);
							targetPath = new String(basePath + File.separator + baseName + File.separator + evalName);
							evalDir = new File(targetPath);
							if(!evalDir.exists()) {
								evalDir.mkdir();
							}
						}
						else if(tag.equals("endEvaluator")) {
							pathPrefix = new String("");
							targetPath = new String(basePath + File.separator + baseName);
						}
					}
					spreadsheet.println(ss_data);
					summary.decreaseIndent();
					summary.printTag("endHit");
					cnt++;
				}
			} catch(EOFException e) {System.out.println(e.toString());}
		}
		spreadsheet.close();
		summary.close();
		for(int i = 0; i < dirs.length; i++) {
			System.out.println("deleting " + dirs[i].getPath());
			ream(dirs[i]);
		}
	}
	
	public void ream(File f) {
		if(f.isDirectory()) {
			File [] list = f.listFiles();
			for(int i = 0; i < list.length;i++) {
				ream(list[i]);
			}
			f.delete();
		}
		else f.delete();
	}

	
	public void serveFiles() throws Exception {
		FileOps.reportln("Sorry. need a sampleParfile for this iterator to serve files");
	}
	
	private void getNextParameters() throws Exception {
		if(sampleMode.equals(RANDOM)) {
			Parameter p;
			for(Enumeration e = paramsToVary.elements(); e.hasMoreElements();) {
				p = (Parameter) e.nextElement();
				Parameters.setParameter(p.getName(),p.getRandomValue());
			}
		}
		else throw new Exception("Iterator.getNextParameters: bad value for sampleMode = " + sampleMode);
	}
	
	public void readNextParameters(AMInputStream in) throws Exception {
		String name;
		double val;
		while(!(name = in.nextTag()).equals("ParameterValues")) {}
		
		while(!(name = in.nextTag()).equals("endParameterValues")) {
			val = in.nextDouble();
			setParameter(name,val);
		}
	}
	
	Parameter [] getParameterArray() {
		Parameter [] pars = new Parameter[paramsToVary.size()];
		int i = 0;
		for(Enumeration e = paramsToVary.elements(); e.hasMoreElements(); ) {
			pars[i++] = (Parameter)e.nextElement();
		}
		return pars;
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
		
							
	public void runAsClient() throws Exception {
		boolean done = false;
		int cnt = 0;
		File ready = new File(basePath + File.separator + "CLIENT_READY");
		File parfile = new File(basePath + File.separator + "PARFILE.par");
		File alldone = new File(basePath + File.separator + "ALL_DONE");
		ready.createNewFile();
		while(!done) {
			while(!alldone.exists() && !parfile.exists())  {
				try {
					Thread.sleep(100);
				}
				catch(Exception e) {}
			}
			if(alldone.exists()) {
				
				done = true;
			}
			else { // parfile exists
				readNextParameters(new AMInputStream(parfile));
				doOneIterate(cnt++);
				parfile.delete();
				ready.createNewFile();
			}
		}
	}
	
	public void load(AMInputStream in)  throws Exception {
		
		String tag = in.nextTag();
		while(!tag.equals("endIterator"))  {
			loadParameter(tag, in);
			tag = in.nextTag();
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
	
	public String paramsToString() {
		Parameter p;
		String s = new String("");;
		try {
			for(Enumeration e = paramsToVary.elements(); e.hasMoreElements(); ) {
				p = (Parameter)e.nextElement();
				s = s + p.getDoubleValue() + "\t";
			}
			return s;
		} catch(Exception e) { return new String("bad");}
	}

	public void loadParameter(String tag, AMInputStream in) throws Exception  {
			
		if(tag.equals("Evaluator"))  {
			String eval_name = "no_name";
			try {
				eval_name = in.nextString();
				Class c = Class.forName("iterators." + eval_name);
				Evaluator eval = (Evaluator)c.newInstance();
				eval.load(in);
				addEvaluator(eval);
			} catch(Exception e) {
				System.out.println("Problem loading evaluator " + eval_name + ":" + e.toString());
				e.printStackTrace();
			}
		}
		else if(tag.equals("ParametersToVary")) {
			String parName = in.nextTag();
			while(!parName.equals("endParametersToVary")) {
				if(parName.equalsIgnoreCase("add")) {
					parName = in.nextString();
					if(Parameters.isValidParameterName(parName)) {
						throw new Exception("Iterator.loadParameter(): Trying to add a new parameter with name that already exists = " + parName);
					}
					Parameter p = new Parameter("Iterator",parName,Parameter.DOUBLE);
					p.loadFromLine(in);
					Parameters.addParameter(p);
					addParameter(p);
				}
				else if(parName.equalsIgnoreCase("CompositeParameter")) {
					CompositeParameter p = new CompositeParameter();
					p.loadFromBlock(in);
					Parameters.addParameter(p);
					addParameter(p);
				}
					
				else {
					if(!Parameters.isValidParameterName(parName)) {
						throw new Exception("Iterator.loadParameter(): got bad name for parameter to vary = " + parName);
					}
					Parameter p = new Parameter("Iterator",parName,Parameter.DOUBLE);
					p.loadFromLine(in);
					addParameter(p);
				}
				parName = in.nextTag();
			}
		}
		else if(tag.equals("NumIteratesToMake")) {
			numIteratesToMake = in.nextInt();
		}
	}

	void addParameter(Parameter p) {
		paramsToVary.add(p);
	}
	
	public void setBasePath(String path) {
		basePath = path;
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator) en.nextElement();
			eval.setBasePath(path);
		}
	}
		
	public void setBaseName(String name) {
		baseName = name;
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator) en.nextElement();
			eval.setBaseName(name);
		}
	}
		
	
	public Vector getEvaluators()  { return theEvaluators; }
	
	public void addEvaluator(Evaluator e) { theEvaluators.add(e); }
	
	public void setRunner(Sim2D runner) {
		theRunner = runner;
	}
	

	
	

				
		
	
	
}
