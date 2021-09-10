package main;

/*
	Sim2D... the file containing "main" and the simulation doLoop for a 2D psuedo-biological simulation
	units are in nm, pN, seconds
*/

import java.text.*;
import java.util.Date;
import java.util.Random;
import java.util.Vector;
import java.util.Enumeration;
import java.io.*;
import java.net.*;
import java.awt.image.*;
import java.awt.*;

import gui.*;
import io.*;
import parameters.*;
import util.Point2D;
import iterators.*;
import analysis.*;
import initializers.*;
import collision.*;

public class Sim2D {
	
	static String className = new String("main.Sim2D");

	/** This code gets executed as soon as the program is launched and calls a
	 * static method called staticInit for each agent class, allowing it to do
	 * some upfront initializations - creating Parameter objects, registering
	 * ParameterSetListeners, and creating and registering AgentProxy objects.
	 */
	static {
		agentProxies = new Vector();
		staticInit();
		Actin.staticInit();
		Constants.staticInit();
		Myosin.staticInit();
		MyosinV.staticInit();
		MyosinVTail.staticInit();
		MyosinMiniFilament.staticInit();
		MyosinSurface.staticInit();
		Crosslinker.staticInit();
		BarrierElement.staticInit();
		MyosinSurface.staticInit();
		Monomer.staticInit();
		MyosinHolder.staticInit();
		ElasticAttachment.staticInit();
		MiscellaneousParameters.staticInit();
	}
	
	/** A list of AgentProxy objects that receive commands from Sim2D and transmit them
	 * through static methods calls to the agent classes they represent. */
	static Vector agentProxies;
	
	/** A list of Initializer objects that do different initialization tasks before
	 * each simulation run. */
	Vector initializers = new Vector();
	
	/** The one instance of a Sim2D that is instantiated at launch time. */
	static public Sim2D MC;
	
	/** The frame that graphics are drawn into. */
	static WorldFrame theFrame;
	
	/** Graphics rendering control panel. */
	static Random generator=new Random();

	/**  Max resolvable X distance in a periodic world. */
	static public double maxUniqueXDist;
	
	/** Max resolvable Y distance in a periodic world. */
	static public double maxUniqueYDist;
	
	/** Currrent simulation tinme in seconds. */
	static public double simulationTime = 0;
	
	/** Number of simulation time steps performed. */
	static public int counter = 0;
	
	/** An arbitrarily large set of Evaluator objects that get called during each timestep.*/
	Vector theEvaluators = new Vector();
	
	String runMode = RUN;
	
	static public String
		ITERATE = new String("iterate"),
		RUN = new String("run"),
		SERVE = new String("serve"),
		CLIENT = new String("client"),
		GATHER = new String("gather");
	

	/************* file writing info ***************/
	
	/** Counts elapsed sim time between paints. */
	static double timeSinceLastPaint = 0;
	
	/** elapsed sim time since last movie frame. */
	static double timeSinceLastFrame = 0;
	
	static String movieFileName = "frm";
	
	/** The root directory that is created at start time to store all output for a run. */
	static public String basePath;
	
	/** The basename used to build directory names for e.g movie and data output distectories. */
	static public String baseName;
	
	/** The name of the directory to which movie data is written. */
	static String moviePath;
			
	static DecimalFormat fileIdFormat = new DecimalFormat ("#000000000.#;#000000000.#");
	static DecimalFormat timeOutFormat = new DecimalFormat ("#0.00#; #0.00#");	// a clean format for simulation time printing
	static DecimalFormat aveLengthFormat = new DecimalFormat ("#0.#; #0.#");	// a format for ave length printing
	
	// for run control
	static public boolean running = false;
	
	/** An iterator object which (if it exists) controls a sequence of simulations for
	 * different choices of parameter values. */
	static public Iterator theIterator;
	
	/** The object that handles writing image frames. */
	static RemoteImageWriter remoteImageWriter;
	
	public static void main (String args[]) {
		new Sim2D (args);
	}
		
	public Sim2D (String [] inputParams) {
		String s;
		try {
			if(inputParams.length == 1) {
				if(inputParams[0].equalsIgnoreCase("parameters")) {
					Parameters.showParameters();
				}
				else if((s = FileOps.getExtension(inputParams[0])) != null && s.equals("par")) {
					setWriteMovieFrames(true);
					setBaseName(extractBaseNameFromParfile(inputParams[0]));
					loadInputFile(inputParams[0]);
					runMode = RUN;
					makeOutputPaths(inputParams[0],runMode);
					run();
				}
			}
			else if(inputParams.length == 2)  {
				if((s = FileOps.getExtension(inputParams[1])) != null && s.equals("par")) {
					setRunMode(inputParams[0]);
					setWriteMovieFrames(false);
					setBaseName(extractBaseNameFromParfile(inputParams[1]));
					loadInputFile(inputParams[1]);
					if(theIterator == null) {
						throw new Exception("Sim2D: inputfile must specify an iterator");
					}
					makeOutputPaths(inputParams[1],runMode);
					theIterator.setBaseName(baseName);
					runIterator();
				}
			}
			else {
				System.out.println("usage:");
			}
		} catch(Exception e) {
			System.out.println("Sim2D runtime error");
			System.out.println(e.toString());
			e.printStackTrace();
			return;
		}
	}
	
	private void runIterator() throws Exception {
		theIterator.init(this, runMode);
		theEvaluators = theIterator.getEvaluators();
		running = true;
		if(runMode.equals(ITERATE)) {
			theIterator.iterate();
		}
		else if(runMode.equals(CLIENT)) {
			theIterator.runAsClient();
		}
		else if(runMode.equals(SERVE)) {
			theIterator.serveFiles();
		}
		else if(runMode.equals(GATHER)) {
			theIterator.gatherFiles();
		}
	}
	
	private void setRunMode(String mode) throws Exception {
		if(mode.equals(RUN) || mode.equals(ITERATE) || mode.equals(SERVE) || mode.equals(CLIENT) || mode.equals(GATHER)) {
			runMode = mode;
		}
		else throw new Exception("Sim2D.setRunMode: got bad runmode = " + mode);
	}
	
	private void loadInputFile(String fileName) throws Exception {
		AMInputStream in = new AMInputStream(new File(fileName));
		String tag;
		try {
			while(true) {
				tag = in.nextTag();
				if(tag.equals("Iterator")) {
					String iter_name = "no_name";
					try {
						iter_name = in.nextString();
						Class c = Class.forName("iterators." + iter_name);
						Iterator iter = (Iterator)c.newInstance();
						iter.load(in);
						theIterator = iter;
					} catch(Exception e) {
						System.out.println("Problem loading iterator " + iter_name + ":" + e.toString());
						throw(e);
					}
				}
				else if(tag.equals("Evaluator"))  {
					String eval_name = "no_name";
					try {
						eval_name = in.nextString();
						Class c = Class.forName("iterators." + eval_name);
						Evaluator eval = (Evaluator)c.newInstance();
						eval.load(in);
						theEvaluators.add(eval);
					} catch(Exception e) {
						System.out.println("Problem loading evaluator " + eval_name + ":" + e.toString());
						e.printStackTrace();
					}
				}
				else if(tag.equals("Parameters")) {
					Parameters.loadInputParams(in);
				}
				else if(tag.equals("Initializers")) {
					loadInitializers(in);
				}
				else if(tag.equals("XLVariant")) {
					Crosslinker.loadXLVariant(in);
				}
				else if(tag.equals("SimulationInstanceName")) {
					setBaseName(in.nextString());
				}
				else throw new Exception("Sim2D.loadInputFile(): unknown tag = " + tag);
			}
		} catch(EOFException e) {
			return;  // means we loaded the file correctly.
		}
	}
	

	private void loadInitializers(AMInputStream in) throws Exception {
		String tag = in.nextTag();
		while(!tag.equals("endInitializers")) {
			if(tag.equals("Initializer")) {
				String init_name = "no_name";
				try {
					init_name = in.nextString();
					Class c = Class.forName("initializers." + init_name);
					Initializer i = (Initializer)c.newInstance();
					i.load(in);
					initializers.add(i);
				} catch(Exception e) {
					System.out.println("Sim2D.loadInitializers: Problem loading initializer " + init_name + ":" + e.toString());
					throw(e);
				}
			}
			tag = in.nextTag();
		}
	}
	
	static public void addProxy(AgentProxy p) {
		agentProxies.add(p);
	}

	
	public void run() throws Exception {
		init();
		doLoop();
	}
		
	public void init() throws Exception {
		Mesh.create();
		initialize();
		initDiagnostics();
		initCollisions();
		
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator)en.nextElement();
			eval.reset();
			eval.resetGraphics(remote);
		}
		
		if (!remote) {
			// make the frame to which we draw
			theFrame = new WorldFrame (true);
			theFrame.showAll();
		}
		if(FileOps.verbose)
			Parameters.printParametersUsed();

	}
	
	public void initialize() throws Exception {
		for(Enumeration e = initializers.elements(); e.hasMoreElements(); ) {
			((Initializer)e.nextElement()).init();
		}
	}
		
	
	
	/** Give all the individual classes a chance to create static variables before first simulition run. */
	public void initDiagnostics() {
		for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
			((AgentProxy)e.nextElement()).initDiagnostics();
		}
	}

	/** Give all the individual classes a chance to create static variables before first simulition run. */
	public void initCollisions() {
		for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
			((AgentProxy)e.nextElement()).init();
		}
	}
	/** Give all the individual classes a chance to reset static variables before next simulation run. */
	public void reset() {
		for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
			((AgentProxy)e.nextElement()).reset();
		}
		Thing.resetAll();
		if (!remote) {
			if(theFrame!=null)
				theFrame.dispose();
		}
		simulationTime = 0;
		counter = 0;
	}
	
	private boolean evaluatorsDontStop(double simulationTime) {
		Evaluator eval;
		for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
			eval = (Evaluator)en.nextElement();
			if(eval.stop(simulationTime)) {
				return false;
			}
		}
		return true;
	}
			
		
	public void doLoop() {
		MyosinHolder.rootMyoHolder = MyosinMiniFilament.theMiniFilaments[0];
				
		FileOps.reportln("Beginning doLoop..");
		
		// this first sleep just helps avoid freeze on cluster startup (thread creation related I think)
		try { Thread.sleep(1000); } catch (InterruptedException e) { System.out.println ("error sleeping while starting do loop"); }
		
		while (simulationTime <= runTime && evaluatorsDontStop(simulationTime) ) {
			if (!running) {
				try { Thread.sleep(1000); } catch (InterruptedException e) { System.out.println ("error sleeping while paused"); }
			} else {
							
				checkTimeStepLimits();
				
				/* CHECK COLLISIONS:
				 Agents that need to check collisions with other agents, do so here.*/
				for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
					((AgentProxy)e.nextElement()).checkCollisions(simulationTime);
				}
				

				/* DO FORCES:
				 Agents (e.g motors and crosslinkers that need to exert/exchange force
				 on/with other agents, do so here.*/
				for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
					((AgentProxy)e.nextElement()).doForces(actualDeltaT);
				}
				
				/*  STEP:
				 Agents that need to move and/or change kinetic state do so here.*/
				for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
					((AgentProxy)e.nextElement()).step(actualDeltaT);
				}
				
				/* BRING OUT YOUR DEAD:
				 Any agents that were marked for removal (e.g. unbound crosslinkers etc)
				 during this timestep .*/
				Thing.removeArrayedDeadThings();
				
				/* FILL COLLISION DETECTION MESH:
				 Agents that need to register themselves in a collision detection
				 mesh do so here.*/
				for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
					((AgentProxy)e.nextElement()).registerForCollisionDetection();
				}

				/* CHECK BARRIER CROSSINGS:
				 Agents that need to detect and undo barrier crossings register themselves in a collision detection
				 mesh do so here.*/
				for(Enumeration e = agentProxies.elements(); e.hasMoreElements(); ) {
					((AgentProxy)e.nextElement()).checkBarrierCrossings();
				}
								
				simulationTime += actualDeltaT;
				if (!remote) {
					if (paintOn) {
						if(timeSinceLastPaint>=paintStep||simulationTime < actualDeltaT){
							theFrame.showAll();
							for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
								((Evaluator)en.nextElement()).showGraphics();
							}
							
							try { Thread.sleep(2); } catch (InterruptedException e) { System.out.println ("error sleeping while paused"); }
							timeSinceLastPaint=0;
							//try { Thread.sleep(1); } catch (InterruptedException e) { System.out.println ("error sleeping after paint"); }
						}
						else{
							timeSinceLastPaint+=actualDeltaT;
							theFrame.showInfoOnly();
						}
					} else { theFrame.showInfoOnly(); }		// only update Thing painting if paintOn
					if (writeMovieFrames) {
						try { Thread.sleep(0); } catch (InterruptedException e) { System.out.println ("error sleeping while paused"); }
						checkPNGWrite();
					}
				} else if(writeMovieFrames) {
					checkWriteRemoteFrame();
				}
				
				counter++;
			}
		}
	}
	
	static public void checkTimeStepLimits() {
		double maxActinDeltaT, maxMyoDeltaT;
		double myoK;
		if(MyosinMiniFilament.theMiniFilaments[0]!= null) {
			myoK = Myosin.myoStrokeSpringConstant*MyosinMiniFilament.theMiniFilaments[0].getNumBoundMyosins();
			maxMyoDeltaT = 0.5*MyosinMiniFilament.theMiniFilaments[0].myosinfilTransGamma/myoK;
			actualDeltaT = Math.min(deltaT,maxMyoDeltaT);
			for(int i = 0; i < Actin.actinCt; i++) {
				maxActinDeltaT = 0.5*Actin.theActins[i].actinfilTransGamma/myoK;
				actualDeltaT = Math.min(actualDeltaT,maxActinDeltaT);
			}
		}
		else {
			actualDeltaT = deltaT;
		}
	}

	
		
	static void glidingFilamentSetup () {
		/** creation of myosin surface for gliding filament assay */
		//testSurface = new MyosinSurface(Constants.xDimension/2,Constants.yDimension/2,Constants.xDimension,Constants.yDimension,0.002);
		//testSurface = new MyosinSurface(xDimension/2,yDimension/2,xDimension,yDimension/16,0.008);

	}

		
	static public String getBasePath() {
		return basePath;
	}
	
	static public String getBaseName() {
		return baseName;
	}
	
	static public void setBaseName(String name) {
		baseName = name;
	}
	
	public void setMoviePath(String path) {
		File moviedir = new File(basePath + File.separator + path);
		moviePath = moviedir.getAbsolutePath();
		
	}
	
	public String getMoviePath() {
		return moviePath;
	}

				
	static private String getNextDirectory(String basePath)  {
		File f;
		String path;
		try {
				path = basePath + "_" + InetAddress.getLocalHost().getHostName();
		} catch(Exception e)  { System.out.println("unkownHost"); path = (new Long(System.currentTimeMillis())).toString();}
		
		int i = 1;
		while((f = new File(path + "_" + i)).exists())  {
			i++;
		}
		if(!f.mkdir())  {
			System.out.println("can't make the directory");
			return null;
		}
		return f.getAbsolutePath();
	}
	
	private String extractBaseNameFromParfile(String path) throws Exception  {
		File f = new File(path);
		String name = f.getName();
		name = name.substring(0,name.lastIndexOf(".par"));
		return name;
	}

	private void makeOutputPaths(String parFileName, String runmode) throws Exception {
		
		File outdir,moviedir,dotAngDir;
		
		if(runmode.equals(RUN)) {
			if(parFileName == null || parFileName.equals("")) {
				basePath = "output";
			}
			else {
				basePath = parFileName.substring(0,parFileName.lastIndexOf(".par"));
				outdir = new File(basePath);
				outdir.mkdir();
				basePath = basePath + File.separator + "output";
			}
			outdir = new File(basePath);
			int j=1;
			while (outdir.isDirectory()) {
				System.out.println (outdir.getName() + " exists.... keeping file name BUT changing directory name");
				outdir = new File(basePath + "." + String.valueOf(j));
				j++;
			}
			outdir.mkdir();
			basePath = outdir.getAbsolutePath();
			moviedir = new File(basePath + File.separator + "movie");
			moviedir.mkdir();
			moviePath = moviedir.getAbsolutePath();
			
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator)en.nextElement();
				eval.init(basePath,baseName);
				eval.setOutputPaths(0);
			}
		}
		else if(runmode.equals(ITERATE) || runmode.equals(CLIENT)) {
			basePath = getNextDirectory(parFileName.substring(0,parFileName.lastIndexOf(".par")));
			try {
				theIterator.setOutputStream(new AMOutputStream(new File(basePath + File.separator + "iteratorOut")));
			} catch(Exception e) {
				System.out.println("can't open iterator output stream for path = " + basePath + ". Aborting itreration...");
				throw(e);
			}
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator)en.nextElement();
				eval.init(basePath,baseName);
			}
		}
		else if(runmode.equals(GATHER)) {
			basePath = new File(parFileName).getParent();
			if(basePath == null) {
				basePath = System.getProperty("user.dir");
			}
		}
		else if(runmode.equals(SERVE)) {
			basePath = System.getProperty("user.dir");
		}
	}
	
	static public void setRemote(Boolean r) {
		remote = running = r;
	}

	static public void setWriteMovieFrames(Boolean r) {
		writeMovieFrames = r;
	}
		
	public void checkPNGWrite () {
		if (timeSinceLastFrame >= frameInterval | simulationTime == 0) {
			File moviedir = new File(moviePath);
			if(!moviedir.exists()) moviedir.mkdir();
			if (paintConnected) {
				Actin.traceAllClusters();
				Cluster.setNLargestClusters(numClustersToPaint);
				//MyosinHolder.startConnectedRecursion();
				try { Thread.sleep(100); } catch (InterruptedException e) { System.out.println ("error sleeping while paused"); }
			}
			theFrame.showAll(); // paint before save
			String fileName = moviePath + File.separator + movieFileName + String.valueOf(fileIdFormat.format(counter)+".png");
			PNGer.pngFromBufferedImage(theFrame.getImage(),fileName);
			
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator)en.nextElement();
				eval.showGraphics();
				eval.pngFromBufferedImage(eval.getMoviePath(),movieFileName + String.valueOf(fileIdFormat.format(counter)+".png"));
			}
						
			timeSinceLastFrame = 0;
		}
		timeSinceLastFrame+=actualDeltaT;
	}
	
	private void checkWriteRemoteFrame() {
		timeSinceLastFrame+=actualDeltaT;
		if (timeSinceLastFrame >= frameInterval | simulationTime == actualDeltaT) {
			File moviedir = new File(moviePath);
			if(!moviedir.exists()) moviedir.mkdir();
			String fileName = moviePath + File.separator + movieFileName + String.valueOf(fileIdFormat.format(counter)+".png");
			writeRemoteFrame(fileName);
			
			Evaluator eval;
			for(Enumeration en = theEvaluators.elements(); en.hasMoreElements(); ) {
				eval = (Evaluator)en.nextElement();
				eval.remoteWriteMovieFrameToPng(eval.getMoviePath(),movieFileName + String.valueOf(fileIdFormat.format(counter)+".png"));
			}
						
			timeSinceLastFrame = 0;
		}
	}
	
	public void writeRemoteFrame(String fileName) {
		if(remoteImageWriter == null) {
			remoteImageWriter = new RemoteImageWriter(WorldFrame.width,WorldFrame.height);
			remoteImageWriter.initialize();
		}
		remoteImageWriter.draw();
		PNGer.pngFromBufferedImage(remoteImageWriter.getImage(),fileName);
	}
	
	
	
	public static String getInfoString() {
		
		double avgLength=0, stdvLength=0, diff;
		double count=0;
		for (int i=0;i<Actin.actinCt;i++) {
			if(!Actin.theActins[i].removeMe){
				avgLength+=Actin.theActins[i].physicalLength;
				count++;
			}
		}
		avgLength=avgLength/count;
		for (int i=0;i<Actin.actinCt;i++) {
			if(!Actin.theActins[i].removeMe){
				diff=avgLength-Actin.theActins[i].physicalLength;
				stdvLength+=diff*diff;
			}
		}
		stdvLength = Math.sqrt(stdvLength/(count-1));
		String infoString = "";
		infoString += " time = " + String.valueOf(timeOutFormat.format(Sim2D.simulationTime)) + "  ";
		infoString += " aveFilLength: " + String.valueOf(aveLengthFormat.format(avgLength)) + "  ";
		infoString += " stdvFilLength: " + String.valueOf(aveLengthFormat.format(stdvLength)) + "  ";
		infoString += " monDens: " + String.valueOf(aveLengthFormat.format(Monomer.getMonomerDensity())) + "µM  ";
		infoString += " nFils: " + String.valueOf(Actin.actinCt) + " ";
		infoString += " nMFils: " + String.valueOf(MyosinMiniFilament.miniFilamentCt) + "  ";
		infoString += " nMyos: " + String.valueOf(Myosin.myosinCt) + "  ";
		infoString += " nXLs: " + String.valueOf(Crosslinker.crosslinkerCt);
		return infoString;
	}
	
	public static String getTimeString() {
		return " time = " + String.valueOf(timeOutFormat.format(Sim2D.simulationTime)) + "  ";
	}
		
	/************************ SIM2D PARAMETERS ********************************/
	/* Define static variables associated with each User-changable parameter. */
	
	/** Dimensions of the world (nm) */
	static public double xDimension;
	
	static public double yDimension;
	
	static public double zDimension;
	
	/** time-step in seconds. */
	static public double deltaT;
	
	static private double actualDeltaT;
		
	/** MAximum simulation time in seconds. */
	static public double runTime;
	
	/** if false, no graphics */
	static public boolean remote;
		
	/** Sim time between paint steps. */
	static double paintStep;
	
	/** If true, write movie frames at intervals of frameInterval.*/
	static public boolean writeMovieFrames;
	
	/** Sim time between movie frames. */
	static public double frameInterval;
	
	// for painting control
	static public boolean paintOn;
	
	static public int paintConnectedLevel = 10000;
		
	/** If true, use brownian forces. */
	static public boolean brownianOn;
		
	static public boolean paintConnected;
	
	static public int numClustersToPaint = 5;
	
	static public Color [] clusterColors = {Color.red,Color.orange,Color.white,Color.green,Color.cyan};
	
	static public double tauForceAverage = 0;
	
	static public double forceAverageFraction = 0.01;
	
	static public double getDeltaT() {
		return actualDeltaT;
	}
	
	private static void staticInit() {
		
		/* Create and register a ParameterSetListener that will update any static variables for
		 * this class whose values depend on changeable parameters.
		 */
		Parameters.addParameterSetListener (
			new ParameterSetListener() {
				public void parametersChanged() {
					maxUniqueXDist = xDimension/2;
					maxUniqueYDist = yDimension/2;
					Point2D.topLBound = new Point2D (-xDimension,0);
					Point2D.topRBound = new Point2D (2*xDimension,0);
					Point2D.botLBound = new Point2D (-xDimension,yDimension);
					Point2D.botRBound = new Point2D (2*xDimension,yDimension);
					Point2D.leftTBound = new Point2D (0,-yDimension);
					Point2D.leftBBound = new Point2D (0,2*yDimension);
					Point2D.rightTBound = new Point2D (xDimension,-yDimension);
					Point2D.rightBBound = new Point2D (xDimension,2*yDimension);
					if(tauForceAverage > 0) {
						forceAverageFraction = deltaT/tauForceAverage;
					}
					actualDeltaT = deltaT;
				}
			}
		);
		
		/* Now create Paraeter objects associated with each Parameter variable defined above */
		Parameters.addParameter(className, "remote",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					setRemote(p.getBooleanValue());
				}
			}
		);
		Parameters.addParameter(className, "brownianOn",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					brownianOn = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className, "paintConnected",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					paintConnected = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className, "numClustersToPaint",5,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					numClustersToPaint = (int)p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "paintOn",false,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					paintOn = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className, "writeMovieFrames",true,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					writeMovieFrames = p.getBooleanValue();
				}
			}
		);
		Parameters.addParameter(className, "xDimension",10000,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					xDimension = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "yDimension",10000,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					yDimension = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "zDimension",200,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					zDimension = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "TimeStep",0.001,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					deltaT = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "TotalTime",300,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					runTime = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "paintStep",0.01,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					paintStep = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "frameInterval",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					frameInterval = p.getDoubleValue();
				}
			}
		);
		Parameters.addParameter(className, "tauForceAverage",0.1,
			new ParameterListener() {
				public void parameterChanged(Parameter p) {
					tauForceAverage = p.getDoubleValue();
				}
			}
		);
	}
}
		
		
		
		
		
		
		
		
		
		
		


