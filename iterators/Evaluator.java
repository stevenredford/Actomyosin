package iterators;


import java.util.*;
import java.io.*;
import java.awt.*;

import main.*;
import io.*;
import util.*;


public class Evaluator  {
	
	String name = null;
	
	double  stopTime;
	
	double interval = 0.1;
	
	double threshold;

	boolean autoPass = true;
	
	static String sepString = ";";

	public double nextEval = 0;

	String basePath = null;
	
	String baseName = null;
	
	String dataPath= null;
	
	String moviePath = null;
	
	File movieDir = null;
	
	File baseDir = null;
	
	boolean firstEval = true;
	
	String dataFileName = null;

	
	PrintWriter dataPW;
	
	
	int evaluationCt;
	
	
	public static class EvaluatorGraphic {
		/** Graphics for Evaluators */
		public static EvaluatorGraphic [] evaluatorGraphics = new EvaluatorGraphic [100];
		public static int evalGraphicCt = 0;
		Point2D pt1 = new Point2D();
		Point2D pt2 = new Point2D();
		Color eGColor;
		
		public EvaluatorGraphic (Point2D point1, Point2D point2, Color color) {
			pt1.copy(point1);
			pt2.copy(point2);
			eGColor = color;
			addEvaluatorGraphic();
		}
		
		public void addEvaluatorGraphic () {
			evaluatorGraphics[evalGraphicCt] = this;
			evalGraphicCt++;
		}
		
		public static void showAllEvaluatorGraphics (Graphics g) {
			for (int i=0;i<evalGraphicCt;i++) {
				System.out.println("drawing evalG " + i);
				g.setColor(evaluatorGraphics[i].eGColor);
				drawLineTwixtPoints(evaluatorGraphics[i].pt1, evaluatorGraphics[i].pt2, g);
			}
		}
		
		public static void drawLineTwixtPoints (Point2D p1, Point2D p2, Graphics g) {
			g.drawLine(p1.getPixX(), p1.getPixY(), p2.getPixX(), p2.getPixY());
		}
		
	}
	
	
	/** This method get's called once by the Iterator at the beginning of an
	 * iterator sequence (i.e. a sequence of simulation runs for different
	 * choices of parameter value */
	public void init(String path, String name) throws Exception {
		basePath = path + File.separator + getName();
		baseDir = new File(basePath);
		baseDir.mkdir();
		baseName = name;
	}
	
	public String getName() {
		if(name == null) {
			name = getClass().getName();
			name = name.substring(name.lastIndexOf(".")+1);
			return name;
		}
		return name;
	}
	
	final public void load(AMInputStream in)  throws Exception {
		
		String tag = in.nextTag();
		while(!tag.equals("endEvaluator"))  {
			loadParameter(tag, in);
			tag = in.nextTag();
		}
	}

	/** Override this method to load in your evaluator specific parameters.  The basic format
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
	 * else super.loadParameter(tag,in);

	 *
	 *
	 * The last call to super.loadParameter() ensures that parameters defined in the
	 * Evaluator class also get loaded...
	 *  */
	public void loadParameter(String tag, AMInputStream in)  throws Exception {
			
		if(tag.equalsIgnoreCase("StopTime"))  {
			stopTime = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("Interval"))  {
			interval = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("Threshold"))  {
			threshold = in.nextDouble();
		}
		else if(tag.equalsIgnoreCase("AutoPass"))  {
			autoPass = in.nextBoolean();
		}
		else if(tag.equalsIgnoreCase("Name"))  {
			name = in.nextString();
		}
		else throw new Exception("Evaluator.loadParameter(): got bad tag = " + tag);
	}
	
	public void remoteWriteMovieFrameToPng(String path,String name) {}
	
	public void showGraphics() {}
	
	public void pngFromBufferedImage(String path, String name) {}

	public void resetGraphics(boolean remote) {}
	
	public String [] getGraphicsInfoStrings() { return null; }
	
	public void mkDataFile () {}

	/** Subclasses override to actually write file header */
	public void writeDataFileHeader() {
		
	}
	
	/** This gets called just before a new iteration begins.  Override to reset any
	 * of your Evaluator-specific data,. Again, make sure your overriding method calls super.reset()
	 * */
	public void reset()  {
		nextEval = 0;
		firstEval = true;
		evaluationCt = 0;
	}
	
	/** Override this method to determine when a specific simulation should stop.  The default
	 * behavior encoded by the base class is to stop when tm == stopTime.
	 */
	public boolean stop(double tm)  {
		if(tm > nextEval) {
			if (firstEval){
				doFirstEval();
				mkDataFile();
				firstEval = false;
			} else {
				evaluate(tm);
				evaluationCt++;
				nextEval += interval;
			}
		}
		return decide(tm);
	}
	
	/** Give evaluators a chance to initialize something the first time through */
	public void doFirstEval() {}
		
	
	public void evaluate(double tm) {}
	
	public boolean decide(double tm) {
		if(tm > stopTime) {
			if(dataPW != null) dataPW.close();
			return true;
		}
		return false;
	}
	
	/** Override this method to implement specific decisions about when to "pass" a simulation run acording to
	 * some criterion that your Evaluator evaluates.  The default behavior is to consult the value
	 * of the parameter autoPass, whose default value is true, and whose value can be specified in the input file.
	 * */
	public boolean gotHit()  {
		return autoPass;
	}
	
	public void recordHit(AMOutputStream out, int cnt) {
			
		out.printTag("Evaluator",name);
		out.increaseIndent();
		out.printTag("Data",new String("\"" + getDataString() + "\""));
		if(hasData()) {
			out.printTag("DataPath",new String("\"" + getDataPath() + "\""));
		}
		if(Sim2D.writeMovieFrames && hasMovie()) {
			out.printTag("MoviePath",new String("\"" + moviePath + "\""));
			String path = new String(basePath + File.separator + "ff_" + cnt + ".png");
			out.printTag("FinalFrame",new String("\"" + path + "\""));
			remoteWriteMovieFrameToPng(basePath,new String("ff_" + cnt + ".png"));
		}
		out.decreaseIndent();
		out.printTag("endEvaluator");
	}


	
	/** Override this method to supply a list of tab-delimited header Strings that define column headers for the data
	 * you supply via the getDataString() method. The number of header Strings should be identical to number of data values.
	 * and the orders should correspond.  Basic format is:
	 *
	 * 		return new String(name + "\t" + "header1" + "\t" + "header2" + "\t" + .... + "headerN");
	 **/
	public String getHeaderString() {
		return new String("");
	}
	
	/** Override this method to supply a list of tab-delimited data values.
	 * Basic format is:
	 *
	 * 		return new String(name + "\t" + "value1" + "\t" + "value2" + ... + "\t" + "valueN");

	 **/
	public String getDataString() {
		return new String("");
	}
	
	/** Override this method to specify if Evaluator subclass writes out a datafile containing
	 * subtype-specific data.
	 * */
	public boolean hasData()  {
		return false;
	}
	
	/** Override this method to specify if Evaluator subclass writes out movie frames containing
	 * subtype-specific data.
	 * */
	public boolean hasMovie()  {
		return false;
	}
	
	public void setOutputPaths(int cnt) throws Exception {
		String path = new String("data" + cnt);
		File dataDir = new File(basePath + File.separator + path);
		if(!dataDir.mkdir()) {
			throw new Exception("Evaluator.setOutputPaths: Can't make new directory for path = " + basePath + File.separator + path);
		}
		dataPath = dataDir.getAbsolutePath();
		
		if(hasMovie()) {
			path = new String("movie" + cnt);
			movieDir = new File(basePath + File.separator + path);
			movieDir.mkdir();
			moviePath = movieDir.getAbsolutePath();
		}

	}

	
	public final String getDataPath() {
		return dataPath;
	}
	
	public final String getMoviePath() {
		return moviePath;
	}
	
	final public void setBasePath(String path) {
		basePath = path;
	}
		
	final public void setBaseName(String name) {
		baseName = name;
	}
	
	static void printDivider (PrintWriter pw) {
		pw.print(sepString);
	}
	
	static void printDivider (PrintWriter pw,String sep) {
		pw.print(sep);
	}
	
	static void printValue (int val, PrintWriter pw) {
		pw.print(String.valueOf(val));
		printDivider(pw);
	}
	
	static void printValue (double val, PrintWriter pw) {
		pw.print(String.valueOf(val));
		printDivider(pw);
	}

	static void printValue (int val, PrintWriter pw,String sep) {
		pw.print(String.valueOf(val));
		printDivider(pw,sep);
	}
	
	static void printValue (double val, PrintWriter pw,String sep) {
		pw.print(String.valueOf(val));
		printDivider(pw,sep);
	}
		

}
