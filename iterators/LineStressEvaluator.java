package iterators;


import java.util.*;
import java.awt.Color;
import java.io.*;

import collision.CollisionDetector;

import main.*;
import util.*;
import io.*;
import parameters.*;

public class LineStressEvaluator extends Evaluator {
	
	LineSegment referenceSegment;
	LineSegment [] lineSegments;
	boolean showLineSegmentGraphics = false;
	Color forceLineColor = Color.orange;
	boolean useSegmentArray = false;
	Point2D arrayIncrement;
	int numberOfSegments;

	
	int maxAs = 1000;
	Actin [] crossingAs = new Actin[maxAs];
	double [] crossingAsArc = new double[maxAs];
	int crossingAsCt = 0;
	
	int maxLs = 20000;
	Crosslinker [] crossingLs = new Crosslinker[maxLs];
	int crossingLsCt = 0;
	
	int maxMs = 1000;
	MyosinMiniFilament [] crossingMs = new MyosinMiniFilament[maxMs];
	int [] crossingMsInd = new int[maxMs];
	int crossingMsCt = 0;
	
	
	public void init(String path, String name) throws Exception  {

		super.init(path,name);
		
		Parameters.setParameter("monitorStress",true);
		
		if(useSegmentArray) {
			
			lineSegments = new LineSegment[numberOfSegments];
			Point2D p1 = new Point2D(referenceSegment.p1);
			Point2D p2 = new Point2D(referenceSegment.p2);
			LineSegment seg;
			for(int i = 0; i < numberOfSegments; i++) {
				seg = new LineSegment(p1,p2);
				seg.init();
				lineSegments[i] = seg;
				p1.inc(arrayIncrement);
				p2.inc(arrayIncrement);
				p1.wrapPoint();
				p2.wrapPoint();
			}
		}
	}
	
	
	public  void setOutputPaths(int cnt) throws Exception {
		super.setOutputPaths(cnt);
		dataFileName = dataPath + File.separator + "LineStresses.dat";
	}

	public void mkDataFile () {
		try {
			dataPW = new PrintWriter(new FileWriter(new File (dataFileName)),true);
			writeDataFileHeader();
			
		} catch (IOException ioe) { System.out.println("LineStressEvaluator.mkDataFile(): error creating File"); }
	}

	public void loadParameter(String tag, AMInputStream in)  throws Exception {
		if(tag.equalsIgnoreCase("LineSegment")) {
			loadLineSegment(in);
		}
		else if(tag.equalsIgnoreCase("LineSegmentArray")) {
			loadLineSegmentArray(in);
		}
		else if(tag.equals("showLineSegmentGraphics")) {
			showLineSegmentGraphics = in.nextBoolean();
		}
		else super.loadParameter(tag,in);
	}
	
	private void loadLineSegment(AMInputStream in) throws Exception  {
		String tag = in.nextTag();
		double x1 = -1, y1 = -1, x2 = -1, y2 = -1;
		while(!tag.equals("endLineSegment")) {
			if(tag.equals("X1"))  {
				x1 = in.nextDouble();
			}
			else if(tag.equals("Y1"))  {
				y1 = in.nextDouble();
			}
			else if(tag.equals("X2"))  {
				x2 = in.nextDouble();
			}
			else if(tag.equals("Y2"))  {
				y2 = in.nextDouble();
			}
		}
		if(x1 > 0 && x2>0 && y1>0 && y2>0) {
			referenceSegment = new LineSegment(x1,y1,x2,y2);
		}
		else throw new Exception("LineStressEvaluator.loadLineSegment(): Incomplete specification of egment endpoints");
		
	}
	
	private void loadLineSegmentArray(AMInputStream in) throws Exception  {
		String tag = in.nextTag();
		numberOfSegments = -1;
		double x1 = -1, y1 = -1, x2 = -1, y2 = -1;
		arrayIncrement = new Point2D(-1,-1);
		while(!tag.equals("endLineSegmentArray")) {
			if(tag.equals("X1"))  {
				x1 = in.nextDouble();
			}
			else if(tag.equals("Y1"))  {
				y1 = in.nextDouble();
			}
			else if(tag.equals("X2"))  {
				x2 = in.nextDouble();
			}
			else if(tag.equals("Y2"))  {
				y2 = in.nextDouble();
			}
			else if(tag.equalsIgnoreCase("xIncrement"))  {
				arrayIncrement.x = in.nextDouble();
			}
			else if(tag.equals("yIncrement"))  {
				arrayIncrement.y = in.nextDouble();
			}
			else if(tag.equals("NumberOfSegments"))  {
				numberOfSegments = in.nextInt();
			}
			tag = in.nextTag();
		}
		if(x1 < 0 || y1 < 0 || x2 < 0 || y2  < 0 ||
		   arrayIncrement.x < 0 || arrayIncrement.y < 0 || numberOfSegments < 0) {
			throw new Exception("LineStressEvaluator.loadLineSegmentArray(): incmplete array specification");
		}
		referenceSegment = new LineSegment(x1,y1,x2,y2);
		useSegmentArray = true;
	}
	
	

	public void evaluate(double tm)  {
		writeStressData();
	}
	
	public void writeDataFileHeader()  {
		dataPW.print("posMyoStress" + "\t" + "negMyoStress" + "\t" + "totalMyoStress" + "\t");
		dataPW.print("posActinStress" + "\t" + "negActinStress" + "\t" + "totalActinStress" + "\t");
		dataPW.print("posLinkerStress" + "\t" + "negLinkerStress" + "\t" + "totalLinkerStress" + "\t");
		dataPW.print("totalStress" + "\n");
	}
		
	private void writeStressData()  {

		double [] stresses = new double[10];
		if(useSegmentArray) {

			double [] s;
			
			for(int i = 0; i < numberOfSegments; i++) {
				lineSegments[i].measureStress();
				s = lineSegments[i].getStresses();
				for(int j = 0; j < 10; j++) {
					stresses[j] += s[j];
				}
			}
			for(int j = 0; j < 10; j++) {
				stresses[j]/= numberOfSegments;
			}
		}
		
		else {
			referenceSegment.measureStress();
			stresses = referenceSegment.getStresses();
		}
		for(int j = 0; j < 10; j++) {
			dataPW.print(stresses[j] + "\t");
		}
		dataPW.print("\n");
		dataPW.flush();
	}
	
	
	
	
	class LineSegment  {
		
		Point2D p1 = new Point2D();
		
		Point2D p2 = new Point2D();
		
		Point2D uVect = new Point2D();
		
		public LineSegment(double x1, double y1, double x2, double y2) {
			p1.x = x1; p1.y = y1; p2.x = x2; p2.y = y2;
		}
		
		public LineSegment(Point2D p1, Point2D p2) {
			this.p1.x = p1.x; this.p1.y = p1.y; this.p2.x = p2.x; this.p2.y = p2.y;
		}
		
		double totalStress = 0,totalMyoStress = 0,totalLinkerStress = 0,totalActinStress = 0;
		double posMyoStress = 0, posLinkerStress = 0, posActinStress = 0;
		double negMyoStress = 0, negActinStress = 0, negLinkerStress = 0;
		double maxMyoStress, maxActinStress, maxLinkerStress;
		
		public void init() {
			// initialize the unit vector defined from p1 to p2
			uVect.sub(p2, p1);
			uVect.uVect();
			uVect.orthogonalize();
			
			// make graphic to show this force line
			if(showLineSegmentGraphics)
				new EvaluatorGraphic(p1,p2,forceLineColor);
		}
		
		public void measureStress() {
			posMyoStress = negMyoStress = totalMyoStress = posActinStress  = negActinStress = 0;
			totalActinStress = posLinkerStress = negLinkerStress = totalLinkerStress = totalStress = 0;
			findCrossingActins();
			findCrossingLinkers();
			findCrossingMinifilaments();
			getOrthogonalForcesAtCrossings();
		}
		
		public double [] getStresses() {
			double [] s = new double[10];
			s[0] = posMyoStress;
			s[1] = negMyoStress;
			s[2] = totalMyoStress;
			s[3] = posActinStress;
			s[4] = negActinStress;
			s[5] = totalActinStress;
			s[6] = posLinkerStress;
			s[7] = negLinkerStress;
			s[8] = totalLinkerStress;
			s[9] = totalStress;
			return s;
		}
		
	
		/** This method sets the value for forceSum as the sum of all orthogonal forces across the specified line segment */
		public void getOrthogonalForcesAtCrossings () {
			Monomer curM;
			double stress;
			double orthoComponent;
			for (int i=0;i<crossingAsCt;i++) {
				curM = crossingAs[i].getMonomerAt(crossingAsArc[i]);
				if(curM == null) {
					System.out.println("curM null");
				}
				else if(curM.myFilament == null) {
					System.out.println("curMmyFilament null");
				}

				else if(curM.myFilament.uVect == null) {
					System.out.println("curMmyFilament.uVect null");
				}

				else if(uVect == null) {
					System.out.println("uVect null");
				}

				orthoComponent = Math.abs(Point2D.dot(curM.myFilament.uVect,uVect));
				stress = curM.avgMonStress*orthoComponent;
				if(stress > 0)
					posActinStress += stress;
				else negActinStress += stress;
				if(stress > maxActinStress) maxActinStress = stress;
				totalActinStress += stress;
			}
			for (int i=0;i<crossingLsCt;i++) {
				orthoComponent = Math.abs(Point2D.dot(crossingLs[i].uVect,uVect));
				stress = crossingLs[i].forceAv*orthoComponent;
				if(stress > 0)
					posLinkerStress += stress;
				else negLinkerStress += stress;
				if(stress > maxLinkerStress) {
					maxLinkerStress = stress;
				}
				totalLinkerStress += stress;
			}
			for (int i=0;i<crossingMsCt;i++) {
				orthoComponent = Math.abs(Point2D.dot(crossingMs[i].uVect,uVect));
				stress = crossingMs[i].getInternalStress(crossingMsInd[i])*orthoComponent;
				if(stress > 0)
					posMyoStress += stress;
				else negMyoStress += stress;
				if(stress > maxMyoStress) maxMyoStress = stress;
				totalMyoStress += stress;
			}
			totalStress = totalMyoStress + totalLinkerStress + totalActinStress;
		}
		
		public double actinCrossed (Actin curA) {
			return CollisionDetector.actinSegCollision(curA, p1, p2, new Point2D());
		}
	
		public void addCrossedA (Actin a, double arcL) {
			crossingAs[crossingAsCt] = a;
			crossingAsArc[crossingAsCt] = (arcL);
			crossingAsCt++;
		}
	
		public void findCrossingActins () {
			crossingAsCt = 0;	// reset
			double crossingArc;
			for (int i=0;i<Actin.actinCt;i++) {
				crossingArc = actinCrossed(Actin.theActins[i]);
				if (crossingArc != -1) {
					addCrossedA (Actin.theActins[i],crossingArc);
				}
			}
		}
		
		public void findCrossingLinkers () {
			crossingLsCt = 0;	// reset
			double crossingArc;
			for (int i=0;i<Crosslinker.crosslinkerCt;i++) {
				if (linkerCrossed(Crosslinker.theCrosslinkers[i])) {
					addCrossedL(Crosslinker.theCrosslinkers[i]);
				}
			}
		}
		
		public boolean linkerCrossed (Crosslinker l) {
			if(!l.isBound())
				return false;
			
			return (CollisionDetector.lineSegmentIntersectionTest(l.site1,l.getSite2IgnorePeriod(),p1,p2,new Point2D()) >= 0);
		}
	
		public void addCrossedL (Crosslinker l) {
			crossingLs[crossingLsCt] = l;
			crossingLsCt++;
		}
	
		public void findCrossingMinifilaments () {
			crossingMsCt = 0;	// reset
			double crossingArc;
			for (int i=0;i<MyosinMiniFilament.miniFilamentCt;i++) {
				minifilamentCrossed(MyosinMiniFilament.theMiniFilaments[i]);
			}
		}
		
		public void minifilamentCrossed (MyosinMiniFilament m) {
			Point2D ep1 = m.getEndPoint1();
			Point2D ep2 = m.getEndPoint2();
			Point2D [][] segs = Point2D.segmentLine(ep1,ep2,m.uVect);
			int segCt = segs.length;
			double [] segLength = new double[segCt];
			for(int i = 0; i < segCt; i++) {
				segLength[i] = Point2D.getDistance(segs[i][0], segs[i][1]);
			}
			
			Point2D intPt = new Point2D();
			double arcL = -1;
			double arcLNormed = 0;
	
			switch (segCt) {
			case 1:
				arcLNormed = CollisionDetector.lineSegmentIntersectionTest(segs[0][0],segs[0][1],p1,p2,intPt);
				if (arcLNormed >=0) {
					arcL = arcLNormed*segLength[0];
					crossingMs[crossingMsCt] = m;
					crossingMsInd[crossingMsCt] = m.getStressIndex(arcL);
					crossingMsCt++;
				}
				return;
			case 2:
				arcLNormed = CollisionDetector.lineSegmentIntersectionTest(segs[0][0],segs[0][1],p1,p2,intPt);
				if (arcLNormed >=0) {
					arcL = arcLNormed*segLength[0];
					crossingMs[crossingMsCt] = m;
					crossingMsInd[crossingMsCt] = m.getStressIndex(arcL);
					crossingMsCt++;
					return;
				}
				arcLNormed = CollisionDetector.lineSegmentIntersectionTest(segs[1][0],segs[1][1],p1,p2,intPt);
				if (arcLNormed >=0) {
					arcL = segLength[0] + arcLNormed*segLength[1];
					crossingMs[crossingMsCt] = m;
					crossingMsInd[crossingMsCt] = m.getStressIndex(arcL);
					crossingMsCt++;
					return;
				}
				return;
			case 3:
				arcLNormed = CollisionDetector.lineSegmentIntersectionTest(segs[0][0],segs[0][1],p1,p2,intPt);
				if (arcLNormed >=0) {
					arcL = arcLNormed*segLength[0];
					crossingMs[crossingMsCt] = m;
					crossingMsInd[crossingMsCt] = m.getStressIndex(arcL);
					crossingMsCt++;
					return;
				}
				arcLNormed = CollisionDetector.lineSegmentIntersectionTest(segs[1][0],segs[1][1],p1,p2,intPt);
				if (arcLNormed >=0) {
					arcL = segLength[0] + arcLNormed*segLength[1];
					crossingMs[crossingMsCt] = m;
					crossingMsInd[crossingMsCt] = m.getStressIndex(arcL);
					crossingMsCt++;
					return;
				}
				arcLNormed = CollisionDetector.lineSegmentIntersectionTest(segs[2][0],segs[2][1],p1,p2,intPt);
				if (arcLNormed >=0) {
					arcL = segLength[0] + segLength[1] + arcLNormed*segLength[2];
					crossingMs[crossingMsCt] = m;
					crossingMsInd[crossingMsCt] = m.getStressIndex(arcL);
					crossingMsCt++;
					return;
				}
				return;
			}
			return;
		}
	}
}
		



