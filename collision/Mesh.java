package collision;


import java.awt.*;

import main.*;
import gui.*;
import io.*;
import util.*;


public class Mesh {
	public static Mesh ACTIN_MESH;
	public static Mesh BARRIER_MESH;
	public static Mesh MFIL_MESH;
	public static Mesh MVTAIL_MESH;
	public static final float SHIFT=0.49f*0;
	
	public static float MIN_LENGTH_FOR_LINE_ALGORITHM=200;
	public static float SIZE=200;
	public static float X_BIN_WIDTH=SIZE;
	public static float Y_BIN_WIDTH=SIZE;
	public static float Z_BIN_WIDTH=SIZE;
	public static int BIN_DEPTH=400;
	public static boolean OVERLAP=true;
	
	public double[][][] meshpoints;
	public int [][] timeStamps;
	public int [][] activeCts;
	
	public static float xBinWidth=X_BIN_WIDTH;
	public static float yBinWidth=Y_BIN_WIDTH;
	public static float zBinWidth=Z_BIN_WIDTH;
	
	public static float nXbinsFractional=(float)(Sim2D.xDimension)/xBinWidth;
	public static float nYbinsFractional=(float)(Sim2D.yDimension)/yBinWidth;
	public static int nXBins=1+(int)Math.ceil(nXbinsFractional);
	public static int nYBins=1+(int)Math.ceil(nYbinsFractional);
	public static int binDepth=BIN_DEPTH;
	
	public static float xMinValue=0;
	public static float xMaxValue=(float)Sim2D.xDimension;
	
	public static float yMinValue=0;
	public static float yMaxValue=(float)Sim2D.yDimension;
	
	public boolean showGraphics=false;
	MeshCellGraphics[][] cells;
	

	public static void create(){
		FileOps.report(
				" "+"\n"
				+"Mesh stats:"+"\n"
				+"nXBins = "+nXBins+";"+" xBinWidth = "+xBinWidth+" nm"+"\n"
				+"nYBins = "+nYBins+";"+" yBinWidth = "+yBinWidth+" nm"+"\n"
				+"binDepth = "+binDepth+"\n"
				+" "+"\n"
		);
		
		ACTIN_MESH=new Mesh(false);
		BARRIER_MESH=new Mesh(false);
		MFIL_MESH= new Mesh(false);
		MVTAIL_MESH= new Mesh(false);
	}
	
		
	// Constructor for the mesh class
	public Mesh(boolean showGraphics){
		this.showGraphics = showGraphics;
		meshpoints = new double[nXBins][nYBins][binDepth];
		timeStamps = new int[nXBins][nYBins];
		activeCts = new int[nXBins][nYBins];
		reset();
		makeGraphics();
	}
	
	//puts -1 in all meshpoints array elements and 0 in timeStamps and activeCts array elements
	public void reset() {
		for(int x=0;x<meshpoints.length;x++){
			for(int y=0;y<meshpoints[x].length;y++){
				activeCts[x][y] = 0;
				timeStamps[x][y] = 0;
				for(int i=0;i<meshpoints[x][y].length;i++){
					meshpoints[x][y][i]=-1;
				}
			}
		}
	}
		
	
	public void makeGraphics(){
		if(showGraphics){
			cells=new MeshCellGraphics[nXBins][nYBins];
			for(int x=0;x<cells.length;x++){
				for(int y=0;y<cells[x].length;y++){
					double cx=xMinValue+xBinWidth/2+x*xBinWidth;
					double cy=yMinValue+yBinWidth/2+y*yBinWidth;
					cells[x][y]=new MeshCellGraphics(cx,cy,xBinWidth,yBinWidth);
					cells[x][y].xNumber=x;
					cells[x][y].yNumber=y;
				}
			}
		}
	}
	
	/** This method adds a line segment to the Mesh.  Assumes that the line segment
	 * does NOT cross periodic boundaries.  If it does, you're in trouble.  If the line segment
	 * is longer than the value of the static variable MIN_LENGTH_FOR_LINE_ALGORITHM
	 * then Bresenhelm's algorithm is used, otherwise, just fills in the smallest
	 * rectangle containing both endpoints.
	 */
	public void addLineSegmentToMesh(int id, Point2D startPt, Point2D stopPt){

		double segLength = Math.sqrt(Math.pow(startPt.x-stopPt.x,2) + Math.pow(startPt.y-stopPt.y,2));
		boolean useLineAlgorithm=segLength>MIN_LENGTH_FOR_LINE_ALGORITHM;
		
		//X AXIS
		int startBinX = getBinX(startPt.x);
		int stopBinX = getBinX(stopPt.x);
		
		if(!useLineAlgorithm&&startBinX>stopBinX){
			int temp=stopBinX;
			stopBinX=startBinX;
			startBinX=temp;
		}
		
		//Y AXIS
		int startBinY = getBinY(startPt.y);
		int stopBinY = getBinY(stopPt.y);
		
		if(!useLineAlgorithm&&startBinY>stopBinY){
			int temp=stopBinY;
			stopBinY=startBinY;
			startBinY=temp;
		}
		
		if(OVERLAP){
			startBinX--;
			stopBinX++;
			if(startBinX<0)startBinX=0;
			if(stopBinX>=nXBins)stopBinX=nXBins-1;
			
			startBinY--;
			stopBinY++;
			if(startBinY<0)startBinY=0;
			if(stopBinY>=nYBins)stopBinY=nYBins-1;
		}
	
		//FILL BINS
		if(useLineAlgorithm){
			lineBresenham(startBinX,startBinY,stopBinX,stopBinY,id);
		} else {
			for(int x=startBinX;x<=stopBinX;x++){
				for(int y=startBinY;y<=stopBinY;y++){
					addToMesh(x,y,id);
				}
			}
		}
	}


	/** Fills in the smallest rectangle containing the circle with the given center and radius.
	 * This method respects periodic boundary conditions.
	 */
	public void addPointToMesh(int id, Point2D center, double radius){
			
		//X AXIS
		int startBinX = getBinX(center.x - radius);
		int stopBinX = getBinX(center.x + radius);
		
		
		//Y AXIS
		int startBinY = getBinY(center.y - radius);
		int stopBinY = getBinY(center.y + radius);
		
		
		if(OVERLAP){
			startBinX--;
			stopBinX++;
			startBinY--;
			stopBinY++;
		}
		
		for(int x=startBinX;x<=stopBinX;x++){
			for(int y=startBinY;y<=stopBinY;y++){
				addToMesh(wrapX(x),wrapY(y),id);
			}
		}
	}

	public void lineBresenham(int x0, int y0, int x1, int y1,int id){
        int dy = y1 - y0;
        int dx = x1 - x0;
        int stepx, stepy;

        if (dy < 0) { dy = -dy;  stepy = -1; } else { stepy = 1; }
        if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }
        dy <<= 1;                                                  // dy is now 2*dy
        dx <<= 1;                                                  // dx is now 2*dx

        fillMeshCell(id,x0,y0);
        if (dx > dy) {
            int fraction = dy - (dx >> 1);                         // same as 2*dy - dx
            while (x0 != x1) {
                if (fraction >= 0) {
                    y0 += stepy;
                    fraction -= dx;                                // same as fraction -= 2*dx
                }
                x0 += stepx;
                fraction += dy;                                    // same as fraction -= 2*dy
                fillMeshCell(id,x0,y0);
            }
        } else {
            int fraction = dx - (dy >> 1);
            while (y0 != y1) {
                if (fraction >= 0) {
                    x0 += stepx;
                    fraction -= dy;
                }
                y0 += stepy;
                fraction += dx;
                fillMeshCell(id,x0,y0);
            }
        }
	}

	public void fillMeshCell(int id,int x0,int y0){
		for(int x=x0-1;x<=x0+1;x++){
			for(int y=y0-1;y<=y0+1;y++){
				addToMesh(wrapX(x),wrapY(y),id);
			}
		}
	}
	
	public void addToMesh (int x, int y, int id) {
		int timeToSet = Sim2D.counter+1;
		if(timeStamps[x][y]!=timeToSet){
			timeStamps[x][y]=timeToSet;
			meshpoints[x][y][0]=id;
			activeCts[x][y]=1;
		} else {
			meshpoints[x][y][activeCts[x][y]]=id;
			activeCts[x][y]++;
			if (activeCts[x][y] >= Mesh.BIN_DEPTH) {
				activeCts[x][y]=Mesh.BIN_DEPTH-1;
			}
		}
	}
	
	/** Returns an integer corresponding to the X bin position of val.  Pays
	 * no attention to periodic boundary conditions.
	 */
	public static int getBinX (double val) {
		return(int)(((val-xMinValue)/xBinWidth)+SHIFT);
	}
	
	/** Returns an integer corresponding to the Y bin position of val.  Pays
	 * no attention to periodic boundary conditions.
	 */
	public static int getBinY (double val) {
		return (int)(((val-yMinValue)/yBinWidth)+SHIFT);
	}
		
			
	public static int wrapX(int x) {
		return x < 0 ? nXBins+x : x >= nXBins ? x - nXBins : x;
	}
		
	public static int wrapY(int y) {
		return y < 0 ? nYBins+y : y >= nYBins ? y - nYBins : y;
	}
		

	public static void reportMeshActives () {
		FileOps.reportln ("Actin_Mesh actives = " + ACTIN_MESH.countActiveEntries());
		FileOps.reportln ("Barrier_Mesh actives = " + BARRIER_MESH.countActiveEntries());
	}
	
	public int countActiveEntries () {
		int numActive = 0;
		for (int x=0;x<nXBins;x++) {
			for (int y=0;y<nYBins;y++) {
				if (timeStamps[x][y] == Sim2D.counter+1) {
					numActive += activeCts[x][y];
				}
			}
		}
		return numActive;
	}
	
	public void print(int x,int y,int z){
		FileOps.reportln("ENV.SIM_TIME="+Sim2D.simulationTime+"**************************************************************");
		
		FileOps.reportln("["+x+"]"+"["+y+"]"+"["+z+"]"+"\n"
				+"\t"+"[TIMESTAMP]"+" = "+timeStamps[x][y]
				+"  [COUNT]"+" = "+activeCts[x][y]);

		int idStop=activeCts[x][y];
		if(idStop>=meshpoints[x][y].length)idStop=meshpoints[x][y].length;
		for(int id=0;id<idStop;id++){
			FileOps.reportln("\t"+"i="+id+"["+x+"]"+"["+y+"]"+"["+z+"]"+"[ID]"+" = "+meshpoints[x][y][id]);
		}
		FileOps.reportln("**************************************************************");
	}
	
	public void print(){
		FileOps.reportln("ENV.SIM_TIME="+Sim2D.simulationTime+"**************************************************************");
		for(int x=0;x<meshpoints.length;x++){
			for(int y=0;y<meshpoints[x].length;y++){
				for(int z=0;z<meshpoints[x][y].length;z++){
					FileOps.reportln("["+x+"]"+"["+y+"]"+"["+z+"]"+"\n"
							+"\t"+"[TIMESTAMP]"+" = "+timeStamps[x][y]
							+"  [COUNT]"+" = "+activeCts[x][y]
					);
					int idStop=activeCts[x][y];
					//idStop=meshpoints[x][y][z].length;
					for(int id=0;id<idStop;id++){
						FileOps.reportln("\t"+"i="+id+"["+x+"]"+"["+y+"]"+"["+z+"]"+"[ID]"+" = "+meshpoints[x][y][id]);
					}
				}
			}
		}
		FileOps.reportln("**************************************************************");
	}
	
	public void drawYourself (Graphics g, double scale, double [] offset) {
		for(int x=0;x<cells.length;x++){
			for(int y=0;y<cells[x].length;y++) {
				if(timeStamps[x][y] == Sim2D.counter+1) {
					cells[x][y].setCount(activeCts[x][y]);
				}
				else cells[x][y].setCount(0);
				cells[x][y].drawYourself(g,scale,offset);
			}
		}
	}

	
	static public void staticInit() {
		
		Sim2D.addProxy(new AgentProxy() {
			public void reset() {
				ACTIN_MESH.reset();
				BARRIER_MESH.reset();
				MFIL_MESH.reset();
				MVTAIL_MESH.reset();
			}
		});
	}
	
	
}
