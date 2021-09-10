package main;

abstract public class AgentProxy {
	
	public void checkCollisions(double simulationTime) {}
	public void registerForCollisionDetection() {}
	public void doForces(double dT) {}
	public void step(double dT) {}
	public void checkBarrierCrossings() {}
	public void init() {}
	public void initDiagnostics() {}
	public void reset() {}
	
		
}
