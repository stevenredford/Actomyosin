package iterators;

import java.util.*;
import java.io.*;
import java.net.*;

import parameters.*;
import io.*;
;

public class TransectIterator extends Iterator {

	Parameter theParameter;
	
	int numSteps;
		
	public TransectIterator()  {}
	
	
	public void iterate() throws Exception  {
		Parameter p;
		for(Enumeration e = paramsToVary.elements(); e.hasMoreElements();) {
			p = (Parameter) e.nextElement();
			Parameters.setParameter(p.getName(),p.getRandomValue());
		}
		for(int i = 0; i < numSteps; i++) {
			Parameters.setParameter(theParameter.getName(),theParameter.setDiscreteLevel(i,numSteps));
			doOneIterate(i);
		}
		gatherFiles();
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
				throw(e);
			}
		}
		else if(tag.equals("Parameter")) {
			theParameter = new Parameter(Parameter.DOUBLE);
			theParameter.loadFromBlock(in);
			if(!Parameters.isValidParameterName(theParameter.getName())) {
				throw new Exception("TransectIterator.loadParameter(): got invalid Parameter name = " + theParameter.getName());
			}
			paramsToVary.add(theParameter);
		}
		else if(tag.equalsIgnoreCase("CompositeParameter")) {
			theParameter = new CompositeParameter();
			theParameter.loadFromBlock(in);
			Parameters.addParameter(theParameter);
			paramsToVary.add(theParameter);
		}
		else if(tag.equals("NumSamples")) {
			numSteps = in.nextInt();
		}
		else super.loadParameter(tag,in);
	}
	
	
	

				
		
	
	
}
