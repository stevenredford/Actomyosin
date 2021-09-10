package main;
import java.awt.*;
import java.awt.geom.GeneralPath;
import java.text.DecimalFormat;
import java.util.*;
import java.awt.Graphics2D;
import sun.tools.tree.ThisExpression;

import analysis.*;
import util.*;
import collision.*;
import io.*;
import iterators.GlidingAssayEvaluator;
import gui.*;
import parameters.*;

public class MiscellaneousParameters {
	/**used for making mixed polarity bundles
	 * 
	 */
	public static double fractionBarbedEndsRight;
	public static int nBarbedEndsRight;
	public static int nBarbedEndsLeft;
	public static double actinLength;
	public static double bundleLength;
	public static double avgFilsXsection;
	public static double bundleWidth;
	public static double bundleSpacing;
	public static double actinPairSpacing;
	public static int barbedEndLocationIndicator;//integer 1-10 indicating different ways of arranging 3 barbed ends
	static void staticInit() {
		Parameters.addParameter("MiscellaneousParameters","fractionBarbedEndsRight",0.5,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						fractionBarbedEndsRight= p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","nBarbedEndsRight",0.5,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						nBarbedEndsRight=(int) (p.getDoubleValue());
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","nBarbedEndsLeft",0.5,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						nBarbedEndsLeft=(int) (p.getDoubleValue());
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","actinLength",1000,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						actinLength= p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","bundleLength",10000,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						bundleLength= p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","bundleWidth",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						bundleWidth= p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","bundleSpacing",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						bundleSpacing= p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","barbedEndLocationIndicator",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						barbedEndLocationIndicator= (int) p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","actinPairSpacing",0,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						actinPairSpacing= p.getDoubleValue();
					}
				}
			);
		Parameters.addParameter("MiscellaneousParameters","avgFilsXsection",10000,
				new ParameterListener() {
					public void parameterChanged(Parameter p) {
						avgFilsXsection= p.getDoubleValue();
					}
				}
			);
	}
}
