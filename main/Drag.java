package main;


public class Drag {
	public static int pValueIndex=0;
	public static int parIndex=1;
	public static int perpIndex=2;
	public static int rotIndex=3;
	
	public static float[] cylinderPerpendicularRotationSigmas={0.229f,
																-0.216f,
																-0.260f,
																-0.303f,
																-0.348f,
																-0.392f,
																-0.436f,
																-0.481f,
																-0.526f,
																-0.571f,
																-0.616f,
																-0.662f};
	public static float[] cylinderParallelTranslationSigmas={0.70f,
															0.25f,
															0.21f,
															0.18f,
															0.14f,
															0.07f,
															0.03f,
															-0.02f,
															-0.05f,
															-0.11f,
															-0.15f,
															-0.20f};
	public static float[] cylinderPerpendicularTranslationSigmas={1.14f,
																0.99f,
																0.97f,
																0.96f,
																0.94f,
																0.91f,
																0.90f,
																0.88f,
																0.87f,
																0.86f,
																0.85f,
																0.84f};
	public static float[] cylinderPValues={0.99f,//we approximated this sigmas for this pValue
											0.50f,
											0.45f,
											0.40f,
											0.35f,
											0.30f,
											0.25f,
											0.20f,
											0.15f,
											0.10f,
											0.05f,
											0.0f};
	
	//Has to have the order specified by the array index variables pValueIndex,parIndex,perpIndex,rotIndex;
	public static float[][] cylinderDragValues={
		{cylinderPValues[0],cylinderParallelTranslationSigmas[0],cylinderPerpendicularTranslationSigmas[0],cylinderPerpendicularRotationSigmas[0]},
		{cylinderPValues[1],cylinderParallelTranslationSigmas[1],cylinderPerpendicularTranslationSigmas[1],cylinderPerpendicularRotationSigmas[1]},
		{cylinderPValues[2],cylinderParallelTranslationSigmas[2],cylinderPerpendicularTranslationSigmas[2],cylinderPerpendicularRotationSigmas[2]},
		{cylinderPValues[3],cylinderParallelTranslationSigmas[3],cylinderPerpendicularTranslationSigmas[3],cylinderPerpendicularRotationSigmas[3]},
		{cylinderPValues[4],cylinderParallelTranslationSigmas[4],cylinderPerpendicularTranslationSigmas[4],cylinderPerpendicularRotationSigmas[4]},
		{cylinderPValues[5],cylinderParallelTranslationSigmas[5],cylinderPerpendicularTranslationSigmas[5],cylinderPerpendicularRotationSigmas[5]},
		{cylinderPValues[6],cylinderParallelTranslationSigmas[6],cylinderPerpendicularTranslationSigmas[6],cylinderPerpendicularRotationSigmas[6]},
		{cylinderPValues[7],cylinderParallelTranslationSigmas[7],cylinderPerpendicularTranslationSigmas[7],cylinderPerpendicularRotationSigmas[7]},
		{cylinderPValues[8],cylinderParallelTranslationSigmas[8],cylinderPerpendicularTranslationSigmas[8],cylinderPerpendicularRotationSigmas[8]},
		{cylinderPValues[9],cylinderParallelTranslationSigmas[9],cylinderPerpendicularTranslationSigmas[9],cylinderPerpendicularRotationSigmas[9]},
		{cylinderPValues[10],cylinderParallelTranslationSigmas[10],cylinderPerpendicularTranslationSigmas[10],cylinderPerpendicularRotationSigmas[10]},
		{cylinderPValues[11],cylinderParallelTranslationSigmas[11],cylinderPerpendicularTranslationSigmas[11],cylinderPerpendicularRotationSigmas[11]}
	};
}
