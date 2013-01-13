#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cstring>
#include <string>
#include <time.h>
#include <vector>
#include <omp.h>
#include <sstream>
//#include "/bubo/home/h20/pelin/3Cproj/bin/boost_1_50_0/boost/unordered/unordered_map.hpp"
//#include "/bubo/home/h20/pelin/3Cproj/bin/boost_1_50_0/boost/unordered/unordered_set.hpp"
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
using namespace std;

#ifdef _CHAR16T //redefined definition problem
#define CHAR16_T
#endif

<<<<<<< HEAD
#define BOOST_DISABLE_ASSERTS
#define NDEBUG

=======
>>>>>>> 95a9108d132a134505bfacc1093538dd278d0ab0
// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
#include "alglibmisc.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "dataanalysis.h"
#include "specialfunctions.h"
#include "solvers.h"
#include "optimization.h"
#include "diffequations.h"
#include "fasttransforms.h"
#include "integration.h"
#include "interpolation.h"
using namespace alglib_impl;
#pragma warning(disable: 4018) // possible loss of data
#pragma warning(disable: 4244) // possible loss of data
#pragma warning(disable: 4715) // not all control paths return a value

#define UNIX
//#define WINDOWS
//#define GraphGC
string dirname = "/bubo/home/h20/pelin/3Cproj/bin/supportingFiles/";
string wdirname = "C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\PeakFiles\\";
<<<<<<< HEAD
string ExpFileName;
=======

>>>>>>> 95a9108d132a134505bfacc1093538dd278d0ab0
int MinimumJunctionDistance; // To be entered by the user
int MinNumberofReads; // To be entered by the user
int CellType; // 0:mES, 1:XEN, 2:TS // read from the experiments file
int TotalNumberofPairs = 0;
const int BUFFERSIZE = 4000000;
const int ClusterPromoters = 1000; //Cluster Promoters of Isoforms that are 500 bp away from each other
const int coreprom_upstream = 2500;
const int coreprom_downstream = 500;
const double SignificanceThreshold = 0.001;
const int AssociateInteractions = 1500; //decide if an interaction come from a feature (used for negativecontrols)
int padding = 700; //For Sequence Capture Probes
double ExpressionThr = 2.0;
const 	int BinSize  = 1000; // Only To Calculate Background Interaction Frequencies
const int NumberofGenes = 35000;
const int IsoformCount = 10; // Maximum number of isoforms

#include "linear.h"
#include "DataStructs_ProcessSAMFiles.h"
#include "GetOptions_ProcessSAMFiles.h"
#include "GenerateFileNames.h"
#include "GCContentNorm.h"
#include "RESitesCount_MemAccess.h"
#include "Mappability.h"
#include "PromoterClass_ProcessSAMFiles.h"
#include "NegativeControls_ProcessSAMFiles.h"
#include "AnnotateProms_with_PeakFiles.h"
#include "MetaPeakCollection.h"
#include "AssociateProbeswithFeatures.h"
#include "ProcessSAMFiles.h"
#include "BackgroundInteractionFrequency.h"
#include "DetectEnrichedBins.h"
#include "CompareExperiments.h"
#include "EnhancerAnalysis.h"
#include "LaminB1.h"

int main (int argc,char* argv[]){

#ifdef UNIX

<<<<<<< HEAD
	if (argc < 4) {
		print_usage();
		return -1;
	}
	ExpFileName = argv[1];
	cout << ExpFileName << endl;
	MinNumberofReads = atoi(argv[2]);
	MinNumberofReads = double(MinNumberofReads);
	MinimumJunctionDistance = atoi(argv[3]);

	cout << "Min Number of Pairs      " << MinNumberofReads << endl;
	cout << "Min Junction Distance    " << MinimumJunctionDistance << endl;
#endif
#ifdef WINDOWS
	ExpFileName = "Experiments.txt";
=======
	if (argc < 3) {
		print_usage();
		return -1;
	}
	MinNumberofReads = atoi(argv[1]);
	MinNumberofReads = double(MinNumberofReads);
	MinimumJunctionDistance = atoi(argv[2]);

	cout << "Min Number of Pairs      " << MinNumberofReads << endl;
	cout << "Min Junction Distance    " << MinimumJunctionDistance << endl;
>>>>>>> 95a9108d132a134505bfacc1093538dd278d0ab0
#endif
/////////////////////////////////////////////////////////One-time use functions
	//	CompareExperiments intersectExp;
	//	intersectExp.readandCompare();
<<<<<<< HEAD

	//	EnhancerPositions();

//	LaminB1Class LaminB1;
=======

	//	EnhancerPositions();

	LaminB1Class LaminB1;
>>>>>>> 95a9108d132a134505bfacc1093538dd278d0ab0
/////////////////////////////////////////////////////////

//   --        INITIALISE CLASSES   --
	PromoterClass Promoters;
	NegCtrlClass NegativeControls;
	ProbeSet mm9probes;
	DetermineBackgroundLevels background;
	ProcessSAM samfile;
	DetectEnrichedBins EnrichedBins;

	MappabilityClass mapp;
	RESitesClass dpnIIsites;
	GCContent mm9GC;
	dpnIIsites.InitialiseVars();
	mapp.InitialiseVars();
//-------------------//------------------------------
	Promoters.InitialiseData();
	Promoters.ReadPromoterAnnotation(dpnIIsites, mapp);
	cout << "Promoters Annotated" << endl;
	NegativeControls.InitialiseData();
	NegativeControls.FillNegativeCtrls(dpnIIsites, mapp);
	cout << "Negative Controls Annotated" << endl;

/////////////////////////////////////////////////////////One time use functions
<<<<<<< HEAD
//	LaminB1.CalculateGeneLaminValues(Promoters); 

/////////////////////////////////////////////////////////

	ReadMetaPeakFile(); // If peaks are already processed.
	ifstream ExpFile(ExpFileName.c_str()); // Experiment Names and details
	string SAMFILENAME, BaseFileName;
	ExpFile >> SAMFILENAME >> BaseFileName >> CellType;
	while (SAMFILENAME != ("END")){ // Reads all the pairs in each experiment and fills the interaction maps
		cout << SAMFILENAME << "     will be read" << endl;
		samfile.ProcessTheSAMFile(Promoters,NegativeControls,mm9probes,SAMFILENAME);
		ExpFile >> SAMFILENAME >> BaseFileName >> CellType;
	}
	cout << "ALL SAMFILES READ " << endl;
=======
	LaminB1.CalculateGeneLaminValues(Promoters); 

/////////////////////////////////////////////////////////

	return 1;

	ReadMetaPeakFile(); // If peaks are already processed.
	ifstream ExpFile("Experiments.txt"); // Experiment Names and details
	string SAMFILENAME, BaseFileName;
	do{ // Reads all the pairs in each experiment and fills the interaction maps
		ExpFile >> SAMFILENAME >> BaseFileName >> CellType;
		cout << SAMFILENAME << "     will be read" << endl;
		samfile.ProcessTheSAMFile(Promoters,NegativeControls,mm9probes,SAMFILENAME);
	}while(!(ExpFile.eof()));
	cout << "SAMFILES READ " << endl;
>>>>>>> 95a9108d132a134505bfacc1093538dd278d0ab0

	background.CalculateMeanandStdRegress(NegativeControls, BaseFileName); 
	cout << " Background Levels Calculated " << endl;

	EnrichedBins.DetectEnrichedInteractionBins(Promoters,background,BaseFileName);
	EnrichedBins.AssociatePeaksWithIntBins(Promoters,BaseFileName,dpnIIsites);

	EnrichedBins.DetectEnrichedInteractionBins_NegCtrls(NegativeControls,background,BaseFileName);
	EnrichedBins.AssociatePeaksWithIntBins_NegCtrls(NegativeControls,BaseFileName,dpnIIsites);

	EnrichedBins.DetectSignicantFeatFeatInteractions(Promoters, NegativeControls, BaseFileName);

	string FileName;
	FileName.append(BaseFileName);
	FileName.append("Summary.txt");
	ofstream outf(FileName.c_str());
	outf << "Base File Name" << '\t' << BaseFileName << endl;
	outf << "Number of promoters" << '\t' << Promoters.NofPromoters << endl;
	outf << "Number of genes" << '\t' << Promoters.NofGenes << endl;
	outf << "Number of pairs processed" << '\t' << TotalNumberofPairs << endl;
	outf << "Minimum number of pair required to call an interaction" << '\t' << MinNumberofReads << endl;
	outf << "Minimum Junction Distance" << '\t' << MinimumJunctionDistance << endl;

	EnrichedBins.SummariseInteractions(Promoters,NegativeControls,BaseFileName,0,outf);

	dpnIIsites.CleanClass();

}
