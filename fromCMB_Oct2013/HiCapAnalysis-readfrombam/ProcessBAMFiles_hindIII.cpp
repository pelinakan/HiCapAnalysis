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

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

#ifdef _CHAR16T //redefined definition problem
#define CHAR16_T
#endif

#define BOOST_DISABLE_ASSERTS
#define NDEBUG

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

#include "api/api_global.h"
#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/SamConstants.h"
#include "api/SamProgramChain.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"

#define UNIX
//#define WINDOWS
//#define GraphGC
#define CMB

#ifdef CMB
string dirname =  "/mnt/crick/pelina/bin/supportingFiles/";
#endif

string wdirname = "C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\PeakFiles\\";
string ExpFileName;

int MinimumJunctionDistance; // To be entered by the user
int MinNumberofReads; // To be entered by the user
int CellType; // 0:mES, 1:XEN, 2:TS // read from the experiments file
int TotalNumberofPairs = 0;
const double SignificanceThreshold = 0.001; // Whether an interaction is significantly above background
const int BUFFERSIZE = 40000000;
const int ClusterPromoters = 1000; //Cluster Promoters of Isoforms that are 500 bp away from each other
const int coreprom_upstream = 3000;
const int coreprom_downstream = 1000;
const int AssociateInteractions = 1500; //decide if an interaction come from a feature (used for negative controls)
const int BinSize  = 1000; // Only To Calculate Background Interaction Frequencies
const int NumberofGenes = 35000;
const int IsoformCount = 10; // Maximum number of isoforms
const int MaxInsertLen = 700; // To exclude pairs associated with very distal RE sites
const int NOFEXPERIMENTS = 2; // Number of Experiments
int padding = 700; //For Sequence Capture Probes
double ExpressionThr = 2.0;
double Mappability_Threshold = 0.5;
string whichchr;

#include "linear.h"
#include "DataStructs_ProcessSAMFiles.h"
#include "GetOptions_ProcessSAMFiles.h"
#include "RESitesCount_MemAccess_HindIII.h"
//#include "RESitesCount_MemAccess.h"
#include "Mappability.h"
#include "PromoterClass_ProcessSAMFiles.h"
#include "NegativeControls_ProcessSAMFiles.h"
#include "AnnotateProms_with_PeakFiles.h"
#include "MetaPeakCollection.h"
#include "AssociateProbeswithFeatures.h"
#include "ProcessBAMFiles.h"
#include "BackgroundInteractionFrequency.h"
#include "DetectEnrichedBins.h"

int main (int argc,char* argv[]){
/*
	RESitesClass dpnIIs;
	dpnIIs.InitialiseVars();
	int *renums;
	renums = new int [2];
	dpnIIs.GettheREPositions("chr1",197195432,renums);
	return 0;
*/
	string BaseFileName;
#ifdef UNIX

	if (argc < 6) {
		print_usage();
		return -1;
	}
	ExpFileName = argv[1];
	cout << ExpFileName << endl;
	MinNumberofReads = atoi(argv[2]);
	MinNumberofReads = double(MinNumberofReads);
	MinimumJunctionDistance = atoi(argv[3]);
	BaseFileName = argv[4];
	whichchr = argv[5];

	cout << "Min Number of Pairs      " << MinNumberofReads << endl;
	cout << "Min Junction Distance    " << MinimumJunctionDistance << endl;
#endif
#ifdef WINDOWS
	ExpFileName = "Experiments.txt";
	if (argc < 3) {
		print_usage();
		return -1;
	}
	MinNumberofReads = atoi(argv[1]);
	MinNumberofReads = double(MinNumberofReads);
	MinimumJunctionDistance = atoi(argv[2]);

	cout << "Min Number of Pairs      " << MinNumberofReads << endl;
	cout << "Min Junction Distance    " << MinimumJunctionDistance << endl;
#endif

//   --        INITIALISE CLASSES   --
	PromoterClass Promoters;
	NegCtrlClass NegativeControls;
	ProbeSet mm9probes;
	DetermineBackgroundLevels background;
	ProcessBAM bamfile;
	DetectEnrichedBins EnrichedBins;

	MappabilityClass mapp;
	RESitesClass dpnIIsites;
	dpnIIsites.InitialiseVars();
	mapp.InitialiseVars();
//-------------------//------------------------------
	Promoters.InitialiseData();
	Promoters.ReadPromoterAnnotation(dpnIIsites, mapp);
	cout << "Promoters Annotated" << endl;
	NegativeControls.InitialiseData();
	NegativeControls.FillNegativeCtrls(dpnIIsites, mapp, whichchr);
	cout << "Negative Controls Annotated" << endl;

	ReadMetaPeakFile(); // If peaks are already processed.

#ifdef WINDOWS
	ifstream ExpFile("Experiments.txt"); // Experiment Names and details
#endif
#ifdef UNIX
	ifstream ExpFile(ExpFileName.c_str());
#endif
	string BAMFILENAME, ExperimentName;
	vector < string > ExperimentNames;
	int ExperimentNo = 0;
	ExpFile >> BAMFILENAME >> ExperimentName >> CellType;
	do{ // Reads all the pairs in each experiment and fills the interaction maps
		ExperimentNames.push_back(ExperimentName);
		cout << BAMFILENAME << "     will be read" << endl;
		bamfile.ProcessTheBAMFile(Promoters,NegativeControls,mm9probes,BAMFILENAME, ExperimentNo, whichchr);
		
		cout << "Detecting Interactions";
		background.CalculateMeanandStdRegress(NegativeControls, ExperimentName, ExperimentNo); 
		cout << ". ";
		
		EnrichedBins.DetectEnrichedInteractionBins(Promoters,background,ExperimentName,ExperimentNo);
		cout << ". ";
		EnrichedBins.AssociatePeaksWithIntBins(Promoters,ExperimentName,dpnIIsites);
		cout << ". ";

		EnrichedBins.DetectEnrichedInteractionBins_NegCtrls(NegativeControls,background,ExperimentName,ExperimentNo);
		cout << ". ";

		EnrichedBins.AssociatePeaksWithIntBins_NegCtrls(NegativeControls,ExperimentName,dpnIIsites);
		cout << "finished" << endl;
		
		cout << "Cleaning the background" << endl;
		background.bglevels.mean_downstream.clear();
		background.bglevels.mean_upstream.clear();
		background.bglevels.stdev_downstream.clear();
		background.bglevels.stdev_upstream.clear();
		++ExperimentNo;
		cout << BAMFILENAME << "     finished" << endl;
		ExpFile >> BAMFILENAME >> ExperimentName >> CellType;
	}while(BAMFILENAME!="END");

	int NofofExperiments = ExperimentNo;

	EnrichedBins.PrintAllInteractions(Promoters,NegativeControls,BaseFileName, NofofExperiments, ExperimentNames); //In one file

}
