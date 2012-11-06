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
#pragma warning(disable: 4018) // possible loss of data

#define UNIX
//#define WINDOWS
//#define GraphGC

string dirname = "/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/supportingFiles/";
string wdirname = "C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\PeakFiles\\";
ifstream SAMFILE;

const int NUM_OF_THREADS=8;
const int BUFFERSIZE=100000;
const int ClusterPromoters=500; //Cluster Promoters of Isoforms that are 5 kb away from each other
const int AssociateInteractions=5000; //decide if an interaction comes from a feature or not
const int NumberofGenes=35000;
const int IsoformCount=10; // Maximum number of isoforms
const int HalfClusterDist = 750;
const int NotDistalEnough = 1000;
const int coreprom_upstream = 2500;
const int coreprom_downstream = 500;
int NumberofPeakFiles;

int NofInteractorBins; // Depends on the maximum allowed junction distance
int MIN_Insert_Size;
int MaxJunctionDistance;
int BinSize;
int NumberofBins;
double ExpressionThr=2.0;
int CellType; // 0:mES, 1:XEN, 2:TS
int padding=700; //For Sequence Capture Probes
int NumberofExperiments=1;


#include "DataStructs_ProcessSAMFiles.h"
#include "GetOptions_ProcessSAMFiles.h"
#include "AssociateProbeswithFeatures.h"
#include "GCContentNorm.h"
#include "RESitesCount.h"
#include "mappability.h"
#include "PromoterClass_ProcessSAMFiles.h"
#include "NegativeControls_ProcessSAMFiles.h"
#include "ProcessSAMFiles.h"
#include "GenerateFileNames.h"
//#include "DetectEnrichedBins.h"
#include "AnnotateProms_with_PeakFiles.h"

vector <MetaPeakMap> metaPeaks;

void JoinPeaks(vector < PeakClass >& ap, boost::unordered::unordered_set< int >& p1, int currentpeak, string chr ){

	boost::unordered::unordered_set< int, string >::const_iterator iter;
	int chrindex;
		string ifpeak; // will keep if the peak of interest present in other files
		for(iter = p1.begin(); iter != p1.end(); ++iter){
			for(int z = 0;z < currentpeak;++z)
				ifpeak.append("0 ");
			ifpeak.append("1 ");
			int key = (*iter);
			for(int j= currentpeak+1; j < ap.size();++j){ // starting from the current peak position
				for (int k = 0; k < ap[j].p.size();++k){
					if(chr == ap[j].p[k].maptochrname){// get the right chromosome name
						chrindex = k; 
						break;
					}
				}
				boost::unordered::unordered_set<int>::iterator it = ap[j].p[chrindex].peaks.find(key);
				if (ap[j].p[chrindex].peaks.find((key)) == ap[j].p[chrindex].peaks.end()){ // if the peak is not found
					ifpeak.append("0 ");
				}
				else{
					ifpeak.append("1 ");	
					ap[j].p[chrindex].peaks.erase(it); // Remove the repeated peak
				}
			}
			metaPeaks.back().metapeaks[(key)] = ifpeak; // Add the key and peak profile
			ifpeak.clear();
		} // went through all peaks
	}

/*
void ProcessPeaks(PromoterClass& Promoters, NegCtrlClass& NegativeControls,PeakClass& PC,DetectEnrichedBins& EnrichedBins,string &PeakFN, string abname, bool ifBED,string INTERACTIONFILENAMEBASE,int CellType,int abindex,int ExperimentIndex){
	ifstream PeakFile;
	PeakFile.open(PeakFN.c_str());
	if (abname.substr(0,2) == "HS")
		ifBED = 1;
	else 
		ifBED = 0;
	PC.ReadPeakFile(PeakFile,ifBED);
	PeakFile.close();

	string FileName1, FileName2,FileName3,FileName4,FileName5;
	PC.AnnotatewithNegCtrls(NegativeControls,abindex);
	PC.AnnotatewithPromoters(Promoters,abindex);
	EnrichedBins.AssociatePeaksWithIntBins(Promoters,INTERACTIONFILENAMEBASE,abname,CellType,abindex,ExperimentIndex);
	EnrichedBins.AssociatePeaksWithIntBins_NegCtrls(NegativeControls,INTERACTIONFILENAMEBASE,abname,abindex,ExperimentIndex);

}
*/
int main (int argc,char* argv[]){

	RESitesClass dpnIIsites;
	dpnIIsites.InitialiseVars();

string SAMFILENAME, INTERACTIONFILENAMEBASE, BaseFileName;
////////////////////////////////////////////////////////
//READ AND PROCESS PEAKS
string PeakFileName;
int abindex=0;
string ext1,ext2,FileName1,FileName2,FileName3,FileName4,fname;
ifstream AbNameFile;
AbNameFile.open("AbNames.txt");
vector < PeakClass > AllPeaks;
vector < string > AbNames;
vector <string> ChrNames;
do{
	AbNameFile >> fname >> ext1;
	cout << fname << "  " << ext1 << endl;
	if(fname.compare("END")==0)
		break;
#ifdef UNIX
	PeakFileName.append(dirname);
#endif
#ifdef WINDOWS
	PeakFileName.append(wdirname);
#endif
	ifstream PeakFile;
	PeakFileName.append(fname);
	bool ifBED;
	AbNames.push_back(ext1);
	PeakClass Peaks;
	//	ProcessPeaks(Promoters,NegativeControls,Peaks,EnrichedBins,PeakFileName,ext1,0,INTERACTIONFILENAMEBASE,CellType,abindex,0);
	PeakFile.open(PeakFileName.c_str());
	if(!PeakFile.is_open())
		cerr << PeakFileName << "   cannot be opened"  << endl;
	if (ext1.substr(0,2) == "HS")
		ifBED = 1;
	else 
		ifBED = 0;
	Peaks.ReadPeakFile(PeakFile, ifBED,dpnIIsites,ChrNames);
	PeakFile.close();
	
	AllPeaks.push_back(Peaks);
//	Peaks.~PeakClass();
	AllPeaks[abindex].abnames.push_back(ext1);
	++abindex;
	PeakFileName.clear();
	fname.clear();
	ext1.clear();
	cout << AbNames[abindex-1] << "   read " << endl;
}while(fname!="END");
cout << INTERACTIONFILENAMEBASE << "          All Peak Files Read" << endl;
NumberofPeakFiles = abindex;

string chr;
for( int j = 0; j< ChrNames.size();++j){ // For every chromosome 
	chr = ChrNames[j];
	metaPeaks.push_back(MetaPeakMap());
	metaPeaks.back().maptochrname = chr;
	for (int k = 0; k<AllPeaks.size(); ++k){
		for (int l = 0; l<AllPeaks[k].p.size(); ++l){
			if(AllPeaks[k].p[l].maptochrname == chr){
				JoinPeaks(AllPeaks, AllPeaks[k].p[l].peaks,k,chr);			
			}
		}
	}
}

cout << "Meta Peak Struct Generated " << endl;

boost::unordered::unordered_map< int, string >::const_iterator iter;

ofstream outf3("MM9.MetaPeaks_All_13_PeakFiles.txt");

outf3 << "chrname" << '\t' << "ClosestRESite" << '\t';
for(int m = 0; m < AbNames.size();++m)
	outf3 << AbNames[m] << '\t';
outf3 << endl;
for (int i = 0; i < metaPeaks.size();++i ){
	for(iter = metaPeaks[i].metapeaks.begin(); iter != metaPeaks[i].metapeaks.end(); ++iter){
		outf3 << metaPeaks[i].maptochrname <<  '\t' << iter->first << '\t' << iter->second << endl;
	}
}



/*
#ifdef UNIX

	if (argc < 4) {
		print_usage();
		return -1;
	}
	SAMFILENAME.append(argv[1]);
	INTERACTIONFILENAMEBASE.append(argv[2]);
	CellType=atoi(argv[3]);
	cout << "INPUTFileName    " <<  SAMFILENAME << endl;
	cout << "INTERACTORFILENAMEBASE" << INTERACTIONFILENAMEBASE << endl;
	cout << "Cell Type        " <<  CellType << endl;
	
	SAMFILE.open(SAMFILENAME.c_str());
	if(!SAMFILE.is_open())
		cerr << "Input File cannot be opened" << endl;
#endif
#ifdef WINDOWS
	SAMFILE.open("F:\\3C\\HiC_mES_10082012\\HiC_mES_10082012_truncated.noduplpairs.sam");
	INTERACTIONFILENAMEBASE.append("test_");
	BinSize=5000; MaxJunctionDistance=100000;
	NofInteractorBins=(MaxJunctionDistance)/BinSize;
	NumberofBins=NofInteractorBins*2;
#endif


//   --        INITIALISE CLASSES   --
	PromoterClass Promoters;
	NegCtrlClass NegativeControls;
	ProcessSAM samfile;
	ProbeSet mm9probes;
	GCContent mm9GC;
	MappabilityClass mapp;
//	DetectEnrichedBins EnrichedBins;
//-------------------//------------------------------
string ext1,ext2,FileName1,FileName2,FileName3,FileName4;
	dpnIIsites.GetRESitesCount("chr7", 52830546, 52831546);
//	mapp.InitialiseVars();
	Promoters.InitialiseData();
	Promoters.ReadPromoterAnnotation(dpnIIsites);
	cout << "Promoters Annotated" << endl;
	NegativeControls.InitialiseData();
	NegativeControls.FillNegativeCtrls(dpnIIsites);
	cout << "Negative Controls Annotated" << endl;
	mm9probes.InitialiseData();
//	mm9probes.ReadProbeCoordinates();
	cout << "Probes Read" << endl;
/////////////////////////////////////////////////////////

	samfile.ProcessTheSAMFile(Promoters,NegativeControls,mm9probes,dpnIIsites,mapp);
	Promoters.PrintInteractions();
	*/
}
