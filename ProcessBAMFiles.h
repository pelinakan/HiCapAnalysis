class ProcessBAM{
	friend class PromoterClass;
	friend class NegCtrlClass;
	friend class ProbeSet;
public:
	void ProcessTheBAMFile(PromoterClass&,NegCtrlClass&,ProbeSet&,string,int);
private:
	void GetPairInformation(BamTools::BamAlignment&, BamTools::BamAlignment&, PairStruct&);
	void ProcessPairs(PairStruct&,RESitesClass&,PromoterClass&,NegCtrlClass&,int, int&,int&, int&);
	boost::unordered::unordered_map<int, string> RefIDtoChrNames;
};

void ProcessBAM::GetPairInformation(BamTools::BamAlignment& al, BamTools::BamAlignment& almate, PairStruct& temppair){
	boost::unordered::unordered_map<int, string>::const_iterator it1 = RefIDtoChrNames.find(al.RefID);
	if(it1 != RefIDtoChrNames.end())
		temppair.chr_1 = it1->second;

	boost::unordered::unordered_map<int, string>::const_iterator it2 = RefIDtoChrNames.find(almate.RefID);
	if(it2 != RefIDtoChrNames.end())
		temppair.chr_2 = it2->second;
	
	temppair.startcoord = al.Position;
	temppair.endcoord = almate.Position;
	
//	cout << it1->second << "   " << it2->second << "   " << al.Position << "   " << almate.Position << endl;
}
void ProcessBAM::ProcessTheBAMFile(PromoterClass& promoters,NegCtrlClass& negctrls,ProbeSet& mm9prs, string BAMFILENAME,int ExperimentNo){

int NofPairsAnnWProms = 0;
int NofPairsAnnWNegCtrls = 0;
int NofPairsNoAnn = 0;

BamReader reader;
if ( !reader.Open(BAMFILENAME.c_str()) ) {
	cerr << "Could not open input BAM files." << endl;
}
// retrieve 'metadata' from BAM files, these are required by BamWriter
const SamHeader header = reader.GetHeader();
const RefVector references = reader.GetReferenceData();

// Make a map of chr names to RefIDs
RefVector::const_iterator chrit;
for (chrit = references.begin(); chrit != references.end(); ++chrit){
	int key = reader.GetReferenceID(chrit->RefName);
	RefIDtoChrNames[key] = chrit->RefName;
//	cout << chrit->RefName << '\t' << key << endl;
}

PairStruct* PairPool = new PairStruct[BUFFERSIZE];
BamAlignment al,almate;
vector < RESitesClass > dpnIIparallel; // Each thread will have its own RESiteClass

int nofthreads, poolsize = 0;

	nofthreads = 1;
	cout << "Number of Threads   " << nofthreads << endl;
	for (int t = 0; t < nofthreads; ++t){
		cout << "Thread No " << t << "    ";
		dpnIIparallel.push_back(RESitesClass());
		dpnIIparallel.back().InitialiseVars(); // Initialize all RESiteClass
	}
//Read the BAM file into PairPool vector
	cout << "Reading BAM file..." << endl;
while ( reader.GetNextAlignmentCore(al) ){
		reader.GetNextAlignmentCore(almate);
		GetPairInformation(al,almate,PairPool[poolsize]);
		if(al.Position == -1 || almate.Position == -1){
			cout << "Mate not correct, exiting  the reading loop" << endl;
			break;
		}
		++poolsize;

	if(poolsize == BUFFERSIZE){
		cout << poolsize << "    Pairs Read" << endl << "Annotating Interactions... ";
		int tid = 0;
		for(int i=0;i<poolsize;++i){
			ProcessPairs(PairPool[i],dpnIIparallel[tid],promoters,negctrls,ExperimentNo, NofPairsAnnWProms, NofPairsAnnWNegCtrls, NofPairsNoAnn); //Annotate the pair
		}
		cout << endl << BAMFILENAME << "   " <<  poolsize << "   pairs finished " << endl;
		delete[] PairPool;
		PairPool = new PairStruct[BUFFERSIZE];
		poolsize = 0;
	}
}
if(poolsize > 0){ // This is to read the last batch of pairs
	cout << poolsize << "    Pairs Read" << endl << "Annotating Interactions...";

	nofthreads = 1;
	int tid = 0;
	for(int i = 0; i < poolsize; ++i){
		ProcessPairs(PairPool[i],dpnIIparallel[tid],promoters,negctrls,ExperimentNo,NofPairsAnnWProms,NofPairsAnnWNegCtrls,NofPairsNoAnn); //Annotate the pair
	}
	cout << endl << BAMFILENAME << "   reading finished" << endl;
	delete[] PairPool;
	PairPool = new PairStruct[BUFFERSIZE];
	poolsize = 0;
}
cout << "Number of Pairs Annotated with Promoters    " << NofPairsAnnWProms << endl;
cout << "Number of Pairs Annotated with NegCtrls     " << NofPairsAnnWNegCtrls << endl;
cout << "Number of Pairs with No Annotation           " << NofPairsNoAnn << endl;

}
void ProcessBAM::ProcessPairs(PairStruct& thispair, RESitesClass& dpnII, PromoterClass& promoters, NegCtrlClass& negctrls,int ExperimentNo, int& annwp, int& annwnc, int& noann){

	bool pairann = false;
	int *renums1, *renums2;
	int resite_firstinpair, resite_secondinpair;
	bool passed1 = 1, passed2 = 1;
	renums1 = new int [2];	renums2 = new int [2];
	
	bool re1found = dpnII.GettheREPositions(thispair.chr_1,thispair.startcoord,renums1);
	bool re2found = dpnII.GettheREPositions(thispair.chr_2,thispair.endcoord,renums2);
	if(re1found)	
		resite_firstinpair = renums1[0];
	else
		resite_firstinpair = thispair.startcoord;
	if(re2found)
		resite_secondinpair = renums2[0];
	else
		resite_secondinpair = thispair.endcoord;

	pairann = promoters.AnnotatewithPromoters(thispair.chr_1,resite_firstinpair,thispair.startcoord,thispair.chr_2,resite_secondinpair,thispair.endcoord,dpnII,ExperimentNo);
	if(pairann)
		++annwp;
	else {// If the pair is not annotated with promoters
		pairann = negctrls.AnnotateWithNegCtrls(thispair.chr_1,resite_firstinpair,thispair.startcoord,thispair.chr_2,resite_secondinpair,thispair.endcoord,dpnII,ExperimentNo);
		if(pairann)
			++annwnc;
		else
			++noann;
	}
/*
	if(passed1 && passed2)
		thispair.flagged = 0;  // Regular Pair
	else // if the REsites are very far away from the read coordinate
		thispair.flagged = 1; // RE sites are far away
	if(!(re1found && re2found))	
		thispair.flagged = 2; // RE sites cannot be found
*/
}