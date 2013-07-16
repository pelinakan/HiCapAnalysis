class ProcessBAM{
	friend class PromoterClass;
	friend class NegCtrlClass;
	friend class ProbeSet;
public:
	void ProcessTheBAMFile(PromoterClass&,NegCtrlClass&,ProbeSet&,string,int);
private:
	void GetPairInformation(BamTools::BamAlignment&, PairStruct&);
	void ProcessPairs(PairStruct&,RESitesClass&,PromoterClass&,NegCtrlClass&,int);
	boost::unordered::unordered_map<int, string> RefIDtoChrNames;
};

void ProcessBAM::GetPairInformation(BamTools::BamAlignment& al, PairStruct& temppair){
	boost::unordered::unordered_map<int, string>::const_iterator it1 = RefIDtoChrNames.find(al.RefID);
	if(it1 != RefIDtoChrNames.end())
		temppair.chr_1 = it1->second;

	boost::unordered::unordered_map<int, string>::const_iterator it2 = RefIDtoChrNames.find(al.MateRefID);
	if(it2 != RefIDtoChrNames.end())
		temppair.chr_2 = it2->second;
	
	temppair.startcoord = al.Position;
	temppair.endcoord = al.MatePosition;

//	cout << it1->second << "  " << it2->second << "   " << al.Position << "   " << al.MatePosition << endl;
//	cout << temppair.chr_1 << "   " << temppair.chr_2 << "    " << temppair.startcoord << "    " << temppair.endcoord << endl;

}
void ProcessBAM::ProcessTheBAMFile(PromoterClass& promoters,NegCtrlClass& negctrls,ProbeSet& mm9prs, string BAMFILENAME,int ExperimentNo){

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
}

PairStruct* PairPool = new PairStruct[BUFFERSIZE];
BamAlignment al,almate;
vector < RESitesClass > dpnIIparallel; // Each thread will have its own RESiteClass

int nofthreads, poolsize = 0;
/*
#pragma omp parallel
{
	nofthreads = omp_get_num_threads(); // Get num of threads
}
*/
	nofthreads = 4;
	cout << "Number of Threads   " << nofthreads << endl;
	for (int t = 0; t < nofthreads; ++t){
		cout << "Thread No " << t << "    ";
		dpnIIparallel.push_back(RESitesClass());
		dpnIIparallel.back().InitialiseVars(); // Initialize all RESiteClass
	}

//	ofstream flagged("flagged_interactions.txt");
//Read the BAM file into PairPool vector
	cout << "Reading BAM file..." << endl;
	PairStruct mate;
while ( reader.GetNextAlignmentCore(al) ){
	GetPairInformation(al,PairPool[poolsize]);
	reader.GetNextAlignmentCore(almate);
	//GetPairInformation(almate,mate);
	++poolsize;
	
	if(poolsize == BUFFERSIZE){
		cout << poolsize << "    Pairs Read" << endl << "Annotating Interactions... ";
#pragma omp parallel for num_threads (nofthreads)
		for(int i=0;i<poolsize;++i){
//			cout << PairPool[i].chr_1 << '\t' << PairPool[i].startcoord << '\t' << PairPool[i].chr_2 << '\t' << PairPool[i].endcoord << endl;
			 int tid = omp_get_thread_num(); // Get thread id
			ProcessPairs(PairPool[i],dpnIIparallel[tid],promoters,negctrls,ExperimentNo); //Annotate the pair
		}
		cout << endl << BAMFILENAME << "   " <<  poolsize << "   pairs finished, printing unusual pairs..." << endl;
/*
		for(int i=0;i<poolsize;++i){
			if (PairPool[i].flagged != 0)
				flagged << PairPool[i].chr_1 << '\t' << PairPool[i].startcoord << '\t' << PairPool[i].chr_2 << '\t' << PairPool[i].endcoord << '\t' << PairPool[i].flagged << endl;
		}
*/
		delete[] PairPool;
		PairPool = new PairStruct[BUFFERSIZE];
		poolsize = 0;
	}
}
if(poolsize > 0){ // This is to read the last batch of pairs
	cout << poolsize << "    Pairs Read" << endl << "Annotating Interactions...";
/*
#pragma omp parallel
	{
		nofthreads = omp_get_num_threads(); // Get num of threads
	}
*/
	nofthreads = 4;
#pragma omp parallel for num_threads (nofthreads)
	for(int i=0;i<poolsize;++i){
	//	cout << "last pool  " << PairPool[i].chr_1 << '\t' << PairPool[i].startcoord << '\t' << PairPool[i].chr_2 << '\t' << PairPool[i].endcoord << endl;
		int tid = omp_get_thread_num(); // Get thread id
		ProcessPairs(PairPool[i],dpnIIparallel[tid],promoters,negctrls,ExperimentNo); //Annotate the pair
	}
	cout << endl << BAMFILENAME << "   reading finished, printing unusual pairs..." << endl;
/*
	for(int i=0;i<poolsize;++i){
		if (PairPool[i].flagged != 0)
			flagged << PairPool[i].chr_1 << '\t' << PairPool[i].startcoord << '\t' << PairPool[i].chr_2 << '\t' << PairPool[i].endcoord << '\t' << PairPool[i].flagged << endl;
	}
*/
	delete[] PairPool;
	PairPool = new PairStruct[BUFFERSIZE];
	poolsize = 0;
}

}



void ProcessBAM::ProcessPairs(PairStruct& thispair, RESitesClass& dpnII, PromoterClass& promoters, NegCtrlClass& negctrls,int ExperimentNo){

	bool pairann = false;
	int *renums1, *renums2;
	int resite_firstinpair, resite_secondinpair;
	bool passed1 = 1, passed2 = 1;
	renums1 = new int [2];	renums2 = new int [2];
	
	bool re1found = dpnII.GettheREPositions(thispair.chr_1,thispair.startcoord,renums1);
	bool re2found = dpnII.GettheREPositions(thispair.chr_2,thispair.endcoord,renums2);
	if(re1found){	
		resite_firstinpair = renums1[1];
		if((abs(thispair.startcoord - renums1[1])) > MaxInsertLen ){
			if(abs((thispair.startcoord - renums1[0])) < MaxInsertLen )
				resite_firstinpair = renums1[0];
			else
				passed1 = 0;
		}
	}
	else
		resite_firstinpair = thispair.startcoord;

	if(re2found){
		resite_secondinpair = renums2[0];
		if(abs((thispair.endcoord - renums2[0])) > MaxInsertLen){ 
			if(abs(renums2[1] - thispair.endcoord) < MaxInsertLen) // if the junction is contained within the pair 
				resite_secondinpair = renums2[1];
			else
				passed2 = 0;
		}
	}
	else
		resite_secondinpair = thispair.endcoord;

	pairann = promoters.AnnotatewithPromoters(thispair.chr_1,resite_firstinpair,thispair.startcoord,thispair.chr_2,resite_secondinpair,thispair.endcoord,dpnII,ExperimentNo);
	if(!pairann)// If both reads are not annotated with promoters
		pairann = negctrls.AnnotateWithNegCtrls(thispair.chr_1,resite_firstinpair,thispair.startcoord,thispair.chr_2,resite_secondinpair,thispair.endcoord,dpnII,ExperimentNo);
	//cout << "..";
	if(passed1 && passed2)
		thispair.flagged = 0;  // Regular Pair
	else // if the REsites are very far away from the read coordinate
		thispair.flagged = 1; // RE sites are far away
	if(!(re1found && re2found))	
		thispair.flagged = 2; // RE sites cannot be found
}