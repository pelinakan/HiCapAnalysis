class ProcessSAM{
	friend class PromoterClass;
	friend class NegCtrlClass;
	friend class ProbeSet;
public:
	void ProcessTheSAMFile(PromoterClass&,NegCtrlClass&,ProbeSet&,string);
private:
	bool GetPairInformation(stringstream&, stringstream&, PairStruct&);
};


bool ProcessSAM::GetPairInformation(stringstream &r1, stringstream &r2, PairStruct& pair){
	
	string chr_1,chr_2,field;
	int start, end;
	bool not_icj=0; //Not a interchromosomal junction
	getline(r1,field,'\t');
	getline(r1,field,'\t'); 
	getline(r1,chr_1,'\t');
	getline(r1,field,'\t');
	start=(atoi(field.c_str())); //start of the pair taken from the first read
	
	getline(r2,field,'\t');
	getline(r2,field,'\t'); 
	getline(r2,chr_2,'\t'); 
	getline(r2,field,'\t');
	end=(atoi(field.c_str())); //End of the pair taken from reverse read
	
	pair.chr_1.append(chr_1);
	pair.startcoord=start; //start of the pair taken from the first read
	pair.chr_2.append(chr_2); 
	pair.endcoord=end; //End of the pair taken from reverse read
	
	if (chr_1 == chr_2)
		return 0;
	else
		return 1;
}


void ProcessSAM::ProcessTheSAMFile(PromoterClass& promoters,NegCtrlClass& negctrls,ProbeSet& mm9prs, string SAMFILENAME){

	ifstream SAMFILE(SAMFILENAME.c_str());
	if(!SAMFILE.is_open())
		cerr << "Input File cannot be opened" << endl;

	PairStruct* PairPool = new PairStruct[BUFFERSIZE];
	int poolsize = 0,nofloops = 0;
	string r;
	bool icj; //not a inter-chromosomal junction
	MappabilityClass mapp;
	do{
		do{
			getline(SAMFILE,r); // Read the forward read
			stringstream read1 ( r );
			getline(SAMFILE,r);
			stringstream read2 ( r ); //Read the reverse read
			icj = GetPairInformation(read1,read2,PairPool[poolsize]);
			++poolsize;
		}while(poolsize < BUFFERSIZE && (!SAMFILE.eof()));
		if(SAMFILE.eof())
			poolsize--;

		cout << BUFFERSIZE  << "    Reads Read" << endl;
		vector < RESitesClass > dpnIIparallel; // Each thread will have its own RESiteClass

		int nofthreads;
#pragma omp parallel
		{
		nofthreads = omp_get_num_threads(); // Get num of threads
		}
		nofthreads = 4;
		cout << "Number of Threads   " << nofthreads << endl;
		for (int t = 0; t < nofthreads; ++t){
			cout << "Thread No " << t << "    ";
			dpnIIparallel.push_back(RESitesClass());
			dpnIIparallel.back().InitialiseVars(); // Initialise all RESiteClass
		}
#pragma omp parallel for num_threads (nofthreads)
		for(int i=0;i<poolsize;++i){
			bool pairann = false;
			bool onprobe_forward = 1;
			bool onprobe_reverse = 1;
			int tid = omp_get_thread_num(); // Get thread id
			//onprobe_forward=mm9prs.AssociateReadwithProbes(PairPool[i].chr_1,PairPool[i].startcoord,PairPool[i].endcoord);
			//onprobe_reverse=mm9prs.AssociateReadwithProbes(PairPool[i].chr_2,PairPool[i].endcoord,PairPool[i].startcoord);
			if (onprobe_forward || onprobe_reverse){
				 pairann = promoters.AnnotatewithPromoters(PairPool[i].chr_1,PairPool[i].startcoord,PairPool[i].chr_2,PairPool[i].endcoord,dpnIIparallel[tid],mapp);
				if(!pairann)// If both reads are not annotated with promoters
					pairann = negctrls.AnnotateWithNegCtrls(PairPool[i].chr_1,PairPool[i].startcoord,PairPool[i].chr_2,PairPool[i].endcoord,dpnIIparallel[tid],mapp);
			}
		}

		for (int t = 0; t < nofthreads; ++t){
			dpnIIparallel[t].posvector.clear();
		}	

		TotalNumberofPairs += poolsize;
		delete[] PairPool;
		PairPool = new PairStruct[BUFFERSIZE];
		poolsize = 0;
		++nofloops;
		cout << BUFFERSIZE*nofloops << "   Pairs Processed" << endl;
//		if(BUFFERSIZE*nofloops == 20000000)	
//			break;
	}while(!SAMFILE.eof());
	cout << SAMFILENAME << "      Read" << endl;
	SAMFILE.close();
}
