class ProcessSAM{
	friend class PromoterClass;
	friend class NegCtrlClass;
	friend class ProbeSet;
public:
	void ProcessTheSAMFile(PromoterClass&,NegCtrlClass&,ProbeSet&,string,int);
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
	end=(atoi(field.c_str())); //start of the pair taken from reverse read
	
	pair.chr_1.append(chr_1);
	pair.startcoord=start; //start of the pair taken from the first read
	pair.chr_2.append(chr_2); 
	pair.endcoord=end; //End of the pair taken from reverse read
	
	if (chr_1 == chr_2)
		return 0;
	else
		return 1;
}


void ProcessSAM::ProcessTheSAMFile(PromoterClass& promoters,NegCtrlClass& negctrls,ProbeSet& mm9prs, string SAMFILENAME,int ExperimentNo){

//	ofstream flagged("flagged_ints2.txt");
	ifstream SAMFILE(SAMFILENAME.c_str());
	if(!SAMFILE.is_open())
		cerr << "Input File cannot be opened" << endl;

	PairStruct* PairPool = new PairStruct[BUFFERSIZE];
	cout << "PairStruct Generated" << endl;
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

		int nofthreads = 1;

#pragma omp parallel
		{
		nofthreads = omp_get_num_threads(); // Get num of threads
		}

		cout << "Number of Threads   " << nofthreads << endl;
		for (int t = 0; t < nofthreads; ++t){
			cout << "Thread No " << t << "    ";
			dpnIIparallel.push_back(RESitesClass());
			dpnIIparallel.back().InitialiseVars(); // Initialise all RESiteClass
		}
#pragma omp parallel for num_threads (nofthreads)
		for(int i=0;i<poolsize;++i){
			bool pairann = false;
			int tid = omp_get_thread_num(); // Get thread id
			int *renums1, *renums2;
			int resite_firstinpair, resite_secondinpair;
			bool passed1 = 1, passed2 = 1;
			renums1 = new int [2];
			renums2 = new int [2];
			bool re1found = dpnIIparallel[tid].GettheREPositions(PairPool[i].chr_1,PairPool[i].startcoord,renums1);
			bool re2found = dpnIIparallel[tid].GettheREPositions(PairPool[i].chr_2,PairPool[i].endcoord,renums2);
			if(re1found && re2found){	
				resite_firstinpair = renums1[1];
				if((abs(PairPool[i].startcoord - renums1[1])) > MaxInsertLen ){
					if(abs((PairPool[i].startcoord - renums1[0])) < MaxInsertLen )
						resite_firstinpair = renums1[0];
					else
						passed1 = 0;
				}
				resite_secondinpair = renums2[0];
				if(abs((PairPool[i].endcoord - renums2[0])) > MaxInsertLen){ 
					if(abs(renums2[1] - PairPool[i].endcoord) < MaxInsertLen) // if the junction is contained within the pair 
						resite_secondinpair = renums2[1];
					else
						passed2 = 0;
				}
				if(passed1 && passed2){
					pairann = promoters.AnnotatewithPromoters(PairPool[i].chr_1,resite_firstinpair,PairPool[i].chr_2,resite_secondinpair,dpnIIparallel[tid],mapp,ExperimentNo);
					if(!pairann)// If both reads are not annotated with promoters
						pairann = negctrls.AnnotateWithNegCtrls(PairPool[i].chr_1,resite_firstinpair,PairPool[i].chr_2,resite_secondinpair,dpnIIparallel[tid],mapp,ExperimentNo);
				}
		//		if(i % 5000 == 0)
		//			cout << i << "pairs processed" << endl;
	/*
				else{ // if the REsites are very far away from the read coordinate
					flagged << PairPool[i].chr_1 << '\t' << PairPool[i].startcoord << '\t' <<  renums1[0] << '\t' << renums1[1] << '\t';
					flagged	<< PairPool[i].chr_2 << '\t' << PairPool[i].endcoord << '\t' <<  renums2[0] << '\t' << renums2[1] << endl;
				}
	*/
			}
		}
		TotalNumberofPairs += poolsize;
		delete[] PairPool;
		PairPool = new PairStruct[BUFFERSIZE];
		poolsize = 0;
		++nofloops;
		cout << BUFFERSIZE*nofloops << "   Pairs Processed" << endl;
	//	break;
	}while(!SAMFILE.eof());
	cout << SAMFILENAME << "      Read" << endl;
	SAMFILE.close();
}
