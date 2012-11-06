class ProcessSAM{
	friend class PromoterClass;
	friend class NegCtrlClass;
	friend class ProbeSet;
public:
	void ProcessTheSAMFile(PromoterClass&,NegCtrlClass&,ProbeSet&,RESitesClass&,MappabilityClass&);
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


void ProcessSAM::ProcessTheSAMFile(PromoterClass& promoters,NegCtrlClass& negctrls,ProbeSet& mm9prs,RESitesClass& dpnII, MappabilityClass& mapp){
	PairStruct *PairPool;
	PairPool = new PairStruct[BUFFERSIZE];

	int poolsize=0,nofloops=0;
	string r;
	bool icj; //not a interchromosomal junction
	do{
		do{
			getline(SAMFILE,r); // Read the forward read
			stringstream read1 ( r );
			getline(SAMFILE,r);
			stringstream read2 ( r ); //Read the reverse read
			icj=GetPairInformation(read1,read2,PairPool[poolsize]);
			++poolsize;
		}while(poolsize<BUFFERSIZE && (!SAMFILE.eof()));
		if(SAMFILE.eof())
			poolsize--;
	//#pragma omp parallel for default(shared)
		for(int i=0;i<poolsize;++i){
			bool pair1ann=false,pair2ann=false;
			bool annotateMate=1;
			bool onprobe_forward = 1;
			bool onprobe_reverse = 1;
			//onprobe_forward=mm9prs.AssociateReadwithProbes(PairPool[i].chr_1,PairPool[i].startcoord,PairPool[i].endcoord);
			//onprobe_reverse=mm9prs.AssociateReadwithProbes(PairPool[i].chr_2,PairPool[i].endcoord,PairPool[i].startcoord);
			if (onprobe_forward || onprobe_reverse){
				pair1ann=promoters.AnnotatewithPromoters(PairPool[i].chr_1,PairPool[i].startcoord,PairPool[i].chr_2,PairPool[i].endcoord,annotateMate,dpnII,mapp,0,0);
				if(annotateMate)
					pair2ann=promoters.AnnotatewithPromoters(PairPool[i].chr_2,PairPool[i].endcoord,PairPool[i].chr_2,PairPool[i].startcoord,annotateMate,dpnII,mapp, 1,pair1ann);		
				/*
				if(!pair1ann && !pair2ann){
					pair1ann=negctrls.AnnotateWithNegCtrls(PairPool[i].chr_1,PairPool[i].startcoord,PairPool[i].endcoord,annotateMate,dpnII, mapp);
					if(annotateMate)
						pair2ann=negctrls.AnnotateWithNegCtrls(PairPool[i].chr_2,PairPool[i].endcoord,PairPool[i].startcoord,annotateMate,dpnII,mapp);
				}
				*/
			}
		}
		delete []PairPool;
		PairPool = new PairStruct[BUFFERSIZE];
		poolsize=0;
		++nofloops;
		cout << BUFFERSIZE*nofloops << "   Pairs Processed" << endl;
	}while(!SAMFILE.eof());



}
