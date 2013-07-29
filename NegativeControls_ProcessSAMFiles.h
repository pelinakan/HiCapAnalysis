
class NegCtrlClass{
public:
	NegativeControlStruct *negctrls;
	vector <string> ChrNames_enhancers;
	vector < vector <int> > enhancer_indexes; //Based on enhancer chromosome
int NofNegCtrls;
void InitialiseData(int);
void FillNegativeCtrls(RESitesClass&, MappabilityClass&);
void FillEnhancers_asNegCtrls(RESitesClass&, MappabilityClass&,string);
bool AnnotateWithNegCtrls(string,int,int,string,int,int,RESitesClass&,int);
bool AnnotateWithEnhancers(string,int,int,string,int,int,RESitesClass&,int);
void BinPeaks(int,int,int);

private:
	void AnnotateDistalInteractor(string, string,int,int,int,int,int);
	void AnnotateFeatFeatInteraction(int, int,int);
	void PopulateInteractions(boost::unordered::unordered_map<int, int* >&, int, int, int);

};
void NegCtrlClass::InitialiseData(int size){
	negctrls = new NegativeControlStruct [size];
	for(int i = 0; i < size;++i)
		negctrls[i].closestREsitenums = new int [2];
}

void NegCtrlClass::FillNegativeCtrls(RESitesClass& dpnIIsites, MappabilityClass& mapp){

#ifdef UNIX
string filename1,filename2;
filename1.append(dirname);
filename1.append("100exons_min100kfromTSS_GATC_150perend.bed");
ifstream infile1(filename1.c_str());
filename2.append(dirname);
filename2.append("100intergenic_min100kfromTSS_GATC_150perend.bed");
ifstream infile2(filename2.c_str());
#endif

#ifdef WINDOWS
	ifstream infile1("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\100exons_min100kfromTSS_GATC_150perend.bed");
	ifstream infile2("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\100intergenic_min100kfromTSS_GATC_150perend.bed");
#endif
int i=0;
string chr;
do{
	infile1 >>	negctrls[i].chr >> negctrls[i].start >> negctrls[i].end;
	negctrls[i].midpoint=(abs((negctrls[i].end-negctrls[i].start))/2)+negctrls[i].start;
	negctrls[i].start = negctrls[i].midpoint - AssociateInteractions;
	negctrls[i].end = negctrls[i].midpoint + AssociateInteractions;
	negctrls[i].type="G";
	++i;
}while(!infile1.eof());

infile1.close();
i--;

do{
	infile2 >> negctrls[i].chr >> negctrls[i].start >> negctrls[i].end;
	negctrls[i].midpoint=(abs((negctrls[i].end-negctrls[i].start))/2)+negctrls[i].start;
	negctrls[i].start = negctrls[i].midpoint - AssociateInteractions;
	negctrls[i].end = negctrls[i].midpoint + AssociateInteractions;
	negctrls[i].type="I";
	++i;

}while(!infile2.eof());
infile2.close();
NofNegCtrls=i-1;

	for (i = 0; i < NofNegCtrls; ++i){
		negctrls[i].nofRESites = dpnIIsites.GetRESitesCount(negctrls[i].chr,negctrls[i].start,negctrls[i].end);
		dpnIIsites.GettheREPositions(negctrls[i].chr,negctrls[i].midpoint,negctrls[i].closestREsitenums);
		negctrls[i].featmappability = mapp.GetMappability(negctrls[i].chr, negctrls[i].start, (negctrls[i].end - negctrls[i].start));
	}
}

void NegCtrlClass::FillEnhancers_asNegCtrls(RESitesClass& dpnIIsites, MappabilityClass& mapp,string whichchr){

#ifdef UNIX
	string filename1,filename2;
	filename1.append(dirname);
	filename1.append("HiCap.merged.enhancers.txt");
	ifstream infile1(filename1.c_str());
#endif

#ifdef WINDOWS
	ifstream infile1("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\100exons_min100kfromTSS_GATC_150perend.bed");
#endif
	int i=0,midpoint, rest,reend;
	string chr;
if(whichchr == "CTX"){
	do{
		infile1 >> chr >> midpoint >> rest >> reend;
		negctrls[i].chr = chr;
		negctrls[i].midpoint = midpoint;
		negctrls[i].closestREsitenums[0] = rest;
		negctrls[i].closestREsitenums[1] = reend;
		negctrls[i].start = midpoint - AssociateInteractions;
		negctrls[i].end = midpoint + AssociateInteractions;
		negctrls[i].type="E";
		++i;
	}while(!infile1.eof());
	infile1.close();
	NofNegCtrls=i-1;
	cout << NofNegCtrls << "    Enhancer Read" << endl;
}

else{
	do{
		infile1 >> chr >> midpoint >> rest >> reend;
		if(chr == whichchr){
			negctrls[i].chr = chr;
			negctrls[i].midpoint = midpoint;
			negctrls[i].closestREsitenums[0] = rest;
			negctrls[i].closestREsitenums[1] = reend;
			negctrls[i].start = midpoint - AssociateInteractions;
			negctrls[i].end = midpoint + AssociateInteractions;
			negctrls[i].type="E";
			++i;
		}
	}while(!infile1.eof());
	infile1.close();
	NofNegCtrls=i-1;
	cout << NofNegCtrls << "    Enhancer Read" << endl;
}
	for (i = 0; i < NofNegCtrls; ++i){
		negctrls[i].nofRESites = dpnIIsites.GetRESitesCount(negctrls[i].chr,negctrls[i].start,negctrls[i].end);
//		dpnIIsites.GettheREPositions(negctrls[i].chr,negctrls[i].midpoint,negctrls[i].closestREsitenums);
//		negctrls[i].featmappability = mapp.GetMappability(negctrls[i].chr, negctrls[i].start, (negctrls[i].end - negctrls[i].start));
		negctrls[i].featmappability = 0;
	}

//Index Enhancers
	int found = 0;
	for(int eindex = 0; eindex < NofNegCtrls; ++eindex){
		found = -1;
		for(i = 0; i < ChrNames_enhancers.size();++i){
			if(negctrls[eindex].chr.compare(ChrNames_enhancers[i]) == 0) 
				found = i;
		}
		if(found == -1){
			ChrNames_enhancers.push_back(negctrls[eindex].chr);
			enhancer_indexes.push_back(vector<int>());
			enhancer_indexes.back().push_back(eindex);
		}
		else
			enhancer_indexes[found].push_back(eindex);

	}
	cout << "Enhancer Chromosome Indexes Generated" << endl;


}
bool NegCtrlClass::AnnotateWithNegCtrls(string p_chr_1, int re_firstinpair, int pos_firstinpair, string p_chr_2, int re_secondinpair, int pos_secondinpair, RESitesClass& dpnIIsites,int ExperimentNo){

bool pann = false;

p_chr_1.erase(p_chr_1.find_last_not_of(" \n\r\t")+1); //trim the string
p_chr_2.erase(p_chr_2.find_last_not_of(" \n\r\t")+1);

int ncidx1 = -1, ncidx2 = -1; // PromoterID
for(int i = 0; i < NofNegCtrls; ++i){ //Iterate over all refseq genes on that chromosome
	if (negctrls[i].chr.compare(p_chr_1) == 0){
		if((negctrls[i].start <= pos_firstinpair && negctrls[i].end >= pos_firstinpair)){ // If the readstart is contained within an negative control
			ncidx1 = i;
			pann = 1; // Read is annotated with a promoter
			if((p_chr_1.compare(p_chr_2) == 0) && (negctrls[i].start <= pos_secondinpair && negctrls[i].end >= pos_secondinpair)) // if the pair of the read is also contained within the core promoter
				return pann;
			// Check if the pair is close to another negative control
			for(int m = 0; m < NofNegCtrls; ++m){ //Iterate over all refseq genes on that chromosome
				if (negctrls[m].chr.compare(p_chr_2) == 0){
					if((negctrls[m].start <= pos_secondinpair && negctrls[m].end >= pos_secondinpair)){ // If the pair is contained within another negative control
						ncidx2 = m; // It is prom-prom interaction
						AnnotateFeatFeatInteraction(ncidx1, ncidx2,ExperimentNo);
						return pann;
					}
				} 
			} // Checked if it is negctrl-negctrl interaction
			AnnotateDistalInteractor(p_chr_1,p_chr_2,re_secondinpair,pos_firstinpair, pos_secondinpair, ncidx1,ExperimentNo);			
			return pann;
		}
	}
}// First read in the pair processed
// Annotate the second read in the pair (this could only be distal-prom interaction)
for(int m = 0; m < NofNegCtrls; ++m){ //Iterate over all refseq genes on that chromosome
	if (negctrls[m].chr.compare(p_chr_2) == 0){
		if((negctrls[m].start <= pos_secondinpair && negctrls[m].end >= pos_secondinpair)){ // If the readstart is contained within an negative control
			ncidx2 = m;
			pann = 1; // Read is annotated with a promoter
			AnnotateDistalInteractor(p_chr_2,p_chr_1,re_firstinpair,pos_secondinpair,pos_firstinpair,ncidx2,ExperimentNo); 
			return pann;
		}
	}			
}
return pann;
}
bool NegCtrlClass::AnnotateWithEnhancers(string p_chr_1, int re_firstinpair, int pos_firstinpair, string p_chr_2, int re_secondinpair, int pos_secondinpair, RESitesClass& dpnIIsites,int ExperimentNo){
	int enh1_chrindex, enh2_chrindex;
	bool pann = 0;
	bool chrfound1 = 0, chrfound2 = 0;

	p_chr_1.erase(p_chr_1.find_last_not_of(" \n\r\t")+1); //trim the string
	for(int k = 0;k < ChrNames_enhancers.size(); ++k){
		if(ChrNames_enhancers[k].compare(p_chr_1.c_str()) == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
			enh1_chrindex = k;
			chrfound1 = 1;
			break;
		}
	}
	p_chr_2.erase(p_chr_2.find_last_not_of(" \n\r\t")+1);
	for(int k = 0;k < ChrNames_enhancers.size(); ++k){
		if(ChrNames_enhancers[k].compare(p_chr_2.c_str()) == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
			enh2_chrindex = k;
			chrfound2 = 1;
			break;
		}
	}

	int eidx1 = -1, eidx2 = -1; // PromoterID
	if(chrfound1){
		for(int i = 0; i < enhancer_indexes[enh1_chrindex].size(); ++i){ //Iterate over all refseq genes on that chromosome
			if((negctrls[enhancer_indexes[enh1_chrindex][i]].start <= pos_firstinpair && negctrls[enhancer_indexes[enh1_chrindex][i]].end >= pos_firstinpair)){ // If the readstart is contained within the core promoter
				eidx1 =	enhancer_indexes[enh1_chrindex][i];
				pann = 1; // Read is annotated with a promoter
				if(p_chr_1.compare(p_chr_2) == 0 && (negctrls[eidx1].start <= pos_secondinpair && negctrls[eidx1].end >= pos_secondinpair))
					return pann; // if the pair of the read is also contained within the core promoter
				if(chrfound2){	// Check if the pair is close to another promoter
					for(int m = 0; m < enhancer_indexes[enh2_chrindex].size(); ++m){ //Iterate over all refseq genes on that chromosome
						if((negctrls[enhancer_indexes[enh2_chrindex][m]].start <= pos_secondinpair && negctrls[enhancer_indexes[enh2_chrindex][m]].end >= pos_secondinpair)){ // If the readstart is contained within the core promoter
							eidx2 =	enhancer_indexes[enh2_chrindex][m]; // It is prom-prom interaction
							AnnotateFeatFeatInteraction(eidx1, eidx2,ExperimentNo);
							AnnotateFeatFeatInteraction(eidx2, eidx1,ExperimentNo);
							return pann;
						}
					}
				} // Checked if it is prom-prom interaction
				AnnotateDistalInteractor(p_chr_1,p_chr_2,re_secondinpair,pos_firstinpair, pos_secondinpair, eidx1,ExperimentNo);
				return pann;
			}
		}
	}// First read in the pair processed
	// Annotate the second read in the pair (this could only be distal-prom interaction)
	if(chrfound2){
		for(int m = 0; m < enhancer_indexes[enh2_chrindex].size(); ++m){ //Iterate over all refseq genes on that chromosome
			if((negctrls[enhancer_indexes[enh2_chrindex][m]].start <= pos_secondinpair && negctrls[enhancer_indexes[enh2_chrindex][m]].end >= pos_secondinpair)){ // If the readstart is contained within the core promoter
				eidx2 =	enhancer_indexes[enh2_chrindex][m];
				//	if (proms[pidx2].genes[0] == "0610005C13Rik")
				//		cout << proms[pidx2].genes[0] << '\t' << p_chr_2 << '\t' << p_chr_1 << '\t' << pos_secondinpair << '\t' << pos_firstinpair << '\t' << resite_secondinpair << '\t' << resite_firstinpair << endl;
				pann = 1;
				AnnotateDistalInteractor(p_chr_2,p_chr_1,re_firstinpair,pos_secondinpair,pos_firstinpair,eidx2,ExperimentNo); 
				return pann;
			}			
		}
	}
	return pann;

}
void NegCtrlClass::PopulateInteractions(boost::unordered::unordered_map<int, int* >& signals, int interactor_resite, int interactor_pos,int ExperimentNo){

	if(signals.find(interactor_resite) == signals.end()){
		signals[interactor_resite] = new int[NOFEXPERIMENTS + 1]; // add a new entry
		for(int z = 0; z < (NOFEXPERIMENTS); ++z)
			signals[interactor_resite][z + 1] = 0;
		signals[interactor_resite][ExperimentNo + 1] = 1;
	}
	else{ // if inserted before	
		(signals[interactor_resite][ExperimentNo + 1]) = signals[interactor_resite][ExperimentNo + 1] + 1.0;
	}
	signals[interactor_resite][0] = interactor_pos; // Actual read pos
}

void NegCtrlClass::AnnotateDistalInteractor(string p_chr_1, string p_chr_2, int resite_2, int pos_1, int pos_2, int ncidx,int ExperimentNo){

	if(p_chr_1.compare(p_chr_2) == 0){ // Intra chromosomal interaction
		if (pos_2 < pos_1) // if upstream 
			PopulateInteractions(negctrls[ncidx].Signals.signal_ups,resite_2, pos_2,ExperimentNo);

		else{
			PopulateInteractions(negctrls[ncidx].Signals.signal_down,resite_2, pos_2,ExperimentNo);
		}		
	}
	else{  // inter chromosomal interaction
		bool chrfound = 0;
		vector < SignalStruct_CTX >::iterator itx;
		for(itx = negctrls[ncidx].Signals_CTX.begin() ; itx < negctrls[ncidx].Signals_CTX.end();++itx){
			if (p_chr_2.compare(itx->maptochrname) == 0){
				if(itx->signal.find(resite_2) == itx->signal.end()){
					itx->signal[resite_2] = new int[NOFEXPERIMENTS + 1]; // add a new entry
					for(int z = 0; z < NOFEXPERIMENTS; ++z)
						itx->signal[resite_2][z + 1] = 0;
					itx->signal[resite_2][ExperimentNo + 1] = 1;
					itx->signal[resite_2][0] = pos_2; // Actual read pos
				}
				else{ // if inserted before
					itx->signal[resite_2][ExperimentNo + 1] += 1;
					itx->signal[resite_2][0] = pos_2; // Actual read pos
				}
				chrfound = 1;
				break;
			}
		}
		if(!chrfound){
			negctrls[ncidx].Signals_CTX.push_back(SignalStruct_CTX());
			negctrls[ncidx].Signals_CTX.back().maptochrname = p_chr_2;
			negctrls[ncidx].Signals_CTX.back().signal[resite_2] = new int[NOFEXPERIMENTS + 1];
			for(int z = 0; z < NOFEXPERIMENTS; ++z)
				negctrls[ncidx].Signals_CTX.back().signal[resite_2][z +1] = 0;
			negctrls[ncidx].Signals_CTX.back().signal[resite_2][ExperimentNo + 1] = 1;
			negctrls[ncidx].Signals_CTX.back().signal[resite_2][0] = pos_2;
		}
	}	 
}

void NegCtrlClass::AnnotateFeatFeatInteraction(int ncidx1, int ncidx2,int ExperimentNo){
	bool foundbefore = 0;	
	for (int x = 0; x < negctrls[ncidx1].NegcNegcSignals.size(); ++x){ // Check if the interaction with that promoter is seen before
		if (negctrls[ncidx1].NegcNegcSignals[x].feat_index == ncidx2){
			negctrls[ncidx1].NegcNegcSignals[x].normsignal[ExperimentNo] += 1.0;
			foundbefore = 1;
			break;
		}
	}
	if(!foundbefore){ // Create a new entry for that promoter
		negctrls[ncidx1].NegcNegcSignals.push_back(FeattoFeatSignalStruct());
		negctrls[ncidx1].NegcNegcSignals.back().feat_index = ncidx2;
		negctrls[ncidx1].NegcNegcSignals.back().normsignal[ExperimentNo] += 1.0;
	}
}
/*
void NegCtrlClass::AssociateProbeswithNegativeControls(ProbeSet& MM9prs){


int i=0,j=0,k=0,x=0;
int diff;
long int startsearch,endsearch;


for(i=0;i<NofNegCtrls;++i){
	startsearch=0;endsearch=-1;
	for(j=0;j<MM9prs.ChrNames.size();++j){
		if(negctrls[i].chr==MM9prs.ChrNames[j]){
			startsearch=MM9prs.ChrRowStartIndexes[j];
			endsearch=MM9prs.ChrRowEndIndexes[j];
			break;
		}
	}
	x=0;
	for(j=startsearch;j<=endsearch;++j){
		diff=MM9prs.MM9Probes[j].center-negctrls[i].midpoint;
		if(abs(diff)<AssociateInteractions){
				negctrls[i].ProbeIDs.push_back(j);

		}
		if(diff>BinSize)
			break;
	}
	if(i%100==0)
		cout << i << "    Probes Associated with Negative Controls" << endl;
}
}


*/