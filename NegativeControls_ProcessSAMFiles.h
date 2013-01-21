

class NegCtrlClass{
public:
	NegativeControlStruct *negctrls;
int NofNegCtrls;
void InitialiseData(void);
void FillNegativeCtrls(RESitesClass&, MappabilityClass&);
bool AnnotateWithNegCtrls(string,int,string,int,RESitesClass&, MappabilityClass&);
void BinPeaks(int,int,int);

private:
	void AnnotateDistalInteractor(string, string,int*,int);
	void AnnotateFeatFeatInteraction(int, int);


};
void NegCtrlClass::InitialiseData(){
	NofNegCtrls=366;
	negctrls = new NegativeControlStruct [400];
	for(int i = 0; i < NofNegCtrls;++i){
		negctrls[i].closestREsitenums = new int [2];
	}
}

void NegCtrlClass::FillNegativeCtrls(RESitesClass& dpnIIsites, MappabilityClass& mapp){

#ifdef UNIX
string filename1,filename2;
string dirname="/bubo/proj/b2011029/bin/3CAnalysis/";
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

bool NegCtrlClass::AnnotateWithNegCtrls(string p_chr_1, int pstcoord, string p_chr_2, int pendcoord,RESitesClass& dpnIIsites, MappabilityClass& mappability){

bool pann = false;

p_chr_1.erase(p_chr_1.find_last_not_of(" \n\r\t")+1); //trim the string
p_chr_2.erase(p_chr_2.find_last_not_of(" \n\r\t")+1);

int ncidx1 = -1, ncidx2 = -1; // PromoterID
for(int i = 0; i < NofNegCtrls; ++i){ //Iterate over all refseq genes on that chromosome
	if (negctrls[i].chr.compare(p_chr_1) == 0){
		if((negctrls[i].start <= pstcoord && negctrls[i].end >= pstcoord)){ // If the readstart is contained within an negative control
			ncidx1 = i;
			pann = 1; // Read is annotated with a promoter
			if((p_chr_1.compare(p_chr_2) == 0) && (negctrls[i].start <= pendcoord && negctrls[i].end >= pendcoord)) // if the pair of the read is also contained within the core promoter
				return pann;
			// Check if the pair is close to another negative control
			for(int m = 0; m < NofNegCtrls; ++m){ //Iterate over all refseq genes on that chromosome
				if (negctrls[m].chr.compare(p_chr_2) == 0){
					if((negctrls[m].start <= pendcoord && negctrls[m].end >= pendcoord)){ // If the pair is contained within another negative control
						ncidx2 = m; // It is prom-prom interaction
						AnnotateFeatFeatInteraction(ncidx1, ncidx2);
						return pann;
					}
				} 
			} // Checked if it is negctrl-negctrl interaction
			int *renums;
			renums = new int [2];
			dpnIIsites.GettheREPositions(p_chr_2,pendcoord,renums); // Get the fragment number of the distal interactor
			AnnotateDistalInteractor(p_chr_1,p_chr_2,renums,ncidx1);
			return pann;
		}
	}
}// First read in the pair processed
// Annotate the second read in the pair (this could only be distal-prom interaction)
for(int m = 0; m < NofNegCtrls; ++m){ //Iterate over all refseq genes on that chromosome
	if (negctrls[m].chr.compare(p_chr_2) == 0){
		if((negctrls[m].start <= pendcoord && negctrls[m].end >= pendcoord)){ // If the readstart is contained within an negative control
			ncidx2 = m;
			pann = 1; // Read is annotated with a promoter
			if((p_chr_1.compare(p_chr_2) == 0) && (negctrls[m].start <= pstcoord && negctrls[m].end >= pstcoord)) // if the pair of the read is also contained within the core promoter
				return pann;
			int *renums;
			renums = new int [2];
			dpnIIsites.GettheREPositions(p_chr_1, pstcoord,renums); // Get the fragment number of the distal interactor (this is it is the first read)
			AnnotateDistalInteractor(p_chr_2,p_chr_1,renums,ncidx2);
			return pann;
		}
	}			
}
return pann;
}

void NegCtrlClass::AnnotateDistalInteractor(string p_chr_1, string p_chr_2, int *renums, int ncidx){

	if(p_chr_1.compare(p_chr_2) == 0){ // Intra chromosomal interaction
		if (renums[1] < negctrls[ncidx].closestREsitenums[0]){ // if upstream 
			if(negctrls[ncidx].Signals.signal_ups.find(renums[0]) == negctrls[ncidx].Signals.signal_ups.end())
				negctrls[ncidx].Signals.signal_ups[renums[0]] = 1; // add a new entry
			else // if inserted before
				negctrls[ncidx].Signals.signal_ups[renums[0]] = negctrls[ncidx].Signals.signal_ups[renums[0]] + 1;
		}
		else{
			if(negctrls[ncidx].Signals.signal_down.find(renums[0]) == negctrls[ncidx].Signals.signal_down.end())
				negctrls[ncidx].Signals.signal_down[renums[0]] = 1; // add a new entry
			else // if inserted before
				negctrls[ncidx].Signals.signal_down[renums[0]] = negctrls[ncidx].Signals.signal_down[renums[0]] + 1;
		}		
	}
	else{ // inter chromosomal interaction
		bool chrfound = 0;	
		vector < SignalStruct_CTX >::iterator it;
		for(it = negctrls[ncidx].Signals_CTX.begin(); it < negctrls[ncidx].Signals_CTX.end(); ++it){
			if (p_chr_2.compare(it->maptochrname) == 0){
				if(it->signal.find(renums[0]) == it->signal.end())
					it->signal[renums[0]] = 1; // add a new entry
				else // if inserted before
					it->signal[renums[0]] = it->signal[renums[0]] + 1;
				chrfound = 1;
				break;
			}
		}
		if(!chrfound){
			negctrls[ncidx].Signals_CTX.push_back(SignalStruct_CTX());
			negctrls[ncidx].Signals_CTX.back().maptochrname = p_chr_2;
			negctrls[ncidx].Signals_CTX.back().signal[renums[0]] = 1;
		}
	}
}

void NegCtrlClass::AnnotateFeatFeatInteraction(int ncidx1, int ncidx2){
	bool foundbefore = 0;	
	for (int x = 0; x < negctrls[ncidx1].NegcNegcSignals.size(); ++x){ // Check if the interaction with that promoter is seen before
		if (negctrls[ncidx1].NegcNegcSignals[x].feat_index == ncidx2){
			negctrls[ncidx1].NegcNegcSignals[x].normsignal += 1.0;
			foundbefore = 1;
			break;
		}
	}
	if(!foundbefore){ // Create a new entry for that promoter
		negctrls[ncidx1].NegcNegcSignals.push_back(FeattoFeatSignalStruct());
		negctrls[ncidx1].NegcNegcSignals.back().feat_index = ncidx2;
		negctrls[ncidx1].NegcNegcSignals.back().normsignal += 1.0;
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