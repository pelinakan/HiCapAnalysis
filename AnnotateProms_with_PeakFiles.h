class PeakClass {
	friend class PromoterClass;
	friend class NegCtrlClass;
public:
	vector < string > abnames;
	boost::unordered::unordered_map< string, int > chr_indexes;
	vector < PeakMap > p;
	void ReadPeakFile(ifstream&, RESitesClass&,vector<string>&);
	void AnnotatewithPromoters(PromoterClass&,int);
//	void AnnotatewithNegCtrls(NegCtrlClass&,int);
//	void AssociateProbeswithPeaks(ProbeSet&);

};
void PeakClass::ReadPeakFile(ifstream &peakfile, RESitesClass& dpnIIsites, vector<string>& cn){
	string chr1,chr2,temp,pchr;
	int st,end;
	int *renums;
	renums = new int [2];

	peakfile >> pchr >> st >> end;
	getline(peakfile,temp);
	p.push_back(PeakMap()); // Create a new PeakMap for a new chromosome
	chr_indexes[pchr] = (p.size()-1); // Map chr name to peak array index
	dpnIIsites.GettheREPositions(pchr,(st + ((end-st)/2)),renums); // Get the closest REsite
	p.back().peaks.insert(renums[0]); // Insert the closest RE site to the map
	chr1 = pchr;
	bool foundbefore = 0;
	for(int c = 0;c < cn.size();++c ){
		if(cn[c] == chr1 ){
			foundbefore = 1;
		}
	}
	if(!foundbefore)
		cn.push_back(chr1);
	do{
		peakfile >> pchr >> st >> end;
		getline(peakfile,temp);
		chr2 = pchr;
		while(chr1==chr2 && (!peakfile.eof())){
			dpnIIsites.GettheREPositions(pchr,(st + ((end-st)/2)),renums); // Get the closest REsite
			p.back().peaks.insert(renums[0]); // Insert the closest RE site to the map
			peakfile >> pchr >> st >> end;
			getline(peakfile,temp);
			chr1 = chr2;
			chr2 = pchr;
		}
		if(peakfile.eof() || (st + ((end-st)/2) == 0))
			break;
		p.push_back(PeakMap()); // Create a new PeakMap for a new chromosome
		chr_indexes[pchr] = (p.size()-1); // Map chr name to peak array index
		dpnIIsites.GettheREPositions(pchr,(st + ((end-st)/2)),renums); // Get the closest REsite
		p.back().peaks.insert(renums[0]); // Insert the closest RE site to the map
		chr1=pchr;
		bool foundbefore = 0;
		for(int c = 0;c < cn.size();++c ){
			if(cn[c] == chr1 ){
				foundbefore = 1;
			}
		}
		if(!foundbefore)
			cn.push_back(chr1);
		cout << chr1 << endl;
	}while(!peakfile.eof());
}

/*
void PeakClass::AnnotatewithPromoters(PromoterClass &Prs,int abindex){
	
	
	for(int i = 0;i < Prs.NofPromoters; ++i){
	for(j=0;j<ChrNames.size();++j){
//		if(Prs.refseq[i].chr==ChrNames[j]){
			startsearch=ChrRowStartIndexes[j];
			endsearch=ChrRowEndIndexes[j];
//			break;
//		}
//	}
	for(k=0;k<Prs.refseq[i].isoformpromotercoords.size();++k)
		PeakBins.push_back(vector <int>());
	for(j=startsearch;j<=endsearch;++j){
		for(k=0;k<Prs.refseq[i].isoformpromotercoords.size();++k){
			diff=peakcenters[j]-Prs.refseq[i].isoformpromotercoords[k];
			if(abs(diff)<MaxJunctionDistance)
//				Prs.BinPeaks(i,peakcenters[j],k,abindex);
				diff = 0;
		}
	}
	if(i%5000==0)
		cout << i << "    Promoters Associated with Peaks" << endl;
}
}

void PeakClass::AnnotatewithNegCtrls(NegCtrlClass &ng,int abindex){
int i=0,j=0,k=0;
long int startsearch,endsearch;

for(i=0;i<ng.NofNegCtrls;++i){
	startsearch=0;endsearch=-1;
//	for(j=1;j<NumberofBins;++j)
//		ng.negctrls[i].AllPeaks_PeakBins[abindex].PeakBins[0][j]=0;
//	for(j=0;j<ChrNames.size();++j){
//		if(ng.negctrls[i].chr==ChrNames[j]){
			startsearch=ChrRowStartIndexes[j];
			endsearch=ChrRowEndIndexes[j];
//			break;
//		}
//	}
	for(j=startsearch;j<=endsearch;++j){
		if((abs(peakcenters[j]-ng.negctrls[i].midpoint))<MaxJunctionDistance)
	//		ng.BinPeaks(i,peakcenters[j],abindex);
		j = 0;
	}
	if(i%100==0)
		cout << i << "    Negative Control Regions Associated with Peaks" << endl;
}
}
void PeakClass::AssociateProbeswithPeaks(ProbeSet& mm9prs){

int i=0,j=0,k=0,x=0;
int diff;

for(i=0;i<ChrRowStartIndexes.size();++i){
	for(j=ChrRowStartIndexes[i];j<=ChrRowEndIndexes[i];++j){
		k=0;
//		while(ChrNames[i]!=mm9prs.ChrNames[k])
			++k;
		for(x=mm9prs.ChrRowStartIndexes[k];x<=mm9prs.ChrRowEndIndexes[k];++x){
			diff=mm9prs.MM9Probes[x].center-peakcenters[j];
			if(abs(diff)<BinSize){
				if(mm9prs.MM9Probes[x].annotated==0)
					mm9prs.MM9Probes[x].annotated=2; // Associate to a peak
			}
			if(diff>BinSize)
				break;

	}
		if(j%1000==0)
			cout << j << "    Probes Associated with a Peak" << endl;
	}
}

}


*/