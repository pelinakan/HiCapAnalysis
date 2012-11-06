struct temppeak{
	string chr;
	int start;
	int end;
	double signal;
};

class PeakClass {
	friend class PromoterClass;
	friend class NegCtrlClass;
public:
	vector < string > abnames;
	vector < PeakMap > p;
 	vector <int> ChrRowStartIndexes;
	vector <int> ChrRowEndIndexes;
	vector <int> peakcenters;
	vector <double> peaksignals;
	vector < bool > activepeak;
	void InitializeData(void);
	void ReadPeakFile(ifstream&,bool, RESitesClass&,vector<string>&);
	void AnnotatewithPromoters(PromoterClass&,int);
	void AnnotatewithNegCtrls(NegCtrlClass&,int);
	void AssociateProbeswithPeaks(ProbeSet&);

	~PeakClass();
private:
	void GetPeakFeats(stringstream&,temppeak&,bool);
	void PrintEnrichments(ofstream&,vector<double>,vector<double>&,vector<int>);
};
PeakClass::~PeakClass(){
	abnames.clear();
//	ChrNames.clear();
	ChrRowStartIndexes.clear();
	ChrRowEndIndexes.clear();
	peakcenters.clear();
	peaksignals.clear();
	activepeak.clear();
}

void PeakClass::InitializeData(void){



}
void PeakClass::ReadPeakFile(ifstream &peakfile, bool ifBEDFile, RESitesClass& dpnIIsites, vector<string>& cn){
	string pline,chr1,chr2;
	temppeak t;
	long int index=0,peakm;
	int *renums;
	renums = new int [2];

	getline(peakfile,pline); // Read the peak line
	stringstream peak ( pline ); 
	GetPeakFeats(peak,t,ifBEDFile);
	p.push_back(PeakMap()); // Create a new PeakMap for a new chromosome
	p.back().maptochrname = t.chr; // Write chr name
	peakm=t.start+((t.end-t.start)/2);
	dpnIIsites.GettheRENums(t.chr,peakm,renums); // Get the closest REsite
	p.back().peaks.insert(renums[0]); // Insert the closest RE site to the map
	chr1=t.chr;
	bool foundbefore = 0;
	for(int c = 0;c < cn.size();++c ){
		if(cn[c] == chr1 ){
			foundbefore = 1;
		}
	}
	if(!foundbefore)
		cn.push_back(chr1);
	do{
		getline(peakfile,pline);
		peak.clear();
		peak.str(pline);
		GetPeakFeats(peak,t,ifBEDFile);
		chr2 = t.chr;
		while(chr1==chr2 && pline!="" && (!peakfile.eof())){
			peakm=t.start+((t.end-t.start)/2);
			dpnIIsites.GettheRENums(t.chr,peakm,renums); // Get the closest REsite
			p.back().peaks.insert(renums[0]); // Insert the closest RE site to the map
			getline(peakfile,pline);
			peak.clear();
			peak.str(pline);
			GetPeakFeats(peak,t,ifBEDFile);
			chr1 = chr2;
			chr2 = t.chr;
		}
		peakm=t.start+((t.end-t.start)/2);
		if(pline == "" || peakfile.eof() || peakm == 0)
			break;
		p.push_back(PeakMap()); // Create a new PeakMap for a new chromosome
		p.back().maptochrname = t.chr; // Write chr name
		dpnIIsites.GettheRENums(t.chr,peakm,renums); // Get the closest REsite
		p.back().peaks.insert(renums[0]); // Insert the closest RE site to the map
		chr1=t.chr;
		bool foundbefore = 0;
		for(int c = 0;c < cn.size();++c ){
			if(cn[c] == chr1 ){
				foundbefore = 1;
			}
		}
		if(!foundbefore)
			cn.push_back(chr1);
		cout << chr1 << endl;
	}while(!peakfile.eof()&& pline!="");
}

void PeakClass::GetPeakFeats(stringstream &line, temppeak &peak, bool BED){
string field;
if(BED){
	getline(line,field,'\t');
	getline(line,peak.chr,'\t');
	getline(line,field,'\t');
	peak.start=atoi(field.c_str());
	getline(line,field,'\t');
	peak.end=atoi(field.c_str());

	getline(line,field,'\t');
	getline(line,field,'\t');
	getline(line,field,'\t');
	getline(line,field,'\t');

	getline(line,field,'\t');
	peak.signal=atof(field.c_str()); 		
}
else{
	getline(line,peak.chr,'\t');
	getline(line,field,'\t');
	peak.start=atoi(field.c_str());
	getline(line,field,'\t');
	peak.end=atoi(field.c_str());

	getline(line,field,'\t');
	peak.signal=atof(field.c_str()); 		
}

}

void PeakClass::AnnotatewithPromoters(PromoterClass &Prs,int abindex){
int i=0,j=0,k=0;
int diff;
long int startsearch,endsearch;

vector < vector < int > > PeakBins;

for(i=0;i<Prs.NofPromoters;++i){
	startsearch=0;endsearch=-1;
//	for(j=0;j<ChrNames.size();++j){
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


