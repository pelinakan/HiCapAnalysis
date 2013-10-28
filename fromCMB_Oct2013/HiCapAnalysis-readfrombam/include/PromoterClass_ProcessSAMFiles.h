struct temppars{
		string chr;
		int start;
		int end;
		double expression[3];
		string strand;
		string name;
		string tr_id;
};
class PromoterClass{ //Probe Clusters Associated with a Promoter
	friend class NegCtrlClass;
	friend class ProbeSet;
public:
	vector <string> ChrNames_refseq;
	vector < vector <int> > refseq_indexes; //Based on refseq gene chromosome
	vector <string> ChrNames_proms;
	vector < vector <int> > prom_indexes; //Based on refseq gene chromosome
public:
GeneStruct *refseq;
PromoterStruct *proms;
temppars *tp;
void InitialiseData(void);
void ReadPromoterAnnotation(RESitesClass&, MappabilityClass&);
bool AnnotatewithPromoters(string,int,int,string,int,int,RESitesClass&,int);
void PrintInteractions(string);
int NofGenes;
int NofPromoters;
private:
	void GetPairInformation(stringstream&, stringstream&, PairStruct&);
	void GetTrFeats(stringstream&,temppars&);
	void ClusterIsoformPromoters(vector<int>,vector<string>,int,int);
	void DealwithSharedPromoters(int);
	void AnnotateDistalInteractor(string, string,int,int,int,int,int);
	void AnnotateFeatFeatInteraction(int, int, int);
	void PopulateInteractions(boost::unordered::unordered_map<int, int* >&, int, int,int);
};

void PromoterClass::InitialiseData(void){

	tp=new temppars [2];
	refseq = new GeneStruct [NumberofGenes];
	proms = new PromoterStruct [NumberofGenes];
	for(int i=0;i<NumberofGenes;i++){
		refseq[i].expression=(double *) calloc(3,sizeof(double));
		proms[i].expression=(double *) calloc(3,sizeof(double));
		proms[i].closestREsitenums = new int [2];
	}	

}
void PromoterClass::GetTrFeats(stringstream &trx, temppars &tpars){

	string field;
	getline(trx,tpars.name,'\t');

	getline(trx,tpars.tr_id,'\t'); // tr id
	getline(trx,tpars.chr,'\t');
	getline(trx,tpars.strand,'\t');
	getline(trx,field,'\t');
	if(tpars.strand=="+"){
		tpars.start=atoi(field.c_str());
		getline(trx,field,'\t');
		tpars.end=atoi(field.c_str());
	}
	else{
		tpars.end=atoi(field.c_str());
		getline(trx,field,'\t');
		tpars.start=atoi(field.c_str());
		}
	getline(trx,field,'\t');
	getline(trx,field,'\t');
	getline(trx,field,'\t');

	getline(trx,field,'\t');
	tpars.expression[0]=atof(field.c_str()); //mES
	getline(trx,field,'\t');
	tpars.expression[1]=atof(field.c_str()); //XEN
	getline(trx,field,'\t');
	tpars.expression[2]=atof(field.c_str()); //TS
	

}
void PromoterClass::ClusterIsoformPromoters(vector <int> isoformprs, vector<string> tr_ids, int gene_idx, int prom_idx){
unsigned int j,k,l,cluster_idx=0;
vector<int> clustercoords;
vector<string> clustered_tr_ids;
vector<bool> clustered;

for(j = 0;j < isoformprs.size(); ++j)
	clustered.push_back(0);

for(j = 0;j < isoformprs.size(); ++j){
	l=0;
	if(!clustered[j]){
		clustercoords.push_back(isoformprs[j]);
		clustered_tr_ids.push_back(tr_ids[j]);
		if(j+1<isoformprs.size()){
			for(k = j+1;k < isoformprs.size();++k){
				if(!clustered[j]){
					if(abs(isoformprs[j] - isoformprs[k]) < ClusterPromoters){
						++l;
						clustered[k] = 1;
						clustercoords.push_back(isoformprs[k]);
						clustered_tr_ids.push_back(tr_ids[k]);
					}
				}
			}
			refseq[gene_idx].isoformpromotercoords.push_back(clustercoords[0] + ((clustercoords[l] - clustercoords[0])/2));
			refseq[gene_idx].TranscriptName.push_back(clustered_tr_ids[0]);
			++cluster_idx;
			clustercoords.clear();
			clustered_tr_ids.clear();
		}
		else{
			refseq[gene_idx].isoformpromotercoords.push_back(isoformprs[j]);
			refseq[gene_idx].TranscriptName.push_back(tr_ids[j]);
		}
	}
}	
}
void PromoterClass::ReadPromoterAnnotation(RESitesClass& dpnIIsites, MappabilityClass& mapp)
{
	string temp,tr1,tr2;
	int geneindex = 0, promindex = 0;

#ifdef UNIX
string RefSeqfilename;
RefSeqfilename.append(dirname);
RefSeqfilename.append("mm9_refseq_sortedbyname_withexpression.txt");
ifstream RefSeq_file(RefSeqfilename.c_str());
#endif

#ifdef WINDOWS
	ifstream RefSeq_file("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\mm9_refseq_sortedbyname_withexpression.txt");
#endif

	vector < int > isoformprs; // to keep isoform promoters
	vector < string > tr_ids;
	getline(RefSeq_file,temp); //get the header row
	
	getline(RefSeq_file,tr1);
	stringstream trx1 ( tr1 ); 
	GetTrFeats(trx1,tp[0]);

	isoformprs.push_back(tp[0].start);
	tr_ids.push_back(tp[0].tr_id);
	do{
		getline(RefSeq_file,tr2);
		stringstream trx2 ( tr2);
		GetTrFeats(trx2,tp[1]);
		while(tp[0].name == tp[1].name){
			isoformprs.push_back(tp[1].start);
			tr_ids.push_back(tp[1].tr_id);
			getline(RefSeq_file,tr2);
			stringstream trx2 ( tr2);
			GetTrFeats(trx2,tp[1]);
		};
		if(isoformprs.size() > 1)
			ClusterIsoformPromoters(isoformprs,tr_ids,geneindex,promindex); //Promoters that are within "ClusterPromoters" of each other are clustered as one
		else{
			refseq[geneindex].isoformpromotercoords.push_back(isoformprs[0]);
			refseq[geneindex].TranscriptName.push_back(tp[0].tr_id);
		}
		refseq[geneindex].RefSeqName.append(tp[0].name);
		refseq[geneindex].chr.append(tp[0].chr);
		refseq[geneindex].end = tp[0].end;
		refseq[geneindex].strand = tp[0].strand;
		refseq[geneindex].expression[0] = tp[0].expression[0];
		refseq[geneindex].expression[1] = tp[0].expression[1];
		refseq[geneindex].expression[2] = tp[0].expression[2];
		
	// Fill promoter struct each unclustered isoform will have its promoter	
		for(int y = 0; y < refseq[geneindex].isoformpromotercoords.size();++y){
			proms[promindex].genes.push_back(refseq[geneindex].RefSeqName);
			proms[promindex].transcripts.push_back(refseq[geneindex].TranscriptName[y]);
			proms[promindex].strand = tp[0].strand;
			proms[promindex].chr.append(tp[0].chr);
			proms[promindex].TSS = refseq[geneindex].isoformpromotercoords[y];	
			proms[promindex].expression[0]=tp[0].expression[0];
			proms[promindex].expression[1]=tp[0].expression[1];
			proms[promindex].expression[2]=tp[0].expression[2];
			if( refseq[promindex].strand == "+" ){
				proms[promindex].start = refseq[geneindex].isoformpromotercoords[y] - coreprom_upstream;
				proms[promindex].end = refseq[geneindex].isoformpromotercoords[y] + coreprom_downstream;
			}
			else{
				proms[promindex].start = refseq[geneindex].isoformpromotercoords[y] - coreprom_downstream;
				proms[promindex].end = refseq[geneindex].isoformpromotercoords[y] + coreprom_upstream;								
			}
			//cout << proms[promindex].chr << "   "  << proms[promindex].start << "  " << proms[promindex].end << endl;
			proms[promindex].nofRESites = dpnIIsites.GetRESitesCount(proms[promindex].chr,proms[promindex].start,proms[promindex].end);
			dpnIIsites.GettheREPositions(proms[promindex].chr, proms[promindex].TSS, proms[promindex].closestREsitenums);
			proms[promindex].featmappability = mapp.GetMappability(proms[promindex].chr,proms[promindex].start, (coreprom_downstream + coreprom_upstream));
			++promindex;
		}

		isoformprs.clear();
		isoformprs.push_back(tp[1].start);
		tr_ids.clear();
		tr_ids.push_back(tp[1].tr_id);
//SWAP TP[0] and TP[1]
		tp[0].name.clear();
		tp[0].name.append(tp[1].name);
		tp[0].tr_id.clear();
		tp[0].tr_id.append(tp[1].tr_id);
		tp[0].chr.clear();
		tp[0].chr.append(tp[1].chr);
		tp[0].end = tp[1].end;
		tp[0].start = tp[1].start;
		tp[0].strand = tp[1].strand;
		tp[0].expression[0] = tp[1].expression[0];
		tp[0].expression[1] = tp[1].expression[1];
		tp[0].expression[2] = tp[1].expression[2];
		
		++geneindex;
		
		if(geneindex%10000 == 0){
		  cout << geneindex << "    Promoters Annotated" << endl;
		  //break;
		}
		//	cout << tp[1].name << endl;
	}while(tp[1].name!="END");

	NofGenes = geneindex;
	NofPromoters = promindex;
	DealwithSharedPromoters(NofGenes);
	cout << "Shared Promoters Determined" << endl;
	
	// Index Genes for faster access
	unsigned int i = 0;
	int found = 0;
	for(geneindex = 0;geneindex < NofGenes; ++geneindex){
		found = -1;
		for(i = 0;i < ChrNames_refseq.size(); ++i){
			if(refseq[geneindex].chr.compare(ChrNames_refseq[i]) == 0)
				found = i;
		}
		if(found == -1){
			ChrNames_refseq.push_back(refseq[geneindex].chr);
			refseq_indexes.push_back(vector<int>());
			refseq_indexes[ChrNames_refseq.size()-1].push_back(geneindex);
		}
		else
			refseq_indexes[found].push_back(geneindex);

	}
	cout << "RefSeq Chromosome Indexes Generated" << endl;
	// Index Promoters for faster access
	for(promindex = 0;promindex < NofPromoters; ++promindex){
		found=-1;
		for(i = 0; i < ChrNames_proms.size();++i){
			if(proms[promindex].chr.compare(ChrNames_proms[i]) == 0) 
				found = i;
		}
		if(found==-1){
			ChrNames_proms.push_back(proms[promindex].chr);
			prom_indexes.push_back(vector<int>());
			prom_indexes.back().push_back(promindex);
		}
		else
			prom_indexes[found].push_back(promindex);
		
	}
	cout << "Promoter Chromosome Indexes Generated" << endl;
}
void PromoterClass::DealwithSharedPromoters(int nofproms){ // If promoters are too close to each other
	
	for (int i = 0; i < nofproms; ++i) {
		proms[i].sharedpromoter = 0;
		for (int j = i + 1; j < nofproms - 1; ++j) {
			proms[j].sharedpromoter = 0;
			if (proms[i].chr.compare(proms[j].chr) == 0){
				if(((proms[i].start < proms[j].start) && (proms[i].end > proms[j].start)) || ((proms[i].start < proms[j].end) && (proms[i].end > proms[j].end))){
					proms[i].sharedpromoter = 1; 
					proms[j].sharedpromoter = 1;
					for( int k = 0; k<proms[i].genes.size();++k)
						proms[i].genes_sharingproms.push_back(proms[j].genes[k]);
					for( int k = 0; k<proms[j].genes.size();++k)
						proms[j].genes_sharingproms.push_back(proms[i].genes[k]);
				}
			}
		}
	}
}
bool PromoterClass::AnnotatewithPromoters(string p_chr_1, int resite_firstinpair, int pos_firstinpair, string p_chr_2, int resite_secondinpair, int pos_secondinpair,RESitesClass& dpnIIsites,int ExperimentNo){

int prom1_chrindex, prom2_chrindex;
bool pann = 0;
bool chrfound1 = 0, chrfound2 = 0;

p_chr_1.erase(p_chr_1.find_last_not_of(" \n\r\t")+1); //trim the string
for(int k = 0;k < ChrNames_proms.size(); ++k){
	if(ChrNames_proms[k].compare(p_chr_1.c_str()) == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
		prom1_chrindex = k;
		chrfound1 = 1;
		break;
	}
}
p_chr_2.erase(p_chr_2.find_last_not_of(" \n\r\t")+1);
for(int k = 0;k < ChrNames_proms.size(); ++k){
	if(ChrNames_proms[k].compare(p_chr_2.c_str()) == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
		prom2_chrindex = k;
		chrfound2 = 1;
		break;
	}
}

int pidx1 = -1, pidx2 = -1; // PromoterID
if(chrfound1){
	for(int i = 0; i < prom_indexes[prom1_chrindex].size(); ++i){ //Iterate over all refseq genes on that chromosome
		if((proms[prom_indexes[prom1_chrindex][i]].start <= pos_firstinpair && proms[prom_indexes[prom1_chrindex][i]].end >= pos_firstinpair)){ // If the readstart is contained within the core promoter
			pidx1 =	prom_indexes[prom1_chrindex][i];
			pann = 1; // Read is annotated with a promoter
			if(p_chr_1.compare(p_chr_2) == 0 && (proms[pidx1].start <= pos_secondinpair && proms[pidx1].end >= pos_secondinpair))
				return pann; // if the pair of the read is also contained within the core promoter
			if(chrfound2){	// Check if the pair is close to another promoter
				for(int m = 0; m < prom_indexes[prom2_chrindex].size(); ++m){ //Iterate over all refseq genes on that chromosome
					if((proms[prom_indexes[prom2_chrindex][m]].start <= pos_secondinpair && proms[prom_indexes[prom2_chrindex][m]].end >= pos_secondinpair)){ // If the readstart is contained within the core promoter
						pidx2 =	prom_indexes[prom2_chrindex][m]; // It is prom-prom interaction
						AnnotateFeatFeatInteraction(pidx1, pidx2,ExperimentNo);
						AnnotateFeatFeatInteraction(pidx2, pidx1,ExperimentNo);
						return pann;
					}
				}
			} // Checked if it is prom-prom interaction
			AnnotateDistalInteractor(p_chr_1,p_chr_2,resite_secondinpair,pos_firstinpair, pos_secondinpair, pidx1,ExperimentNo);
			return pann;
		}
	}
}// First read in the pair processed
// Annotate the second read in the pair (this could only be distal-prom interaction)
if(chrfound2){
	for(int m = 0; m < prom_indexes[prom2_chrindex].size(); ++m){ //Iterate over all refseq genes on that chromosome
		if((proms[prom_indexes[prom2_chrindex][m]].start <= pos_secondinpair && proms[prom_indexes[prom2_chrindex][m]].end >= pos_secondinpair)){ // If the readstart is contained within the core promoter
			pidx2 =	prom_indexes[prom2_chrindex][m];
		//	if (proms[pidx2].genes[0] == "0610005C13Rik")
		//		cout << proms[pidx2].genes[0] << '\t' << p_chr_2 << '\t' << p_chr_1 << '\t' << pos_secondinpair << '\t' << pos_firstinpair << '\t' << resite_secondinpair << '\t' << resite_firstinpair << endl;
			pann = 1;
			AnnotateDistalInteractor(p_chr_2,p_chr_1,resite_firstinpair,pos_secondinpair,pos_firstinpair,pidx2,ExperimentNo); 
			return pann;
		}			
	}
}
return pann;
}  
void PromoterClass::PopulateInteractions(boost::unordered::unordered_map<int, int* >& signals, int interactor_resite, int interactor_pos,int ExperimentNo){
	if(signals.find(interactor_resite) == signals.end()){
		signals[interactor_resite] = new int[NOFEXPERIMENTS + 1]; // add a new entry
		for(int z = 0; z < NOFEXPERIMENTS; ++z)
			signals[interactor_resite][z + 1] = 0;
		signals[interactor_resite][ExperimentNo + 1] = 1;
	}
	else{ // if inserted before	
		(signals[interactor_resite][ExperimentNo + 1]) = signals[interactor_resite][ExperimentNo + 1] + 1.0;
	}
	signals[interactor_resite][0] = interactor_pos; // Actual read pos
}
void PromoterClass::AnnotateDistalInteractor(string promoter_chr, string interactor_chr, int interactor_resite, int promoter_pos, int interactor_pos, int pidx, int ExperimentNo){
 
// Intra chromosomal interaction 
  if(promoter_chr.compare(interactor_chr) == 0){
	  if( (abs(interactor_pos - promoter_pos)) > MinimumJunctionDistance){ // at least minjunctdist away
		  if ( proms[pidx].strand == "+"){
			  if (interactor_pos < promoter_pos) // upstream 
				  PopulateInteractions(proms[pidx].Signals.signal_ups,interactor_resite,interactor_pos,ExperimentNo);
			  else // downstream
				  PopulateInteractions(proms[pidx].Signals.signal_down,interactor_resite,interactor_pos,ExperimentNo);
		  }
		  else{ // (-) strand
			  if (interactor_pos < promoter_pos) // downstream
				  PopulateInteractions(proms[pidx].Signals.signal_down,interactor_resite,interactor_pos,ExperimentNo);
			  else
				  PopulateInteractions(proms[pidx].Signals.signal_ups,interactor_resite,interactor_pos,ExperimentNo);
		  }
	  }
  }
  else{  // inter chromosomal interaction
	  bool chrfound = 0;
	  vector < SignalStruct_CTX >::iterator itx;
	  for(itx = proms[pidx].Signals_CTX.begin() ; itx < proms[pidx].Signals_CTX.end();++itx){
		  if (interactor_chr.compare(itx->maptochrname) == 0){
			  if(itx->signal.find(interactor_resite) == itx->signal.end()){
				  itx->signal[interactor_resite] = new int[NOFEXPERIMENTS + 1]; // add a new entry
				  for(int z = 0; z < NOFEXPERIMENTS; ++z)
					  itx->signal[interactor_resite][z + 1] = 0;
				  itx->signal[interactor_resite][ExperimentNo + 1] = 1;
				  itx->signal[interactor_resite][0] = interactor_pos; // Actual read pos
			  }
			  else{ // if inserted before
				  itx->signal[interactor_resite][ExperimentNo + 1] += 1;
				  itx->signal[interactor_resite][0] = interactor_pos; // Actual read pos
			  }
			  chrfound = 1;
			  break;
		  }
	  }
	  if(!chrfound){
		  proms[pidx].Signals_CTX.push_back(SignalStruct_CTX());
		  proms[pidx].Signals_CTX.back().maptochrname.append(interactor_chr);
		  proms[pidx].Signals_CTX.back().signal[interactor_resite] = new int[NOFEXPERIMENTS + 1];
		  for(int z = 0; z < NOFEXPERIMENTS; ++z)
			  proms[pidx].Signals_CTX.back().signal[interactor_resite][z + 1] = 0;
		  proms[pidx].Signals_CTX.back().signal[interactor_resite][ExperimentNo + 1] = 1;
		  proms[pidx].Signals_CTX.back().signal[interactor_resite][0] = interactor_pos; // Actual read pos
	  }
	}	 
}
void PromoterClass::AnnotateFeatFeatInteraction(int pidx1, int pidx2,int ExperimentNo){
	bool foundbefore = 0;	
	vector < FeattoFeatSignalStruct >::iterator it;
	for (it = proms[pidx1].promPromSignals.begin(); it < proms[pidx1].promPromSignals.end(); ++it){ // Check if the interaction with that promoter is seen before
		if (it->feat_index == pidx2){
			it->normsignal[ExperimentNo] += 1.0;
			foundbefore = 1;
			break;
		}
	}
	if(!foundbefore){ // Create a new entry for that promoter
		proms[pidx1].promPromSignals.push_back(FeattoFeatSignalStruct());
		proms[pidx1].promPromSignals.back().feat_index = pidx2;
		proms[pidx1].promPromSignals.back().normsignal[ExperimentNo] += 1.0;
	}
}
/*
void PromoterClass::AssociateProbeswithPromoters(ProbeSet& MM9prs){

int i=0,j=0,k=0,x=0;
int diff;
long int startsearch,endsearch;

for(i=0;i<NofPromoters;++i){
	startsearch=0;endsearch=-1;
	for(j=0;j<MM9prs.ChrNames.size();++j){
		if(refseq[i].chr==MM9prs.ChrNames[j]){
			startsearch=MM9prs.ChrRowStartIndexes[j];
			endsearch=MM9prs.ChrRowEndIndexes[j];
			break;
		}
	}
	x=0;
	for(j=startsearch;j<=endsearch;++j){
		for(k=0;k<refseq[i].isoformpromotercoords.size();++k){
			diff=MM9prs.MM9Probes[j].center-refseq[i].isoformpromotercoords[k];
			if(abs(diff)<AssociateInteractions){
					refseq[i].ProbeIDs.push_back(j);
			}
		}
		if(diff>BinSize)
			break;
	}
	if(i%1000==0)
		cout << i << "    Probes Associated with Promoters" << endl;
}

}
*/

