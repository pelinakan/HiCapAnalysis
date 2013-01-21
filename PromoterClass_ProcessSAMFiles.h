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
GeneStruct *refseq;
PromoterStruct *proms;
temppars *tp;
void InitialiseData(void);
void ReadPromoterAnnotation(RESitesClass&, MappabilityClass&);
bool AnnotatewithPromoters(string,int,string,int,RESitesClass&, MappabilityClass&);
void PrintInteractions(string);
int NofGenes;
int NofPromoters;
private:
	void GetPairInformation(stringstream&, stringstream&, PairStruct&);
	void GetTrFeats(stringstream&,temppars&);
	void ClusterIsoformPromoters(vector<int>,vector<string>,int,int);
	void DealwithSharedPromoters(int);
	void AnnotateDistalInteractor(string, string,int*,int);
	void AnnotateFeatFeatInteraction(int, int);
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
		
		if(geneindex%10000 == 0)
		  cout << geneindex << "    Promoters Annotated" << endl;
		//	cout << tp[1].name << endl;
	}while(!(RefSeq_file.eof()));

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
				found=i;
		}
		if(found==-1){
			ChrNames_proms.push_back(proms[promindex].chr);
			prom_indexes.push_back(vector<int>());
			prom_indexes[ChrNames_proms.size()-1].push_back(promindex);
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
bool PromoterClass::AnnotatewithPromoters(string p_chr_1, int pstcoord, string p_chr_2, int pendcoord, RESitesClass& dpnIIsites, MappabilityClass& mappability){

int prom1_chrindex, prom2_chrindex;
bool pann = 0;
bool chrfound = 0;

p_chr_1.erase(p_chr_1.find_last_not_of(" \n\r\t")+1); //trim the string
for(int k = 0;k < ChrNames_proms.size(); ++k){
	if(ChrNames_proms[k].compare(p_chr_1.c_str()) == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
		prom1_chrindex = k;
		chrfound = 1;
		break;
	}
}
if (!chrfound)
	return 0;
chrfound = 0;
p_chr_2.erase(p_chr_2.find_last_not_of(" \n\r\t")+1);
for(int k = 0;k < ChrNames_proms.size(); ++k){
	if(ChrNames_proms[k].compare(p_chr_2.c_str()) == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
		prom2_chrindex = k;
		chrfound = 1;
		break;
	}
}
if (!chrfound)
	return 0;

int pidx1 = -1, pidx2 = -1; // PromoterID
for(int i = 0; i < prom_indexes[prom1_chrindex].size(); ++i){ //Iterate over all refseq genes on that chromosome
	if((proms[prom_indexes[prom1_chrindex][i]].start <= pstcoord && proms[prom_indexes[prom1_chrindex][i]].end >= pstcoord)){ // If the readstart is contained within the core promoter
		pidx1 =	prom_indexes[prom1_chrindex][i];
		pann = 1; // Read is annotated with a promoter
		if(p_chr_1.compare(p_chr_2) == 0 && (proms[pidx1].start <= pendcoord && proms[pidx1].end >= pendcoord))
		  return pann; // if the pair of the read is also contained within the core promoter
		// Check if the pair is close to a promoter
		for(int m = 0; m < prom_indexes[prom2_chrindex].size(); ++m){ //Iterate over all refseq genes on that chromosome
			if((proms[prom_indexes[prom2_chrindex][m]].start <= pendcoord && proms[prom_indexes[prom2_chrindex][m]].end >= pendcoord)){ // If the readstart is contained within the core promoter
				pidx2 =	prom_indexes[prom2_chrindex][m]; // It is prom-prom interaction
				AnnotateFeatFeatInteraction(pidx1, pidx2);
				AnnotateFeatFeatInteraction(pidx2, pidx1);
				return pann;
			}
		} // Checked if it is prom-prom interaction
		int *renums;
		renums = new int [2];
		dpnIIsites.GettheREPositions(p_chr_2,pendcoord,renums); // Get the fragment number of the distal interactor
		AnnotateDistalInteractor(p_chr_1,p_chr_2,renums,pidx1);
		return pann;
	}
}// First read in the pair processed
// Annotate the second read in the pair (this could only be distal-prom interaction)
for(int m = 0; m < prom_indexes[prom2_chrindex].size(); ++m){ //Iterate over all refseq genes on that chromosome
	if((proms[prom_indexes[prom2_chrindex][m]].start <= pendcoord && proms[prom_indexes[prom2_chrindex][m]].end >= pendcoord)){ // If the readstart is contained within the core promoter
		pann = 1;
		pidx2 =	prom_indexes[prom2_chrindex][m];
		if(p_chr_1.compare(p_chr_2) == 0 && (proms[pidx2].start <= pstcoord && proms[pidx2].end >= pstcoord))
			return pann; // if the pair of the read is also contained within the core promoter
		int *renums;
		renums = new int [2];
		dpnIIsites.GettheREPositions(p_chr_1,pstcoord,renums); // Get the fragment number of the distal interactor (this is it is the first read)
		AnnotateDistalInteractor(p_chr_2,p_chr_1,renums,pidx2); 
		return pann;
	}			
}
return pann;
}  
void PromoterClass::AnnotateDistalInteractor(string p_chr_1, string p_chr_2, int *renums, int pidx){
 
// Intra chromosomal interaction 
  if(p_chr_1.compare(p_chr_2) == 0){ 
	  if ( proms[pidx].strand == "+"){
		  if (renums[1] < proms[pidx].closestREsitenums[0]){ // upstream 
			  if( (abs(renums[1] - proms[pidx].closestREsitenums[0])) > MinimumJunctionDistance){ // at least minjunctdist away
				  if(proms[pidx].Signals.signal_ups.find(renums[1]) == proms[pidx].Signals.signal_ups.end())
					  proms[pidx].Signals.signal_ups[renums[1]] = 1; // add a new entry
				  else // if inserted before
					  proms[pidx].Signals.signal_ups[renums[1]] = proms[pidx].Signals.signal_ups[renums[1]] + 1;
			  }
		  }
		  if (renums[0] > proms[pidx].closestREsitenums[1]){ // downstream
			  if ((abs(renums[0] - proms[pidx].closestREsitenums[1])) > MinimumJunctionDistance ){ // at least minjunctdist away
				  if(proms[pidx].Signals.signal_down.find(renums[0]) == proms[pidx].Signals.signal_down.end())
					  proms[pidx].Signals.signal_down[renums[0]] = 1; // add a new entry
				  else // if inserted before
					  proms[pidx].Signals.signal_down[renums[0]] = proms[pidx].Signals.signal_down[renums[0]] + 1;
			  }
		  }
	  }
	  else{ // (-) strand
		  if (renums[1] < proms[pidx].closestREsitenums[0]){ // downstream
			  if( (abs(renums[1] - proms[pidx].closestREsitenums[0])) > MinimumJunctionDistance){ // at least minjunctdist away
				  if(proms[pidx].Signals.signal_down.find(renums[1]) == proms[pidx].Signals.signal_down.end())
					  proms[pidx].Signals.signal_down[renums[1]] = 1; // add a new entry
				  else // if inserted before
					  proms[pidx].Signals.signal_down[renums[1]] = proms[pidx].Signals.signal_down[renums[1]] + 1;
			  }
		  }
		  if (renums[0] > proms[pidx].closestREsitenums[1]){ // upstream
			  if ((abs(renums[0] - proms[pidx].closestREsitenums[1])) > MinimumJunctionDistance ){ // at least minjunctdist away
				  if(proms[pidx].Signals.signal_ups.find(renums[0]) == proms[pidx].Signals.signal_ups.end())
					  proms[pidx].Signals.signal_ups[renums[0]] = 1; // add a new entry
				  else // if inserted before
					  proms[pidx].Signals.signal_ups[renums[0]] = proms[pidx].Signals.signal_ups[renums[0]] + 1;
			  }
		  }
	  }
  }
  
  else{  // inter chromosomal interaction
	  bool chrfound = 0;
	  vector < SignalStruct_CTX >::iterator it;
	  for(it = proms[pidx].Signals_CTX.begin() ; it < proms[pidx].Signals_CTX.end();++it){
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
		  proms[pidx].Signals_CTX.push_back(SignalStruct_CTX());
		  proms[pidx].Signals_CTX.back().maptochrname.append(p_chr_2);
		  proms[pidx].Signals_CTX.back().signal[renums[0]] = 1;
	  }
	}	 
 
}
void PromoterClass::AnnotateFeatFeatInteraction(int pidx1, int pidx2){
	bool foundbefore = 0;	
	vector < FeattoFeatSignalStruct >::iterator it;
	for (it = proms[pidx1].promPromSignals.begin(); it < proms[pidx1].promPromSignals.end(); ++it){ // Check if the interaction with that promoter is seen before
		if (it->feat_index == pidx2){
			it->normsignal += 1.0;
			foundbefore = 1;
			break;
		}
	}
	if(!foundbefore){ // Create a new entry for that promoter
		proms[pidx1].promPromSignals.push_back(FeattoFeatSignalStruct());
		proms[pidx1].promPromSignals.back().feat_index = pidx2;
		proms[pidx1].promPromSignals.back().normsignal += 1.0;
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
void PromoterClass::PrintInteractions(string BaseFileName){

	string FileName;
	FileName.append(BaseFileName);
	FileName.append("MM9.RefSeqProms.UpsteamSignals.txt");
	ofstream outf1(FileName.c_str());

	boost::unordered::unordered_map< int, int >::const_iterator iter;
	for(int i = 0; i < NofPromoters; ++i){
		outf1 << proms[i].chr << ":" << proms[i].start << "-" << proms[i].end << '\t' << proms[i].nofRESites << '\t' << proms[i].interactions.size() << '\t';
		for(iter = proms[i].Signals.signal_ups.begin(); iter != proms[i].Signals.signal_ups.end(); ++iter){
			if(iter->second > MinNumberofReads)
			outf1 << (proms[i].closestREsitenums[0] - iter->first ) << " , " << iter->second << '\t';
		}
		outf1 << endl;
	}
	outf1.close();

	FileName.clear();
	FileName.append(BaseFileName);
	FileName.append("MM9.RefSeqProms.DownsteamSignals.txt");
	ofstream outf2(FileName.c_str());

	for(int i = 0; i < NofPromoters; ++i){
		outf2 << proms[i].chr << ":" << proms[i].start << "-" << proms[i].end << '\t' << proms[i].interactions.size() << '\t';
		for( iter = proms[i].Signals.signal_down.begin(); iter != proms[i].Signals.signal_down.end(); ++iter){
			if(iter->second > MinNumberofReads)
				outf2 << (iter->first - proms[i].closestREsitenums[1]) << " , " << iter->second << '\t';
		}
		outf2 << endl;
	}
	outf2.close();
}
