class DetectEnrichedBins{
public:
	void DetectEnrichedInteractionBins(PromoterClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins(PromoterClass&,string,RESitesClass&);
	void DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass&,string,RESitesClass&);
	void PrintAllInteractions(PromoterClass&,NegCtrlClass&,string, int,vector < string >);
	void PrintNearestPromoter(PromoterClass&,int);
private:
	bool CheckSupportingPairs(int*);
	double CalculatepVal(boost::unordered::unordered_map< int, double >,boost::unordered::unordered_map< int, double >, int,int);
	void FillInteractionStruct(vector<InteractionStruct>&,string, int, int, int, int*,double, char);
};
bool DetectEnrichedBins::CheckSupportingPairs(int* supppairs){
	bool recordit = 0;
	for (int t = 0; t < NOFEXPERIMENTS; ++t){
		if (supppairs[t + 1] > MinNumberofReads){
			recordit = 1;
			break;
		}
	}
	return recordit;
}
double DetectEnrichedBins::CalculatepVal(boost::unordered::unordered_map< int, double > mean,boost::unordered::unordered_map< int, double > std, int bin,int sig){
	double q = 0, p_value = 1;
	if(mean.find(bin) == mean.end()){
		mean[bin] = 0;
		std[bin] = 1;
		p_value = 0.0;
	}
	else{
		q = alglib::normaldistribution(( sig - mean[bin]) / std[bin]); //z-score
		p_value = 1- q;
	}
	return p_value;
}

void DetectEnrichedBins::FillInteractionStruct(vector<InteractionStruct>& interactionstr,string chr, int dist, int ExperimentNo, int resite, int* values, double p_value, char inttype){
	vector<InteractionStruct>::iterator interactorit;
	bool addedbefore = 0;
	for(interactorit = interactionstr.begin(); interactorit < interactionstr.end(); ++interactorit){
		if(interactorit->resites[0] == resite){ //Added before
			addedbefore = 1;
			interactorit->supp_pairs[ExperimentNo] = values[ExperimentNo + 1];
			interactorit->p_val[ExperimentNo] = p_value;						
		}
	}
	if(!addedbefore){
		interactionstr.push_back(InteractionStruct());
		interactionstr.back().peakprofile.clear();
		interactionstr.back().p_val[ExperimentNo] = p_value;
		interactionstr.back().supp_pairs[ExperimentNo] = values[ExperimentNo + 1];
		interactionstr.back().chr = chr;
		interactionstr.back().distance = dist; 
		interactionstr.back().resites[0] = resite;
		interactionstr.back().pos = values[0];
		interactionstr.back().type = inttype; // U: upstream , D: downstream , X: inter-chromosomal
	}

} 
void DetectEnrichedBins::DetectEnrichedInteractionBins(PromoterClass &prs, DetermineBackgroundLevels background, string BaseFileName,int ExperimentNo){
//	MappabilityClass mapp;
//	mapp.InitialiseVars();

	for(int i = 0; i < prs.NofPromoters; ++i){
		boost::unordered::unordered_map<int, int* >::const_iterator it; //key:REpos, value[0]:genomic position 
		for (it = prs.proms[i].Signals.signal_ups.begin(); it != prs.proms[i].Signals.signal_ups.end(); ++it){
			bool enoughpairs = 0;
			int dist = 0;
			enoughpairs = CheckSupportingPairs(it->second);
			if (enoughpairs){
				if (prs.proms[i].strand =="+")
					dist = it->second[0] - prs.proms[i].closestREsitenums[0]; // negative
				else
					dist = prs.proms[i].closestREsitenums[1] - it->second[0]; // negative
				int bin = abs(dist) / BinSize; 
				double p_value = CalculatepVal(background.bglevels.mean_upstream,background.bglevels.stdev_upstream, bin,it->second[ExperimentNo + 1]);
				FillInteractionStruct(prs.proms[i].interactions,prs.proms[i].chr, dist, ExperimentNo, it->first, it->second, p_value,'U');
			}
		}
	//Downstream
		for (it = prs.proms[i].Signals.signal_down.begin(); it != prs.proms[i].Signals.signal_down.end(); ++it){
			bool enoughpairs = 0;
			int dist = 0;
			enoughpairs = CheckSupportingPairs(it->second);
			if (enoughpairs){
					if (prs.proms[i].strand =="+")
						dist = it->second[0] - prs.proms[i].closestREsitenums[1]; // positive
					else
						dist = prs.proms[i].closestREsitenums[0] - it->second[0]; // positive
					int bin = abs(dist) / BinSize;
					double p_value = CalculatepVal(background.bglevels.mean_downstream,background.bglevels.stdev_downstream, bin,it->second[ExperimentNo + 1]);
					FillInteractionStruct(prs.proms[i].interactions,prs.proms[i].chr, dist, ExperimentNo, it->first, it->second, p_value,'D');
			}
		}
	// Inter-chromosomal
		for( int j = 0; j < prs.proms[i].Signals_CTX.size(); ++j){ // For each chromosome
			for (it = prs.proms[i].Signals_CTX[j].signal.begin(); it != prs.proms[i].Signals_CTX[j].signal.end(); ++it){
				bool enoughpairs = 0;
				int dist = 0;
				enoughpairs = CheckSupportingPairs(it->second);
				if (enoughpairs){
					FillInteractionStruct(prs.proms[i].interactions,prs.proms[i].Signals_CTX[j].maptochrname, -1, ExperimentNo, it->first, it->second, -1,'X');
				}
			}
		}
	}
}
void DetectEnrichedBins::AssociatePeaksWithIntBins(PromoterClass& prs,string BaseFileName,RESitesClass& dpnII){

	cout << "Associating with peaks" << endl;
	for(int i = 0; i < prs.NofPromoters; ++i){
		for(int j = 0; j < prs.proms[i].interactions.size(); ++j){
			bool peakfound = 0, refound = 0;
			int *renums;
			renums = new int [2];
			refound = dpnII.GettheREPositions(prs.proms[i].chr,prs.proms[i].interactions[j].resites[0],renums); // Interactor RE fragment
			if(refound)
				prs.proms[i].interactions[j].resites[1] = renums[1];
			else
				prs.proms[i].interactions[j].peakoverlap = false;
			boost::unordered::unordered_map<string, int>::const_iterator citer = MetaPeakChrMap.find(prs.proms[i].interactions[j].chr);
			if(refound && !(MetaPeakChrMap.find(prs.proms[i].interactions[j].chr) == MetaPeakChrMap.end()) ){ // Find the right peak map wrt chromosome
				if(metaPeaks[citer->second].metapeaks.find(prs.proms[i].interactions[j].resites[0]) != metaPeaks[citer->second].metapeaks.end()){
					//If there is peak
					boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(prs.proms[i].interactions[j].resites[0]);
					prs.proms[i].interactions[j].peakoverlap = true;
					prs.proms[i].interactions[j].peakprofile = it->second;
					peakfound = 1;
				}
				else{
					int *repos;
					repos = new int [2];
					refound = dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(prs.proms[i].interactions[j].resites[0] - 1),repos);
					if(refound){
						do{
							if(metaPeaks[citer->second].metapeaks.find(repos[0]) != metaPeaks[citer->second].metapeaks.end()){
								//If there is peak
								boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[0]);
								prs.proms[i].interactions[j].peakoverlap = true;
								prs.proms[i].interactions[j].peakprofile = it->second;
								peakfound = 1;
								break;
							}
							else{
								refound = dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(repos[0] - 1),repos);
								if(!refound){
									prs.proms[i].interactions[j].peakoverlap = false;
									break;
								}
							}
						}while((abs(repos[0] - prs.proms[i].interactions[j].resites[0]) < 1000));
					}
					else
						prs.proms[i].interactions[j].peakoverlap = false;
					
					if(!peakfound){
						refound = dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(prs.proms[i].interactions[j].resites[0]),repos);
						if(refound){
							do{
								if(metaPeaks[citer->second].metapeaks.find(repos[1]) != metaPeaks[citer->second].metapeaks.end()){
									//If there is peak
									boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[1]);
									prs.proms[i].interactions[j].peakoverlap = true;
									prs.proms[i].interactions[j].peakprofile = it->second;
									peakfound = 1;
									break;
								}
								else{
									refound = dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(repos[1]),repos);
									if(!refound){
										prs.proms[i].interactions[j].peakoverlap = false;
										break;
									}
								}
							}while ( abs(repos[1] - prs.proms[i].interactions[j].resites[0]) < 1000 );
						}
						else
							prs.proms[i].interactions[j].peakoverlap = false;
					}
					if (!peakfound)
						prs.proms[i].interactions[j].peakoverlap = false;
				}
			}
		}
	}
}
void DetectEnrichedBins::DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass &ngs,DetermineBackgroundLevels background,string BaseFileName,int ExperimentNo){
// Upstream
	for(int i = 0; i < ngs.NofNegCtrls; ++i){
//		outf1 << ngs.negctrls[i].chr << ":" << ngs.negctrls[i].start << "-" << ngs.negctrls[i].end << '\t';
		boost::unordered::unordered_map<int, int* >::const_iterator it;
		for (it = ngs.negctrls[i].Signals.signal_ups.begin(); it != ngs.negctrls[i].Signals.signal_ups.end(); ++it){
			bool enoughpairs = 0;
			enoughpairs = CheckSupportingPairs(it->second);
			if (enoughpairs){
				double dist = ngs.negctrls[i].closestREsitenums[0] - it->first;
				FillInteractionStruct(ngs.negctrls[i].interactions,ngs.negctrls[i].chr, dist, ExperimentNo, it->first, it->second, -1,'U');
			}
		}
//Downstream
		for (it = ngs.negctrls[i].Signals.signal_down.begin(); it != ngs.negctrls[i].Signals.signal_down.end(); ++it){
			bool enoughpairs = 0;
			enoughpairs = CheckSupportingPairs(it->second);
			if (enoughpairs){
				double dist = ngs.negctrls[i].closestREsitenums[0] - it->first;
				FillInteractionStruct(ngs.negctrls[i].interactions,ngs.negctrls[i].chr, dist, ExperimentNo, it->first, it->second, -1,'D');
			}
		}
// Inter-chromosomal
		for( int j = 0; j < ngs.negctrls[i].Signals_CTX.size(); ++j){ // For each chromosome
			for (it = ngs.negctrls[i].Signals_CTX[j].signal.begin(); it != ngs.negctrls[i].Signals_CTX[j].signal.end(); ++it){
				bool enoughpairs = 0;
				enoughpairs = CheckSupportingPairs(it->second);
				if (enoughpairs){
					double dist = ngs.negctrls[i].closestREsitenums[0] - it->first;
					FillInteractionStruct(ngs.negctrls[i].interactions,ngs.negctrls[i].Signals_CTX[j].maptochrname, -1, ExperimentNo, it->first, it->second, -1,'X');
				}
			}
		}
	}
}
void DetectEnrichedBins::AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass& ngs,string BaseFileName,RESitesClass& dpnII){
	
//cout << "Associating with peaks negative controls" << endl;
for(int i = 0; i < ngs.NofNegCtrls; ++i){
//	cout << i << '\t' << ngs.negctrls[i].interactions.size() << '\t';
	for(int j = 1; j < ngs.negctrls[i].interactions.size(); ++j){
		int *renums;
		bool peakfound = 0, refound;
		renums = new int [2];
		refound = dpnII.GettheREPositions(ngs.negctrls[i].chr,ngs.negctrls[i].interactions[j].resites[0],renums); // Interactor RE fragment
		if(refound)
			ngs.negctrls[i].interactions[j].resites[1] = renums[1];				
		else
			ngs.negctrls[i].interactions[j].peakoverlap = "NoPeakOverlap";				
		boost::unordered::unordered_map<string, int>::const_iterator citer = MetaPeakChrMap.find(ngs.negctrls[i].chr);
		if (refound && !(MetaPeakChrMap.find(ngs.negctrls[i].chr) == MetaPeakChrMap.end())){
			if(metaPeaks[citer->second].metapeaks.find(ngs.negctrls[i].interactions[j].resites[0]) != metaPeaks[citer->second].metapeaks.end()){
				//If there is peak
				boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(ngs.negctrls[i].interactions[j].resites[0]);
				ngs.negctrls[i].interactions[j].peakoverlap = true;
				ngs.negctrls[i].interactions[j].peakprofile = it->second;
				peakfound = 1;
			}
			else{
				int *repos;
				repos = new int [2];
				refound = dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,(ngs.negctrls[i].interactions[j].resites[0] - 1),repos);
				if(refound){
					do{
						if(metaPeaks[citer->second].metapeaks.find(repos[0]) != metaPeaks[citer->second].metapeaks.end()){
							//If there is peak
							boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[0]);
							ngs.negctrls[i].interactions[j].peakoverlap = true;
							ngs.negctrls[i].interactions[j].peakprofile = it->second;
							peakfound = 1;
							break;
						}
						else{
							refound = dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,(repos[0] - 1),repos);
							if(!refound){
								ngs.negctrls[i].interactions[j].peakoverlap = false;
								break;
							}
						}
					}while((abs(repos[0] - ngs.negctrls[i].interactions[j].resites[0]) < 1000));
				}
				else
					ngs.negctrls[i].interactions[j].peakoverlap = false;
				
				if(!peakfound){
					refound = dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,ngs.negctrls[i].interactions[j].resites[0],repos);
					if(refound){
						do{
							if(metaPeaks[citer->second].metapeaks.find(repos[1]) != metaPeaks[citer->second].metapeaks.end()){
								//If there is peak
								boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[1]);
								ngs.negctrls[i].interactions[j].peakoverlap = true;
								ngs.negctrls[i].interactions[j].peakprofile = it->second;
								peakfound = 1;
								break;
							}
							else{
								refound = dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,(repos[1]),repos);
								if(!refound){
									ngs.negctrls[i].interactions[j].peakoverlap = false;
									break;					
								}
							}
						}while ( abs(repos[1] - ngs.negctrls[i].interactions[j].resites[0]) < 1000 );
					}
					else
						ngs.negctrls[i].interactions[j].peakoverlap = false;
				}
				if (!peakfound)
					ngs.negctrls[i].interactions[j].peakoverlap = false;
			}			
		}
	}
//	cout << "done" << endl;
}
}

void DetectEnrichedBins::PrintAllInteractions(PromoterClass& prs, NegCtrlClass& ngs, string BaseFileName, int NumberofExperiments,vector < string > ExperimentNames){
//Promoter Interactions
	string FileName, FileName5;
	FileName.append(BaseFileName);
	FileName.append("SignificantInteractions");
	FileName.append("_Promoters");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());
	int s = 0;

	FileName5.append(BaseFileName);
	FileName5.append("NearestPromoterDistance.txt");
	ofstream outf5(FileName5.c_str());

	int dist;
	outf1 << "Gene Name" << '\t' << "Representative Transcript Name" << '\t' << "Gene Expression" << '\t' << "Shared Promoter" << '\t' << "Number of RE Sites within core promoter" << '\t' << "Mappability of Core Promoter" << '\t' 
		<< "Promoter chr" << '\t' << "Promoter TSS"  << '\t' << "Strand" << '\t' <<"Promoter UCSC" << '\t' << "Interactor chr" << '\t' << "Interactor REsite1" << '\t' << "Interactor REsite2" << '\t'
		<< "Read Position" << '\t' << "Interactor Mappability (100 bases around the RE site)" << '\t' << "Distance" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf1 << ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_pval" << '\t';
	outf1 << "Peak Profile" << endl;

	for(int i = 0; i < prs.NofPromoters; ++i){
		for(int j = 0; j < prs.proms[i].interactions.size(); ++j){
			string nearestgenename;
			int nearestgeneTSS;
			outf1 << prs.proms[i].genes[0] << '\t' << prs.proms[i].transcripts[0] << '\t' << prs.proms[i].expression[0] << '\t';
			outf5 << prs.proms[i].genes[0] << '\t' << prs.proms[i].transcripts[0] << '\t' << prs.proms[i].expression[0] << '\t' << prs.proms[i].TSS << '\t';
//			if(prs.proms[i].chr == prs.proms[i].interactions[j].chr){
				dist = prs.FindNearestPromoter(prs.proms[i].interactions[j].chr, prs.proms[i].interactions[j].pos,nearestgenename,nearestgeneTSS);
				outf5 << prs.proms[i].interactions[j].chr << '\t' << prs.proms[i].interactions[j].distance  << '\t' << nearestgenename << '\t' << dist << '\t' << nearestgeneTSS << '\t';
				for (int e = 0; e < NumberofExperiments; ++e)
					outf5 << prs.proms[i].interactions[j].supp_pairs[e] << '\t' << prs.proms[i].interactions[j].p_val[e] << '\t';
				outf5 << endl;
//			}
//			else outf5 << endl;

			if (prs.proms[i].sharedpromoter){
				for (s = 0; s < prs.proms[i].genes_sharingproms.size() - 1; ++s)
					outf1 << prs.proms[i].genes_sharingproms[s] << ',';
				outf1 << prs.proms[i].genes_sharingproms[s];
			}
			else
				outf1 << prs.proms[i].sharedpromoter;
			outf1 << '\t' << prs.proms[i].nofRESites << '\t' << prs.proms[i].featmappability << '\t' 
				<< prs.proms[i].chr << '\t' << prs.proms[i].TSS << '\t' << prs.proms[i].strand << '\t' << prs.proms[i].chr << ":" << prs.proms[i].start << "-" << prs.proms[i].end << '\t'; // Promoter Details
			outf1 << prs.proms[i].interactions[j].chr << '\t' << prs.proms[i].interactions[j].resites[0] << '\t' << prs.proms[i].interactions[j].resites[1] << '\t';
			outf1 << prs.proms[i].interactions[j].pos << '\t';
			outf1 << prs.proms[i].interactions[j].mappability << '\t' << prs.proms[i].interactions[j].distance << '\t';
			for (int e = 0; e < NumberofExperiments; ++e)
				outf1 << prs.proms[i].interactions[j].supp_pairs[e] << '\t' << prs.proms[i].interactions[j].p_val[e] << '\t';
			if(prs.proms[i].interactions[j].peakoverlap)
				outf1 << prs.proms[i].interactions[j].peakprofile << endl;
			else
				outf1 << "NoPeakOverlap" << endl;
		}
	}
	outf1.close();
//Negative Control Interactions
	string FileName2;
	FileName2.append(BaseFileName);
	FileName2.append("SignificantInteractions");
	FileName2.append("_NegCtrls");
	FileName2.append(".txt");
	ofstream outf2(FileName2.c_str());

	outf2 << "NegCtrl ID" << '\t' << "NegCtrl UCSC" << '\t' << "Number of RE sites within core feature" << '\t' << "Interactor Chr " << '\t' << "Interactor RE site1" << '\t' << "Interactor RE site2" << '\t'
		  << "Read Position" << '\t'  << "Distance" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf2 << ExperimentNames[e] << "_SuppPairs" << '\t';
	outf2 << "Peak Profile" << endl;

	for(int i = 0; i < ngs.NofNegCtrls; ++i){
		for(int j = 1; j < ngs.negctrls[i].interactions.size(); ++j){
			outf2 << i << '\t' << ngs.negctrls[i].chr << ":" << ngs.negctrls[i].start << "-" << ngs.negctrls[i].end << '\t' << ngs.negctrls[i].nofRESites << '\t'; // Feat Details
			outf2 << ngs.negctrls[i].interactions[j].chr << '\t' << ngs.negctrls[i].interactions[j].resites[0] << '\t' << ngs.negctrls[i].interactions[j].resites[1] << '\t';
			outf2 << ngs.negctrls[i].interactions[j].pos << '\t';
			outf2 << ngs.negctrls[i].interactions[j].distance << '\t';
			for (int e = 0; e < NumberofExperiments; ++e)
				outf2 << ngs.negctrls[i].interactions[j].supp_pairs[e] << '\t';			
			if(ngs.negctrls[i].interactions[j].peakoverlap)
				outf2 << ngs.negctrls[i].interactions[j].peakprofile << endl;
			else
				outf2 << "NoPeakOverlap" << endl;
		}
	}
	outf2.close();
//Promoter-promoter Interactions
	string FileName3;
	FileName3.append(BaseFileName);
	FileName3.append("Interactions_PromProm");
	FileName3.append(".txt");
	ofstream outf3(FileName3.c_str());

	outf3 << "Gene1 name" << '\t' << "Gene1 promoter UCSC" << '\t' << "Gene1 chr" << '\t' << "Gene1 TSS" << '\t' << "Gene1 expression" << '\t'
		<< "Gene1 number of RE Sites within core promoter" << '\t' << "Gene1 core promoter mappability" << '\t' << "Promoter Shared" << '\t'
		<< "Gene2 name" << '\t' << "Gene2 Promoter UCSC" << '\t' << "Gene2 chr" << '\t' << "Gene2 TSS" << '\t' << "Gene2 expression" << '\t'
		<< "Gene2 number of RE Sites within core promoter" << '\t' << "Gene2 core promoter mappability" << '\t' << "Promoter Shared" << '\t';
	for(int e=0;e < NumberofExperiments;++e)
		outf3 << ExperimentNames[e] << "_Supp_Pairs" << '\t';
	outf3 << endl;

	for(int i = 0; i < prs.NofPromoters; ++i){
		int s = 0;
		vector< FeattoFeatSignalStruct >::const_iterator it; //first: REpos, second: signal 
		for (it = prs.proms[i].promPromSignals.begin(); it != prs.proms[i].promPromSignals.end(); ++it){
			//double normsignal = (double(it->normsignal)) / prs.proms[i].nofRESites;
			int pindex = it->feat_index;
			outf3 << prs.proms[i].genes[0] << '\t' << prs.proms[i].chr << ":" << prs.proms[i].start << "-" << prs.proms[i].end << '\t'
				<< prs.proms[i].chr << '\t' << prs.proms[i].TSS << '\t' << prs.proms[i].expression[0] << '\t'
				<< prs.proms[i].nofRESites << '\t' << prs.proms[i].featmappability << '\t'; 
			if (prs.proms[i].sharedpromoter){
				for (s = 0; s < prs.proms[i].genes_sharingproms.size() - 1; ++s)
					outf3 << prs.proms[i].genes_sharingproms[s] << ',';
				outf3 << prs.proms[i].genes_sharingproms[s];
			}
			else
				outf3 << prs.proms[i].sharedpromoter;						
			outf3 <<  '\t' << prs.proms[pindex].genes[0] << '\t' << prs.proms[pindex].chr << ":" << prs.proms[pindex].start << "-" << prs.proms[pindex].end << '\t'
				<< prs.proms[pindex].chr << '\t' << prs.proms[pindex].TSS << '\t' << prs.proms[pindex].expression[0] << '\t'
				<< prs.proms[pindex].nofRESites << '\t' << prs.proms[pindex].featmappability << '\t';
			if (prs.proms[pindex].sharedpromoter){
				for (s = 0; s < prs.proms[pindex].genes_sharingproms.size() - 1; ++s)
					outf3 << prs.proms[pindex].genes_sharingproms[s] << ',';
				outf3 << prs.proms[pindex].genes_sharingproms[s];
			}
			else
				outf3 << prs.proms[pindex].sharedpromoter;	
			for (int e = 0; e < NumberofExperiments; ++e)
				outf3 << '\t' << it->normsignal[e];
			outf3 << endl;
		}
	}
	outf3.close();
//Negative control - Negative control interactions
	string FileName4;
	FileName4.append(BaseFileName);
	FileName4.append("SignificantInteractions_NegCtrlstoNegCtrls");
	FileName4.append(".txt");
	ofstream outf4(FileName4.c_str());

	outf4 << "Negctrl1 UCSC" << '\t' << "Number of RE sites within negative ctrl (3000 bases)" << '\t' << "Mappability of negative controls (3000 bases)" << '\t'
		<< "Negctrl2 UCSC" << '\t' << "Number of RE sites within negative ctrl (3000 bases)" << '\t' << "Mappability of negative controls (3000 bases)" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf4 << ExperimentNames[e] << "_Supp_Pairs" << '\t';
	outf4 << endl;

	for(int i = 0; i < ngs.NofNegCtrls; ++i){
		vector< FeattoFeatSignalStruct >::const_iterator it; //first: REpos, second: signal 
		for (it = ngs.negctrls[i].NegcNegcSignals.begin(); it != ngs.negctrls[i].NegcNegcSignals.end(); ++it){
//			double normsignal = (double(it->normsignal)) / prs.proms[i].nofRESites;
//			if (it->normsignal > MinNumberofReads){
			int nindex = it->feat_index;
			outf4 << ngs.negctrls[i].chr << ":" << ngs.negctrls[i].start << "-" << ngs.negctrls[i].end << '\t' << ngs.negctrls[i].nofRESites << '\t' << ngs.negctrls[i].featmappability << '\t';
			outf4 << ngs.negctrls[nindex].chr << ":" << ngs.negctrls[nindex].start << "-" << ngs.negctrls[nindex].end << '\t' << ngs.negctrls[nindex].nofRESites << '\t' << ngs.negctrls[nindex].featmappability << '\t';
			for (int e = 0; e < NumberofExperiments; ++e)
				outf4 << it->normsignal[e] << '\t';
			outf4 << endl;
		}		
	}
outf4.close();
}

void DetectEnrichedBins::PrintNearestPromoter(PromoterClass& prs, int NumberofExperiments){

}