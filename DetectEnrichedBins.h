class DetectEnrichedBins{
public:
	void DetectEnrichedInteractionBins(PromoterClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins(PromoterClass&,string,RESitesClass&);
	void DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass&,string,RESitesClass&);
	void PrintAllInteractions(PromoterClass&,NegCtrlClass&,string, int,vector < string >);
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
}
}

void DetectEnrichedBins::PrintAllInteractions(PromoterClass& prs, NegCtrlClass& ngs, string BaseFileName, int NumberofExperiments,vector < string > ExperimentNames){
//Negative Control Interactions
	string FileName2;
	FileName2.append(BaseFileName);
	FileName2.append("SignificantInteractions");
	FileName2.append("_Enhancers");
	FileName2.append(".txt");
	ofstream outf2(FileName2.c_str());

	outf2 << "Enhancer ID" << '\t' << "Enhancer chr" << '\t' << "Enhancer pos" << '\t' << "Enhancer RE site1" << '\t' << "Enhancer RE site2" << '\t' << "Number of RE sites within core feature" << '\t' 
		  << "Interactor Chr" << '\t' << "Read Position" << '\t'  << "Interactor RE site1" << '\t' << "Interactor RE site2" << '\t' <<  "Distance" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf2 << ExperimentNames[e] << "_SuppPairs" << '\t';
	outf2 << "Peak Profile" << endl;

	for(int i = 0; i < ngs.NofNegCtrls; ++i){
		for(int j = 1; j < ngs.negctrls[i].interactions.size(); ++j){
			outf2 << i << '\t' << ngs.negctrls[i].chr << '\t' << ngs.negctrls[i].midpoint << '\t' << ngs.negctrls[i].closestREsitenums[0] << '\t'
				  << ngs.negctrls[i].closestREsitenums[1] << '\t' << ngs.negctrls[i].nofRESites << '\t'; // Feat Details
			outf2 << ngs.negctrls[i].interactions[j].chr << '\t' << ngs.negctrls[i].interactions[j].pos << '\t' << ngs.negctrls[i].interactions[j].resites[0] << '\t' << ngs.negctrls[i].interactions[j].resites[1] << '\t';
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
//Negative control - Negative control interactions
	string FileName4;
	FileName4.append(BaseFileName);
	FileName4.append("SignificantInteractions_EnhtoEnh");
	FileName4.append(".txt");
	ofstream outf4(FileName4.c_str());

	outf4 << "Enhancer 1 chr" << '\t' << "Enhancer 1 pos" << '\t' << "Enhancer 1 resite1" << '\t' <<  "Enh1 resite2" << '\t' << "Number of RE sites within Enh1 (3000 bases)" << '\t' << "Mappability of Enh1 (3000 bases)" << '\t'
		  << "Enhancer 2 chr" << '\t' << "Enhancer 2 pos" << '\t' << "Enhancer 2 resite1" << '\t' <<  "Enh2 resite2" << '\t' << "Number of RE sites within Enh2 (3000 bases)" << '\t' << "Mappability of Enh2 (3000 bases)" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf4 << ExperimentNames[e] << "_Supp_Pairs" << '\t';
	outf4 << endl;

	for(int i = 0; i < ngs.NofNegCtrls; ++i){
		vector< FeattoFeatSignalStruct >::const_iterator it; //first: REpos, second: signal 
		for (it = ngs.negctrls[i].NegcNegcSignals.begin(); it != ngs.negctrls[i].NegcNegcSignals.end(); ++it){
			int nindex = it->feat_index;
			outf4 << ngs.negctrls[i].chr << '\t' << ngs.negctrls[i].midpoint << '\t' << ngs.negctrls[i].closestREsitenums[0] << '\t' << ngs.negctrls[i].closestREsitenums[1] << '\t' << ngs.negctrls[i].nofRESites << '\t' << ngs.negctrls[i].featmappability << '\t';
			outf4 << ngs.negctrls[nindex].chr << '\t' << ngs.negctrls[nindex].midpoint << '\t' << ngs.negctrls[nindex].closestREsitenums[0] << '\t' << ngs.negctrls[nindex].closestREsitenums[1] << '\t' << ngs.negctrls[nindex].nofRESites << '\t' << ngs.negctrls[nindex].featmappability << '\t';
			for (int e = 0; e < NumberofExperiments; ++e)
				outf4 << it->normsignal[e] << '\t';
			outf4 << endl;
		}		
	}
outf4.close();
}

