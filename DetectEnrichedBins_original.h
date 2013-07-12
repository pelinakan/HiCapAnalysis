class DetectEnrichedBins{
public:
	void DetectEnrichedInteractionBins(PromoterClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins(PromoterClass&,string, RESitesClass&);
	void DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass&,string,RESitesClass&);
	void SummariseInteractions(PromoterClass&, NegCtrlClass&, string, int,ofstream&); //Promoters, basefilename, cell type
	void PrintAllInteractions(PromoterClass&,NegCtrlClass&,string, int,vector < string >);
};

void DetectEnrichedBins::DetectEnrichedInteractionBins(PromoterClass &prs, DetermineBackgroundLevels background, string BaseFileName,int ExperimentNo){
/*
	string FileName;
	FileName.append(BaseFileName);
	FileName.append("InteractionsPerPromoter.txt");
	ofstream outf1(FileName.c_str());
*/
	MappabilityClass mapp;
	mapp.InitialiseVars();
	
//	ofstream tempf("Mappabilities.txt");
	// Upstream
	for(int i = 0; i < prs.NofPromoters; ++i){
//		outf1 << prs.proms[i].genes[0] << '\t' << prs.proms[i].chr << ":" << prs.proms[i].start << "-" << prs.proms[i].end << '\t' << prs.proms[i].strand << '\t' << prs.proms[i].nofRESites << '\t';
		boost::unordered::unordered_map<int, int* >::const_iterator it; //first: REpos, second: signal 
		for (it = prs.proms[i].Signals.signal_ups.begin(); it != prs.proms[i].Signals.signal_ups.end(); ++it){
//			double normsignal = (double(it->second)) / prs.proms[i].nofRESites;
			double normsignal = (double(it->second[ExperimentNo]));
			if (it->second[0] > MinNumberofReads){
				double m = mapp.GetMappability(prs.proms[i].chr, (it->first - 100), 200);
	//			tempf << prs.proms[i].chr << ":" << (it->first - 100) << "-" << (it->first + 100) << '\t' << m << endl;
				if (m > Mappability_Threshold){
					double q, p_value;
					int dist;
					if (prs.proms[i].strand =="+")
						dist = it->first - prs.proms[i].closestREsitenums[0]; // negative
					else
						dist = prs.proms[i].closestREsitenums[1] - it->first; // negative
					int bin = abs(dist) / BinSize; 
//					if((background.bglevels.mean_upstream).find(bin) == background.bglevels.mean_upstream.end()){
//						background.bglevels.mean_upstream[bin] = 0;
//						background.bglevels.stdev_upstream[bin] = 1;
//						p_value = 0.0;
//					}
//					else{
						q = alglib::normaldistribution(( normsignal - background.bglevels.mean_upstream[bin]) / background.bglevels.stdev_upstream[bin]); //z-score
						p_value = 1- q;
//					}
//					cout << i << '\t' << dist << '\t' << bin << '\t' << normsignal << '\t' << background.bglevels.mean_downstream[bin] << '\t' << background.bglevels.stdev_upstream[bin] << '\t' << q << '\t' << p_value << endl;
					if (p_value < SignificanceThreshold){ 
						prs.proms[i].interactions.push_back(InteractionStruct());
						prs.proms[i].interactions.back().peakprofile.clear();
						prs.proms[i].interactions.back().p_val[ExperimentNo] = p_value;
						prs.proms[i].interactions.back().mappability = m;
						prs.proms[i].interactions.back().supp_pairs[ExperimentNo] = it->second[ExperimentNo];
						prs.proms[i].interactions.back().chr = prs.proms[i].chr;
						prs.proms[i].interactions.back().distance = dist; 
//						boost::unordered::unordered_map<int, float >::const_iterator itd = prs.proms[i].Signals.distance.find(it->first);
//						prs.proms[i].interactions.back().distance = itd->second / prs.proms[i].interactions.back().norm_signal;
						prs.proms[i].interactions.back().pos[0] = it->first;
						prs.proms[i].interactions.back().type = 'U'; // U: upstream , D: downstream , X: inter-chromosomal
//						outf1 << prs.proms[i].interactions.back().distance << ", " << prs.proms[i].interactions.back().norm_signal << '\t';
					}
				}
			}
		}
	//Downstream
		for (it = prs.proms[i].Signals.signal_down.begin(); it != prs.proms[i].Signals.signal_down.end(); ++it){
		//	double normsignal = (double(it->second)) / prs.proms[i].nofRESites;
			double normsignal = (double(it->second[ExperimentNo]));		
			if( normsignal > MinNumberofReads){ //signal is greater than the min number of reads
				double m = mapp.GetMappability(prs.proms[i].chr, (it->first - 100), 200);
				if (m > Mappability_Threshold){				
					double q, p_value;
					int dist;
					if (prs.proms[i].strand =="+")
						dist = it->first - prs.proms[i].closestREsitenums[1]; // positive
					else
						dist = prs.proms[i].closestREsitenums[0] - it->first; // positive
					int bin = abs(dist) / BinSize; 
	//				if((background.bglevels.mean_downstream).find(bin) == background.bglevels.mean_downstream.end()){
	//					background.bglevels.mean_downstream[bin] = 0;
	//					background.bglevels.stdev_downstream[bin] = 1;
	//					p_value = 0.0;
	//				}
	//				else{
					q = alglib::normaldistribution(( normsignal - background.bglevels.mean_downstream[bin]) / background.bglevels.stdev_downstream[bin]); //z-score
					p_value = 1- q;
	//				}
//					cout << i << '\t' << dist << '\t' << bin << '\t' << normsignal << '\t' << background.bglevels.mean_downstream[bin] << '\t' << background.bglevels.stdev_upstream[bin] << '\t' << q << '\t' << p_value << endl;
					if (p_value < SignificanceThreshold){
						prs.proms[i].interactions.push_back(InteractionStruct());
						prs.proms[i].interactions.back().peakprofile.clear();
						prs.proms[i].interactions.back().p_val[ExperimentNo] = p_value;
						prs.proms[i].interactions.back().mappability = m;
						prs.proms[i].interactions.back().supp_pairs[ExperimentNo] = it->second[ExperimentNo];
						prs.proms[i].interactions.back().chr = prs.proms[i].chr;
						prs.proms[i].interactions.back().distance = dist; 
//						boost::unordered::unordered_map<int, float >::const_iterator itd = prs.proms[i].Signals.distance.find(it->first);
//						prs.proms[i].interactions.back().distance = -1 * (itd->second / prs.proms[i].interactions.back().norm_signal);
						prs.proms[i].interactions.back().pos[0] = it->first;
						prs.proms[i].interactions.back().type = 'D'; // U: upstream , D: downstream , X: inter-chromosomal
//						outf1 << prs.proms[i].interactions.back().distance << ", " << prs.proms[i].interactions.back().norm_signal << '\t';
					}
				}
			}
		}
	// Inter-chromosomal
		for( int j = 0; j < prs.proms[i].Signals_CTX.size(); ++j){ // For each chromosome
			for (it = prs.proms[i].Signals_CTX[j].signal.begin(); it != prs.proms[i].Signals_CTX[j].signal.end(); ++it){
			//	double normsignal = (double(it->second)) / prs.proms[i].nofRESites;
				double normsignal = (double(it->second[ExperimentNo]));
				if( normsignal > MinNumberofReads){ //signal is greater than twice the min number of reads required
					double m = mapp.GetMappability(prs.proms[i].Signals_CTX[j].maptochrname, (it->first - 100), 200);
					if (m > Mappability_Threshold){
						prs.proms[i].interactions.push_back(InteractionStruct());
						prs.proms[i].interactions.back().peakprofile.clear();
						prs.proms[i].interactions.back().p_val[ExperimentNo] = -1;
						prs.proms[i].interactions.back().mappability = m;
						prs.proms[i].interactions.back().chr.append(prs.proms[i].Signals_CTX[j].maptochrname);
						prs.proms[i].interactions.back().distance = -1;
						prs.proms[i].interactions.back().supp_pairs[ExperimentNo] = it->second[ExperimentNo];
						prs.proms[i].interactions.back().pos[0] = it->first;
						prs.proms[i].interactions.back().type = 'X'; // U: upstream , D: downstream , X: inter-chromosomal
//						outf1 << prs.proms[i].interactions.back().chr << ":" << prs.proms[i].interactions.back().pos << ", " << prs.proms[i].interactions.back().norm_signal << '\t';
					}
				}
			}
		}
//		outf1 << endl;
	}
}
void DetectEnrichedBins::AssociatePeaksWithIntBins(PromoterClass& prs,string BaseFileName, RESitesClass& dpnII){

	for(int i = 0; i < prs.NofPromoters; ++i){
		for(int j = 0; j < prs.proms[i].interactions.size(); ++j){
			bool peakfound = 0;
			int *renums;
			renums = new int [2];
			if(prs.proms[i].interactions[j].type == 'D'){
				dpnII.GettheREPositions(prs.proms[i].chr,prs.proms[i].interactions[j].pos[0],renums); // Interactor RE fragment
				prs.proms[i].interactions[j].pos[1] = renums[1];				
			}
			else{
				dpnII.GettheREPositions(prs.proms[i].chr,(prs.proms[i].interactions[j].pos[0] - 1),renums); // Interactor RE fragment
				prs.proms[i].interactions[j].pos[1] = renums[0];				
			}
			boost::unordered::unordered_map<string, int>::const_iterator citer = MetaPeakChrMap.find(prs.proms[i].interactions[j].chr);
			if (MetaPeakChrMap.find(prs.proms[i].interactions[j].chr) == MetaPeakChrMap.end()) // Find the right peak map wrt chromosome
				break;
			if(metaPeaks[citer->second].metapeaks.find(prs.proms[i].interactions[j].pos[0]) != metaPeaks[citer->second].metapeaks.end()){
				//If there is peak
				boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(prs.proms[i].interactions[j].pos[0]);
				prs.proms[i].interactions[j].peakoverlap = true;
				prs.proms[i].interactions[j].peakprofile = it->second;
				peakfound = 1;
			}
			else{
				int *repos;
				repos = new int [2];
				dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(prs.proms[i].interactions[j].pos[0] - 1),repos);
				do{
					if(metaPeaks[citer->second].metapeaks.find(repos[0]) != metaPeaks[citer->second].metapeaks.end()){
						//If there is peak
						boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[0]);
						prs.proms[i].interactions[j].peakoverlap = true;
						prs.proms[i].interactions[j].peakprofile = it->second;
						peakfound = 1;
						break;
					}
					else
						dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(repos[0] - 1),repos);
				}while((abs(repos[0] - prs.proms[i].interactions[j].pos[0]) < 1000));
				if(!peakfound){
					dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(prs.proms[i].interactions[j].pos[0]),repos);
					do{
						if(metaPeaks[citer->second].metapeaks.find(repos[1]) != metaPeaks[citer->second].metapeaks.end()){
							//If there is peak
							boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[1]);
							prs.proms[i].interactions[j].peakoverlap = true;
							prs.proms[i].interactions[j].peakprofile = it->second;
							peakfound = 1;
							break;
						}
						else
							dpnII.GettheREPositions(prs.proms[i].interactions[j].chr,(repos[1]),repos);
					}while ( abs(repos[1] - prs.proms[i].interactions[j].pos[0]) < 1000 );
				}
				if (!peakfound)
					prs.proms[i].interactions[j].peakoverlap = false;
			}
		}
	}
}
void DetectEnrichedBins::DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass &ngs,DetermineBackgroundLevels background,string BaseFileName,int ExperimentNo){
/*
	string FileName;
	FileName.append(BaseFileName);
	FileName.append("InteractionsPerNegCtrl");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());
*/
// Upstream
	for(int i = 0; i < ngs.NofNegCtrls; ++i){
//		outf1 << ngs.negctrls[i].chr << ":" << ngs.negctrls[i].start << "-" << ngs.negctrls[i].end << '\t';
		boost::unordered::unordered_map<int, int* >::const_iterator it;
		for (it = ngs.negctrls[i].Signals.signal_ups.begin(); it != ngs.negctrls[i].Signals.signal_ups.end(); ++it){
			//double normsignal = (double(it->second)) / ngs.negctrls[i].nofRESites;
			double normsignal = (double(it->second[ExperimentNo]));
			if( normsignal > MinNumberofReads){ //signal is greater than the min number of reads
				double q, p_value;
				double dist = ngs.negctrls[i].closestREsitenums[0] - it->first;
				int bin = dist / BinSize; 
				if((background.bglevels.mean_upstream).find(bin) == background.bglevels.mean_upstream.end()){
					background.bglevels.mean_upstream[bin] = 0;
					background.bglevels.stdev_upstream[bin] = 1;
					p_value = 0.0;
				}
				else{
					q = alglib::normaldistribution(( normsignal - background.bglevels.mean_upstream[bin])/background.bglevels.stdev_upstream[bin]); //z-score
					p_value=1-q;
				}
				if (p_value < SignificanceThreshold){
					ngs.negctrls[i].interactions.push_back(InteractionStruct());
					ngs.negctrls[i].interactions.back().peakprofile.clear();
					ngs.negctrls[i].interactions.back().chr = ngs.negctrls[i].chr;
					ngs.negctrls[i].interactions.back().distance = it->first - ngs.negctrls[i].closestREsitenums[0];
					ngs.negctrls[i].interactions.back().supp_pairs[ExperimentNo] = it->second[ExperimentNo];
					ngs.negctrls[i].interactions.back().type = 'U';
					ngs.negctrls[i].interactions.back().pos[0] = it->first;
//					outf1 << ngs.negctrls[i].interactions.back().distance << ", " << ngs.negctrls[i].interactions.back().norm_signal << '\t';
				}
			}
		}
//Downstream
		for (it = ngs.negctrls[i].Signals.signal_down.begin(); it != ngs.negctrls[i].Signals.signal_down.end(); ++it){
			//double normsignal = (double(it->second)) / ngs.negctrls[i].nofRESites;
			double normsignal = (double(it->second[ExperimentNo]));
			if( normsignal > MinNumberofReads){ //signal is greater than the min number of reads
				double q, p_value;
				double dist = ngs.negctrls[i].closestREsitenums[0] - it->first;
				int bin = dist / BinSize; 
				if((background.bglevels.mean_downstream).find(bin) == background.bglevels.mean_downstream.end()){
					background.bglevels.mean_downstream[bin] = 0;
					background.bglevels.stdev_downstream[bin] = 1;
					p_value = 0.0;
				}
				else{
					q = alglib::normaldistribution(( normsignal - background.bglevels.mean_downstream[bin])/background.bglevels.stdev_downstream[bin]); //z-score
					p_value=1 - q;
				}
				if (p_value < SignificanceThreshold){
					ngs.negctrls[i].interactions.push_back(InteractionStruct());
					ngs.negctrls[i].interactions.back().peakprofile.clear();
					ngs.negctrls[i].interactions.back().chr = ngs.negctrls[i].chr;
					ngs.negctrls[i].interactions.back().distance = dist;
					ngs.negctrls[i].interactions.back().supp_pairs[ExperimentNo] = it->second[ExperimentNo];
					ngs.negctrls[i].interactions.back().type = 'D';
					ngs.negctrls[i].interactions.back().pos[0] = it->first;
//					outf1 << ngs.negctrls[i].interactions.back().distance << ", " << ngs.negctrls[i].interactions.back().norm_signal << '\t';
				}
			}
		}
// Inter-chromosomal
		for( int j = 0; j < ngs.negctrls[i].Signals_CTX.size(); ++j){ // For each chromosome
			for (it = ngs.negctrls[i].Signals_CTX[j].signal.begin(); it != ngs.negctrls[i].Signals_CTX[j].signal.end(); ++it){
			//	double normsignal = it->second / ngs.negctrls[i].nofRESites;
				double normsignal = (double(it->second[ExperimentNo]));
				//	if( normsignal > MinNumberofReads){ //signal is greater than the min number of reads
				if (it->second[ExperimentNo] > MinNumberofReads ){
					ngs.negctrls[i].interactions.push_back(InteractionStruct());
					ngs.negctrls[i].interactions.back().peakprofile.clear();
					ngs.negctrls[i].interactions.back().chr = ngs.negctrls[i].Signals_CTX[j].maptochrname;
					ngs.negctrls[i].interactions.back().distance = -1;
					ngs.negctrls[i].interactions.back().supp_pairs[ExperimentNo] = it->second[ExperimentNo];
					ngs.negctrls[i].interactions.back().type = 'X';
					ngs.negctrls[i].interactions.back().pos[0] = it->first;
//					outf1 << ngs.negctrls[i].interactions.back().chr << ":" << ngs.negctrls[i].interactions.back().pos << ", " << ngs.negctrls[i].interactions.back().norm_signal << '\t';
				}
			}
		}
//		outf1 << endl;
	}
}
void DetectEnrichedBins::AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass& ngs,string BaseFileName,RESitesClass& dpnII){

	for(int i = 0; i < ngs.NofNegCtrls; ++i){
		for(int j = 1; j < ngs.negctrls[i].interactions.size(); ++j){
			int *renums;
			renums = new int [2];
			dpnII.GettheREPositions(ngs.negctrls[i].chr,ngs.negctrls[i].interactions[j].pos[0],renums); // Interactor RE fragment
			ngs.negctrls[i].interactions[j].pos[1] = renums[1];				

			bool peakfound = 0;
			boost::unordered::unordered_map<string, int>::const_iterator citer = MetaPeakChrMap.find(ngs.negctrls[i].chr);
			if (MetaPeakChrMap.find(ngs.negctrls[i].chr) == MetaPeakChrMap.end())
				break;
			if(metaPeaks[citer->second].metapeaks.find(ngs.negctrls[i].interactions[j].pos[0]) != metaPeaks[citer->second].metapeaks.end()){
				//If there is peak
				boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(ngs.negctrls[i].interactions[j].pos[0]);
				ngs.negctrls[i].interactions[j].peakoverlap = true;
				ngs.negctrls[i].interactions[j].peakprofile = it->second;
				peakfound = 1;
			}
			else{
				int *repos;
				repos = new int [2];
				dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,(ngs.negctrls[i].interactions[j].pos[0] - 1),repos);
				do{
					if(metaPeaks[citer->second].metapeaks.find(repos[0]) != metaPeaks[citer->second].metapeaks.end()){
						//If there is peak
						boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[0]);
						ngs.negctrls[i].interactions[j].peakoverlap = true;
						ngs.negctrls[i].interactions[j].peakprofile = it->second;
						peakfound = 1;
						break;
					}
					else
						dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,(repos[0] - 1),repos);
				}while((abs(repos[0] - ngs.negctrls[i].interactions[j].pos[0]) < 1000));
				if(!peakfound){
					dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,ngs.negctrls[i].interactions[j].pos[0],repos);
					do{
						if(metaPeaks[citer->second].metapeaks.find(repos[1]) != metaPeaks[citer->second].metapeaks.end()){
							//If there is peak
							boost::unordered::unordered_map<int, string >::const_iterator it = metaPeaks[citer->second].metapeaks.find(repos[1]);
							ngs.negctrls[i].interactions[j].peakoverlap = true;
							ngs.negctrls[i].interactions[j].peakprofile = it->second;
							peakfound = 1;
							break;
						}
						else
							dpnII.GettheREPositions(ngs.negctrls[i].interactions[j].chr,(repos[1]),repos);
					}while ( abs(repos[1] - ngs.negctrls[i].interactions[j].pos[0]) < 1000 );
				}
				if (!peakfound)
					ngs.negctrls[i].interactions[j].peakoverlap = false;
			}			
		}
	}
}

void DetectEnrichedBins::SummariseInteractions(PromoterClass& prs, NegCtrlClass& ngs, string BaseFileName, int CellType, ofstream& SummaryFile){

	string FileName;
	FileName.append(BaseFileName);
	FileName.append("Promoters_Summary");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	double T = 0, TwP = 0, iT = 0, xT = 0, iTwP = 0, xTwP = 0;

	outf1 << "Promoter Name" << '\t' << "Representative Transcript Name" << '\t' << "Promoter chr"  << '\t' << "Promoter TSS"  << '\t' << "Promoter UCSC" << '\t' 
		  << "Gene Expression" << '\t' << "Number of RE sites within core promoter" <<  '\t' << "Mappability of core promoter" << '\t' << "If it is a shared promoter" << '\t'
		  << "Number of upstream interactions overlapping at least one peak" << '\t' << "Number of upstream interactions not overlapping any peak" << '\t' 
		  << "Total Number of upstream interactions" << '\t'
		  << "Number of downstream interactions overlapping at least one peak" << '\t' << "Number of downstream interactions not overlapping any peak" << '\t' 
		  << "Total Number of downstream interactions" << '\t'
		  << "Number of intra-chr interactions overlapping at least one peak" << '\t' << "Number of intra-chr interactions not overlapping any peak" << '\t' 
		  << "Total Number of intra-chr interactions" << '\t'
		  << "Number of inter-chr interactions overlapping at least one peak" << '\t' << "Number of inter-chr interactions not overlapping any peak" << '\t' 
		  << "Total Number of inter-chr interactions" << '\t' 
          << "Number of prom-prom interactions" <<  endl;

	for (int i = 0 ; i < prs.NofPromoters; ++i){
		outf1 << prs.proms[i].genes[0] << '\t' << prs.proms[i].transcripts[0] << '\t' << prs.proms[i].chr << '\t' << prs.proms[i].TSS << '\t' << prs.proms[i].chr << ":" << prs.proms[i].start << "-" << prs.proms[i].end << '\t' 
			  << prs.proms[i].expression[CellType] << '\t' << prs.proms[i].nofRESites << '\t' << prs.proms[i].featmappability << '\t' << prs.proms[i].sharedpromoter << '\t';
		double iwp = 0.0, iwop = 0.0; // intra-chr significant interactions with or without peaks
		double iwpd = 0.0, iwopd = 0.0; // intra-chr downstream significant interactions with or without peaks
		double iwpu = 0.0, iwopu = 0.0; // intra-chr downstream significant interactions with or without peaks
		double iwpx = 0.0, iwopx = 0.0; // intra-chr significant interactions with or without peaks
		for(int j = 0; j < prs.proms[i].interactions.size();++j){
			++T;
			if(prs.proms[i].interactions[j].type == 'D'){
				++iT;
				if(prs.proms[i].interactions[j].peakprofile.empty()){
					++iwop;
					++iwopd;
				}
				else{
					++iwp;
					++iwpd;
					++TwP;
					++iTwP;
				}
			}
			if(prs.proms[i].interactions[j].type == 'U'){
				++iT;
				if(prs.proms[i].interactions[j].peakprofile.empty()){
					++iwop;
					++iwopu;
				}
				else{
					++iwp;
					++iwpu;
					++TwP;
					++iTwP;
				}
			} 
			if(prs.proms[i].interactions[j].type == 'X'){
				++xT;
				if(prs.proms[i].interactions[j].peakprofile.empty())
					++iwopx;
				else{
					++iwpx;
					++TwP;
					++xTwP;
				}
			}
		}
		outf1 << iwpu << '\t' << iwopu << '\t' <<  (iwpu + iwopu) << '\t'
			  << iwpd << '\t' << iwopd << '\t' <<  (iwpd + iwopd) << '\t'
			  << iwp  << '\t' << iwop  << '\t' <<  (iwp + iwop)   << '\t'
			  << iwpx << '\t' << iwopx << '\t' <<  (iwpx + iwopx) << '\t';

		outf1 << prs.proms[i].promPromSignals.size() << endl;
	}
	outf1.close();


	SummaryFile << "Total number of intra-chromosomal interactions of promoters detected (prom-prom interactions excluded)" << '\t' << iT << endl
				<< "Total number of intra-chromosomal interactions overlapping at least one peak" << '\t' << iTwP << endl 
				<< "Total number of intra-chromosomal interactions not overlapping any peak" << '\t' << (iT-iTwP) << endl
				<< "Intra-chromosomal interactions peak overlap ratio" << '\t' << (iTwP/iT) << endl << endl
				<< "Total number of inter-chromosomal interactions of promoters detected (prom-prom interactions excluded)" << '\t' << xT << endl
				<< "Total number of inter-chromosomal interactions overlapping at least one peak" << '\t' << xTwP << endl 
				<< "Total number of inter-chromosomal interactions not overlapping any peak" << '\t' << (xT-xTwP) << endl
				<< "Inter-chromosomal interactions peak overlap ratio" << '\t' << (xTwP/xT) << endl << endl	
				<< "Total number of all distal interactions of promoters detected (prom-prom interactions excluded)" << '\t' << T << endl
				<< "Total number of all distal interactions overlapping at least one peak" << '\t' << TwP << endl 
				<< "Total number of all distal interactions not overlapping any peak" << '\t' << (T-TwP) << endl
				<< "All distal interactions peak overlap ratio" << '\t' << (TwP/T) << endl << endl;		


	string FileName2;
	FileName2.append(BaseFileName);
	FileName2.append("NegCtrls_Summary");
	FileName2.append(".txt");
	ofstream outf2(FileName2.c_str());

	outf2 << "NegCtrl ID" << '\t' << "NegCtrl UCSC" << '\t' << "Number of RE sites within core region" << '\t' << "Mappability of the feature (3000 bases)" << '\t' 
	      << "Number of upstream interactions overlapping at least one peak" << '\t' << "Number of upstream interactions not overlapping any peak" << '\t' 
	  	<< "Total Number of upstream interactions" << '\t'
		<< "Number of downstream interactions overlapping at least one peak" << '\t' << "Number of downstream interactions not overlapping any peak" << '\t' 
		<< "Total Number of downstream interactions" << '\t'
		<< "Number of intra-chr interactions overlapping at least one peak" << '\t' << "Number of intra-chr interactions not overlapping any peak" << '\t' 
		<< "Total Number of intra-chr interactions" << '\t'
		<< "Number of inter-chr interactions overlapping at least one peak" << '\t' << "Number of inter-chr interactions not overlapping any peak" << '\t' 
		<< "Total Number of inter-chr interactions" <<  '\t' 
		<< "Number of feat-feat interactions"  << endl;

	T = 0.0;
	TwP = 0.0; 
	iT = 0.0; xT = 0.0;
	iTwP = 0.0; xTwP = 0.0;
	for (int i = 0 ; i < ngs.NofNegCtrls; ++i){
		outf2 << i << '\t' << ngs.negctrls[i].chr << ":" << ngs.negctrls[i].start << "-" << ngs.negctrls[i].end << '\t' << ngs.negctrls[i].nofRESites << '\t' << ngs.negctrls[i].featmappability << '\t';
		double iwp = 0.0, iwop = 0.0; // intra-chr significant interactions with or without peaks
		double iwpd = 0.0, iwopd = 0.0; // intra-chr downstream significant interactions with or without peaks
		double iwpu = 0.0, iwopu = 0.0; // intra-chr downstream significant interactions with or without peaks
		double iwpx = 0.0, iwopx = 0.0; // intra-chr significant interactions with or without peaks
		for(int j = 0; j < ngs.negctrls[i].interactions.size();++j){
			++T;
			if(ngs.negctrls[i].interactions[j].type == 'D'){
				++iT;
				if(ngs.negctrls[i].interactions[j].peakprofile.empty()){
					++iwop;
					++iwopd;
				}
				else{
					++iwp;
					++iwpd;
					++TwP;
					++iTwP;
				}
			} 
			if(ngs.negctrls[i].interactions[j].type == 'U'){
				++iT;
				if(ngs.negctrls[i] .interactions[j].peakprofile.empty()){
					++iwop;
					++iwopu;
				}
				else{
					++iwp;
					++iwpu;
					++TwP;
					++iTwP;
				}
			} 
			if(ngs.negctrls[i].interactions[j].type == 'X'){
				++xT;
				if(ngs.negctrls[i].interactions[j].peakprofile.empty())
					++iwopx;
				else{
					++iwpx;
					++TwP;
					++xTwP;
				}				
			}
		}
		outf2 << iwpu << '\t' << iwopu << '\t' <<  (iwpu + iwopu) << '\t'
			  << iwpd << '\t' << iwopd << '\t' <<  (iwpd + iwopd) << '\t'
			  << iwp  << '\t' << iwop  << '\t' <<  (iwp + iwop)   << '\t'
			  << iwpx << '\t' << iwopx << '\t' <<  (iwpx + iwopx) << '\t';
		outf2 << ngs.negctrls[i].NegcNegcSignals.size() << endl;


	}
	outf2.close();

	SummaryFile << "Total number of intra-chromosomal interactions of negative controls detected (negctrl-negctrl interactions excluded)" << '\t' << iT << endl
		<< "Total number of intra-chromosomal interactions overlapping at least one peak" << '\t' << iTwP << endl 
		<< "Total number of intra-chromosomal interactions not overlapping any peak" << '\t' << (iT-iTwP) << endl
		<< "Intra-chromosomal interactions peak overlap ratio" << '\t' << (iTwP/iT) << endl << endl
	    << "Total number of inter-chromosomal interactions of negative controls detected (negctrl-negctrl interactions excluded)" << '\t' << xT << endl
		<< "Total number of inter-chromosomal interactions overlapping at least one peak" << '\t' << xTwP << endl 
		<< "Total number of inter-chromosomal interactions not overlapping any peak" << '\t' << (xT-xTwP) << endl
		<< "Inter-chromosomal interactions peak overlap ratio" << '\t' << (xTwP/xT) << endl << endl		
	    << "Total number of all distal interactions of negative controls detected (negctrl-negctrl interactions excluded)" << '\t' << T << endl
		<< "Total number of all distal interactions overlapping at least one peak" << '\t' << TwP << endl 
		<< "Total number of all distal interactions not overlapping any peak" << '\t' << (T-TwP) << endl
		<< "All distal interactions peak overlap ratio" << '\t' << (TwP/T) << endl << endl;		

}
void DetectEnrichedBins::PrintAllInteractions(PromoterClass& prs, NegCtrlClass& ngs, string BaseFileName, int NumberofExperiments,vector < string > ExperimentNames){
//Promoter Interactions
	string FileName;
	FileName.append(BaseFileName);
	FileName.append("SignificantInteractions");
	FileName.append("_Promoters");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());
	int s = 0;

	outf1 << "Gene Name" << '\t' << "Representative Transcript Name" << '\t' << "Gene Expression" << '\t' << "Shared Promoter" << '\t' << "Number of RE Sites within core promoter" << '\t' << "Mappability of Core Promoter" << '\t' 
		<< "Promoter chr" << '\t' << "Promoter TSS"  << '\t' << "Strand" << '\t' <<"Promoter UCSC" << '\t' << "Interactor chr" << '\t' << "Interactor REsite1" << '\t' << "Interactor REsite2" << '\t'
		<< "Interactor Mappability (100 bases around the RE site)" << '\t' << "Distance" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf1 << ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_pval" ;
	outf1 << '\t' << "Peak Profile" << endl;

	for(int i = 0; i < prs.NofPromoters; ++i){
		for(int j = 0; j < prs.proms[i].interactions.size(); ++j){
			outf1 << prs.proms[i].genes[0] << '\t' << prs.proms[i].transcripts[0] << '\t' << prs.proms[i].expression[0] << '\t';
			if (prs.proms[i].sharedpromoter){
				for (s = 0; s < prs.proms[i].genes_sharingproms.size() - 1; ++s)
					outf1 << prs.proms[i].genes_sharingproms[s] << ',';
				outf1 << prs.proms[i].genes_sharingproms[s];
			}
			else
				outf1 << prs.proms[i].sharedpromoter;
			outf1 << '\t' << prs.proms[i].nofRESites << '\t' << prs.proms[i].featmappability << '\t' 
				<< prs.proms[i].chr << '\t' << prs.proms[i].TSS << '\t' << prs.proms[i].strand << '\t' << prs.proms[i].chr << ":" << prs.proms[i].start << "-" << prs.proms[i].end << '\t'; // Promoter Details
			outf1 << prs.proms[i].interactions[j].chr << '\t' << prs.proms[i].interactions[j].pos[0] << '\t' << prs.proms[i].interactions[j].pos[1] << '\t';
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
		<< "Distance" << '\t';
	for (int e = 0; e < NumberofExperiments; ++e)
		outf2 << ExperimentNames[e] << "_SuppPairs";
	outf2 << '\t' << "Peak Profile" << endl;

	for(int i = 0; i < ngs.NofNegCtrls; ++i){
		for(int j = 1; j < ngs.negctrls[i].interactions.size(); ++j){
			outf2 << i << '\t' << ngs.negctrls[i].chr << ":" << ngs.negctrls[i].start << "-" << ngs.negctrls[i].end << '\t' << ngs.negctrls[i].nofRESites << '\t'; // Feat Details
			outf2 << ngs.negctrls[i].interactions[j].chr << '\t' << ngs.negctrls[i].interactions[j].pos[0] << '\t' << ngs.negctrls[i].interactions[j].pos[1] << '\t';
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
	ofstream outf4(FileName2.c_str());

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
				outf3 << it->normsignal[e] << '\t';
			outf3 << endl;
		}		
	}
outf4.close();
}

