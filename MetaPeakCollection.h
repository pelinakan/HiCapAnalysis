vector <MetaPeakMap> metaPeaks;
vector < string > AbNames;
vector <string> ChrNames;

void JoinPeaks(vector < PeakClass >& ap, boost::unordered::unordered_set< int >& p1, int currentpeak, string chr ){

	boost::unordered::unordered_set< int >::const_iterator iter;
	int chrindex;
	string ifpeak; // will keep if the peak of interest present in other files
	for(iter = p1.begin(); iter != p1.end(); ++iter){
		for(int z = 0; z < currentpeak; ++z)
			ifpeak.append("0 ");
		ifpeak.append("1 ");
		int key = (*iter);
		for(int j = currentpeak + 1; j < ap.size();++j){ // starting from the current peak position
			boost::unordered::unordered_map< string, int >::iterator iter1 = ap[j].chr_indexes.find(chr);
			if (ap[j].chr_indexes.find(chr) == ap[j].chr_indexes.end()){
				ifpeak.append("0 ");
			}
			else{
				chrindex = iter1->second; // Get the right index vector
				boost::unordered::unordered_set< int >::iterator it = ap[j].p[chrindex].peaks.find(key);
				if (ap[j].p[chrindex].peaks.find((key)) == ap[j].p[chrindex].peaks.end()) // if the peak is not found
					ifpeak.append("0 ");
				else{
					ifpeak.append("1 ");	
					ap[j].p[chrindex].peaks.erase(it); // Remove the repeated peak
				}
			}
		}
		metaPeaks.back().metapeaks[(key)] = ifpeak; // Add the key and peak profile
		ifpeak.clear();
	} // went through all peaks
}



void FillMetaPeakMap(vector<PeakClass>& AllPeaks){

string chr;
cout << "Chr Names   " << ChrNames[0] << "  " << ChrNames.back() << "  " << ChrNames.size() << endl;
for( int j = 0; j< ChrNames.size();++j){ // For every chromosome 
	chr = ChrNames[j];
	metaPeaks.push_back(MetaPeakMap());
	MetaPeakChrMap[chr] = metaPeaks.size(); 
	for (int k = 0; k < AllPeaks.size(); ++k){
		boost::unordered::unordered_map< string, int >::iterator it = AllPeaks[k].chr_indexes.find(chr);
		if(AllPeaks[k].chr_indexes.find(chr) != AllPeaks[k].chr_indexes.end()){
			int l = it->second; // Get the right index vector
			JoinPeaks(AllPeaks, AllPeaks[k].p[l].peaks,k,chr);
		}
	}
	int nofpeaks = 0;
	for (int n = 0; n < metaPeaks.size();++n){
		nofpeaks += metaPeaks[n].metapeaks.size();
	}
	cout << "Number of Peaks in the Meta Peak Struct " << nofpeaks << endl;
}
cout << "Meta Peak Struct Generated " << endl;

boost::unordered::unordered_map< int, string >::const_iterator iter;
ofstream outf3("MM9.MetaPeaks_All_13_PeakFiles.txt");

outf3 << "chrname" << '\t' << "ClosestRESite" << '\t';
for(int m = 0; m < AbNames.size();++m)
outf3 << AbNames[m] << '\t';
outf3 << endl;
for (int i = 0; i < metaPeaks.size();++i ){
	for(iter = metaPeaks[i].metapeaks.begin(); iter != metaPeaks[i].metapeaks.end(); ++iter){
		outf3 <<  '\t' << iter->first << '\t' << iter->second << endl;
	}
}
}


void ReadMetaPeakFile(void){

	ifstream file("MM9.MetaPeaks_All_13_PeakFiles.txt");
	string temp, chr, peak;
	int repos, index;
	
	cout << "Reading Meta Peak Profile..." << endl;
	file >> temp >> temp;

	for (int i = 0; i < 13; ++i){
		file >> temp;
		AbNames.push_back(temp);
	}
	do{
		file >> chr >> repos;
		for (int i = 0; i < 13; ++i){
			file >> temp;
			peak.append(temp);
		}
		if(MetaPeakChrMap.find(chr) == MetaPeakChrMap.end()){
			index = MetaPeakChrMap.size();
			MetaPeakChrMap[chr] = index;
			metaPeaks.push_back(MetaPeakMap());
			metaPeaks[index].metapeaks[repos] = peak;
		}	
		else{
			boost::unordered::unordered_map< string, int >::const_iterator it = MetaPeakChrMap.find(chr);
			metaPeaks[it->second].metapeaks[repos] = peak;
		}
		peak.clear();
	}while(!file.eof());

cout << "Meta Peak Profile Read" << endl;


}

////////////////////////////////////////////////////////
//READ AND PROCESS PEAKS  ,should be moved to the main body
	/*
ifstream AbNameFile;
AbNameFile.open("AbNames.txt");
vector < PeakClass > AllPeaks;

string fname;
do{
	int abindex = 0;
	string PeakFileName, ext1;

	AbNameFile >> fname >> ext1;
	cout << fname << "  " << ext1 << endl;
	if(fname.compare("END") == 0)
		break;

#ifdef UNIX
	PeakFileName.append(dirname);
#endif
#ifdef WINDOWS
	PeakFileName.append(wdirname);
#endif

	ifstream PeakFile;
	PeakFileName.append(fname);
	AbNames.push_back(ext1);
	PeakClass Peaks;
	PeakFile.open(PeakFileName.c_str());

	if(PeakFile.is_open())
		Peaks.ReadPeakFile(PeakFile,dpnIIsites,ChrNames);
	else
		cerr << PeakFileName << "   cannot be opened"  << endl;

	PeakFile.close();
	
	AllPeaks.push_back(Peaks);
	AllPeaks[abindex].abnames.push_back(ext1);
	cout << AbNames[abindex] << "   read " << endl;
	++abindex;
	PeakFileName.clear();
	fname.clear();
	ext1.clear();
}while(fname != "END");
cout << INTERACTIONFILENAMEBASE << "          All Peak Files Read" << endl;
*/
//FillMetaPeakMap(AllPeaks); // generates the meta peak file, this is generated once then later the program reads the 
//metafile using the function "ReadMetaPeakFile()" above.
//PEAKS ARE PROCESSED METAPEAK MAP GENERATED
/////////////////////////////////////////////////////