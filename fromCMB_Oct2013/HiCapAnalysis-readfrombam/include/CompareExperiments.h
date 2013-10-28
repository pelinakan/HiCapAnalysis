struct Int_Struct{
	string promchr;
	int prompos;
	string intchr;
	int intpos;
	bool present[6];
	int paircounts[6];
	string allrow;
};

class CompareExperiments{
public:
	vector < Int_Struct > allinteractions;
	vector < string > ExpName;
	int intindex;
	string header;
	boost::unordered::unordered_map< int, int > intmap; // it->first interaction REpos, it->second : index of the Int_Struct
	int NumberofExperiments;
	void readandCompare(void);
private:
	void getInteractionDetails(stringstream&,Int_Struct&, int&);
	bool compareInteractions(int, Int_Struct&);
	void readFile(string,int);
	void PrintAll(string);
};

void CompareExperiments::readandCompare(void){

ifstream filenames("BaseFileNames.txt");
	string fname;
	int fileindex = 0;
	intindex = 0;
	
	NumberofExperiments = 3;

	filenames >> fname;
do{
	ExpName.push_back(fname);
	readFile(fname,fileindex);
	cout << fname << "   read" << endl;
	++fileindex;
	filenames >> fname;
}while(!(fname.compare("END") == 0));
cout << "All Interactions Files Read" << endl;
PrintAll(header);
}
void CompareExperiments::readFile(string filename, int fileindex){

	Int_Struct tempstruct;
	string s;
	int temppaircount = 0;

	ifstream file(filename.c_str());

	getline(file,header);
	getline(file,s); 
	stringstream ss ( s );
	while(!(s.empty())){
		getInteractionDetails(ss, tempstruct, temppaircount);
		if (intmap.find(tempstruct.intpos) == intmap.end()){
			allinteractions.push_back(tempstruct);
			for ( int i = 0; i< NumberofExperiments; ++i){
				allinteractions.back().present[i] = false;
				allinteractions.back().paircounts[i] = 0;
			}
			allinteractions.back().allrow.append(ss.str()); 
			allinteractions.back().present[fileindex] = true;
			allinteractions.back().paircounts[fileindex] = temppaircount;
			intmap[allinteractions.back().intpos] = (allinteractions.size() - 1);
		}
		else{
			boost::unordered::unordered_map< int, int >::iterator it = intmap.find(tempstruct.intpos);
			if(allinteractions[it->second].present[fileindex] != true){
				allinteractions[it->second].present[fileindex] = compareInteractions(it->second, tempstruct);
				if (allinteractions[it->second].present[fileindex])
					allinteractions[it->second].paircounts[fileindex] = temppaircount;
			}
		}
		ss.str("");
		getline(file,s); 
		ss.str(s);	
		++intindex;
	}
file.close();
}
void CompareExperiments::getInteractionDetails(stringstream& ss,Int_Struct& tempstruct, int& temppaircount){

	string chr, field;

	getline(ss,field,'\t'); // gene name
	getline(ss,field,'\t'); // tr name
	getline(ss,field,'\t'); // exp
	getline(ss,field,'\t'); // sharedproms
	getline(ss,field,'\t'); // nof re sites
	getline(ss,field,'\t'); // mapp
	getline(ss,tempstruct.promchr,'\t'); // 
	getline(ss,field,'\t'); // 
	tempstruct.prompos = (atoi(field.c_str()));
	getline(ss,field,'\t'); // ignore all until interaction type
	getline(ss,field,'\t'); // 
	getline(ss,field,'\t'); // 
	getline(ss,field,'\t'); // type
	getline(ss,tempstruct.intchr,'\t');
	getline(ss,field,'\t');
	tempstruct.intpos = (atoi(field.c_str())); 
	getline(ss,field,'\t');
	getline(ss,field,'\t');
	getline(ss,field,'\t');
	temppaircount = (atoi(field.c_str()));
	//	Pecam1	NM_008816	9.92781	0	7	0.96273	chr11	106576595	-	chr11:106576095-106579095	chr5:3837639-3837839	CTX	chr5	3837739	1	-1	8	NoPeakOverlap

}
bool CompareExperiments::compareInteractions(int pidx, Int_Struct& tempstruct){ // current idx and past idx

	if ( tempstruct.intchr.compare(allinteractions[pidx].intchr) == 0){
		if (tempstruct.promchr.compare(allinteractions[pidx].promchr) == 0){
			if ( tempstruct.prompos == allinteractions[pidx].prompos)
				return 1;
		}
	}
	return 0;
}

void CompareExperiments::PrintAll(string header){

	ofstream ofile("Comparison.txt");
	ofile << header << '\t';
	for (int i = 0; i < NumberofExperiments; ++i)
		ofile << "Present in experiment " << ExpName[i] << '\t';
	for (int i = 0; i < NumberofExperiments; ++i)
		ofile << "Supp. pair count in experiment " << ExpName[i] << '\t';
	ofile << endl;

	vector < Int_Struct >::iterator it;
	for (it = allinteractions.begin(); it < allinteractions.end(); ++it ){
		ofile << it->allrow << '\t';
		for (int i = 0; i < NumberofExperiments; ++i)
			ofile << it->present[i] << '\t';
		for (int i = 0; i < NumberofExperiments; ++i)
			ofile << it->paircounts[i] << '\t';
		ofile << endl;
	}

	ofile.close();

}