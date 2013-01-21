class LaminB1Class{
	friend class PromoterClass;
public : 
	vector <string> chr_names; // name of chrs
	boost::unordered::unordered_map< string, int > chr_indexes; // the start index of chr in the laminvals struct
	struct LaminVals{
		int start;
		int end;
		double laminval;
	};

	void CalculateGeneLaminValues(PromoterClass&);
	LaminVals* lv;
private:
};

void LaminB1Class::CalculateGeneLaminValues(PromoterClass& proms){

	LaminVals* lv = new LaminVals[2000000];

#ifdef UNIX
	string s;
	s.append(dirname);
	s.append("laminB1_ES.txt");
	ifstream infile(s.c_str());
	s.clear();
#endif
#ifdef WINDOWS
	ifstream infile("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\laminB1_ES.txt");
#endif

	//	    
	//	607	chr1	3000298	3000358	chr1.0	60	1	0	/gbdb/mm9/wib/laminB1_ES.wib	1.398	0	1	1.398	1.9544
	//	607	chr1	3001360	3002620	chr1.1	60	21	1	/gbdb/mm9/wib/laminB1_ES.wib	0.006	0.026	2	0.038	0.00106
	//	607	chr1	3003790	3003850	chr1.2	60	1	22	/gbdb/mm9/wib/laminB1_ES.wib	-0.32	0	1	-0.32	0.1024
	//	607	chr1	3004360	3004420	chr1.3	60	1	23	/gbdb/mm9/wib/laminB1_ES.wib	-0.48	0	1	-0.48	0.2304
	//	607	chr1	3011410	3012610	chr1.4	60	20	24	/gbdb/mm9/wib/laminB1_ES.wib	1.529	0.092	2	3.15	4.96548
	string temp, chr, chrname;
	int st, end,index = 0;
	double sum, count;
	//For indexing
	cout << "Reading LaminB1_ES file..." << endl;

	infile >> temp >> chr >> st >> end >> temp >> temp >> count >> temp >> temp >> temp >> temp >> temp;
	infile >> sum >> temp; // Read the first line
	chrname = chr; // get the chr name outside the loop
	chr_names.push_back(chr);
	chr_indexes[chr] = (0);
	lv[index].start = st;
	lv[index].end = end;
	lv[index].laminval  = sum/count;
	++index;

	while(!(infile.eof())){
		infile >> temp >> chr >> st >> end >> temp >> temp >> count >> temp >> temp >> temp >> temp >> temp;
		infile >> sum >> temp; // Read the first line
		lv[index].start = st;
		lv[index].end = end;
		lv[index].laminval  = sum/count;
		++index;

		if((chrname != chr) ){
			cout << chrname << endl;
			chrname = chr; // get the chr name outside the loop
			chr_names.push_back(chr);
			chr_indexes[chr] = (index - 1);
		}
	}
	cout << "LaminB1 values read " << endl;
//////////////////////////////////////////////

	ofstream ofile("LaminB1_mES_RefSeqGenes.txt");
	int gstart, gend;
	double genelamin = 0;
	for(int i = 0; i < proms.NofGenes; ++i){
		genelamin = 0.0;		count = 0.0;
		if (proms.refseq[i].strand == "+"){
			gstart = proms.refseq[i].isoformpromotercoords[0];
			gend = proms.refseq[i].end;
		}
		else{
			gstart = proms.refseq[i].end;
			gend = proms.refseq[i].isoformpromotercoords[0];
		}
		boost::unordered::unordered_map<string, int>::iterator citer = chr_indexes.find(proms.refseq[i].chr);
		if (chr_indexes.find(proms.refseq[i].chr) != chr_indexes.end()){ // Find the right peak map wrt chromosome
			index  = citer->second;
			while (lv[index].start < lv[index+1].start){ // within the same chromosome
				if ((lv[index].start >= gstart && lv[index].end <= gend)){
					genelamin += lv[index].laminval;
					++count;
				}
				++index;
			}
		}
		ofile << proms.refseq[i].RefSeqName << '\t' << proms.refseq[i].chr << ":" << gstart << "-" << gend << '\t'
			  << proms.refseq[i].expression[0] << '\t' << genelamin/count << endl;

	}
}

