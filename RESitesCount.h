class RESitesClass{
	friend class PromoterClass;
public :
	int span;
	int window;

	vector <string> chr_names;
	vector <int> chr_offsets; //offset bit to start at the right chr
	vector < int > chrstarts;
	vector < int > chrends;
	int MappTextBinaryFileSize;
	void InitialiseVars(void);
	void WriteandIndex_RESitesBinaryFile(void);
	void WriteRESitesText_toBinaryFile(void);
	bool GettheREPositions(std::string, int, int*);
	int GetRESitesCount(string, int , int);

private:
	void FindBinaryPos(int, int, int&, int&);
};

void RESitesClass::WriteandIndex_RESitesBinaryFile(void){
	
#ifdef UNIX
	string s;
	s.append(dirname);
	s.append("mm9.GATCpos.txt");
	ifstream RESitesf(s.c_str());
	s.clear();

	s.append(dirname);
	s.append("MM9.GATCPos.Index.txt");
	ofstream REindexf(s.c_str());
	s.clear();

	s.append(dirname);
	s.append("MM9.GATCPos.bin");
	ofstream REbinaryf(s.c_str(), ios::binary); // start end mappability
	s.clear();

	s.append(dirname);
	s.append("MM9.chrLengths.txt");
	ofstream REchrbounds(s.c_str()); // start end mappability
#endif

#ifdef WINDOWS
	ifstream RESitesf("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCpos.txt");

	ofstream REindexf("MM9.GATCPos.index.txt");// Format : chr '\t' start '\t' offset
	ofstream REbinaryf("MM9.GATCPos.bin", ios::binary); // start end mappability
	ofstream REchrbounds("MM9.chrLengths.txt");
#endif

	string chrname,chrp,temp;
	int pos, chrstart; 
	int span = 1000000; // Window size
	int offset;
	vector < int > posvector;

	//For indexing

	bool f=0;
	RESitesf >> chrp >> pos >> temp >> temp; // Read the first line
	chrname = chrp; // get the chr name outside the loop
	chrstart = pos; // chromosome start
	
	while(!(RESitesf.eof())){
		int binstart = pos;
		int binend = binstart + span;
		while( chrname == chrp && (pos >= binstart && pos <= binend)){
			posvector.push_back(pos);
			offset = REbinaryf.tellp(); // To index the binary file that contains number and pos of RE sites
			RESitesf >> chrp >> pos >> temp >> temp; 
			if(RESitesf.eof()){
				f = 1; //End of file
				break;
			}
		}
		for(int i=0; i < posvector.size();++i){
			REbinaryf.write((char *)(&(posvector[i])), sizeof((posvector[i])));	
		}	
		REindexf << chrname << '\t' << binstart << '\t' << binend << '\t'
			<< span << '\t' << (posvector.size()) << '\t' << offset << endl;
		
		if((chrname != chrp) || f ){			
			REchrbounds << chrname << '\t' << chrstart << '\t' << posvector.back() << endl;
			chrname = chrp;
			chrstart = pos;
			cout << chrname <<  endl;
		}
		posvector.clear();
	}
}
void RESitesClass::WriteRESitesText_toBinaryFile(void){ //Index the index file itself
string temp, chr, chrtemp;

int start, end, count, offset, fsize, fsize_temp;
#ifdef WINDOWS
	ifstream MappF("MM9.GATCPos.Index.txt");
	ofstream MappFBin("MM9.GATCPos.Index.binary", ios::binary);
	ofstream ofile("MM9.GATCPos.Index.chr_offsets.txt");

#endif
#ifdef UNIX
string s;
s.append(dirname);
s.append("MM9.GATCPos.Index.txt");
ifstream MappF(s.c_str());
s.clear();
s.append(dirname);
s.append("MM9.GATCPos.Index.binary");
ofstream MappFBin(s.c_str(),ios::binary);
s.clear();
s.append(dirname);
s.append("MM9.GATCPos.Index.chr_offsets.txt");
ofstream ofile(s.c_str());
#endif

// Mapptext_MM9.binary file contains start, end count and offset information

	fsize = MappFBin.tellp();
	chr_offsets.push_back(fsize);
	MappF >> chr >> start >> end >> temp >> count >> offset;
	chr_names.push_back(chr);
	getline(MappF,temp);
	MappFBin.write((char *)(&start), sizeof(start));
	MappFBin.write((char *)(&end), sizeof(end));
	MappFBin.write((char *)(&offset),sizeof(offset));
	MappFBin.write((char *)(&count), sizeof(count));
	do{
		do{
			fsize_temp = MappFBin.tellp();
			MappF >> chrtemp >> start >> end >> temp >> count >> offset;
			getline(MappF,temp);
			if(chrtemp != chr){
				chr_names.push_back(chrtemp);
				chr_offsets.push_back(fsize_temp);
				chr = chrtemp;
			}		
			fsize = MappFBin.tellp();
			MappFBin.write((char *)(&start), sizeof(start));
			MappFBin.write((char *)(&end), sizeof(end));
			MappFBin.write((char *)(&offset),sizeof(offset));
			MappFBin.write((char *)(&count), sizeof(count));
		}while(chr == chrtemp && (!MappF.eof()));
		chr_names.push_back(chrtemp);
		chr_offsets.push_back(fsize);
		chr = chrtemp;
	}while(!MappF.eof());
	MappFBin.close();

	for(int i=0;i<chr_names.size();++i)
		ofile << chr_names[i] << '\t' << chr_offsets[i] << endl;
}
void RESitesClass::InitialiseVars(void){
	span = 200;
	window = 200000;

#ifdef WINDOWS
	ifstream ifile("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCPos.Index.chr_offsets.txt");
	ifstream Mapptextbin("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCPos.Index.binary", ios::binary);
#endif

#ifdef UNIX
string s;
s.append(dirname);
s.append("MM9.GATCPos.Index.chr_offsets.txt");
ifstream ifile(s.c_str());
s.clear();
s.append(dirname);
s.append("MM9.GATCPos.Index.binary");
ifstream Mapptextbin(s.c_str(), ios::binary);
#endif

	Mapptextbin.seekg(0, ios::end);
	MappTextBinaryFileSize = Mapptextbin.tellg();
	Mapptextbin.close();

	string chrname;
	int chroffset;
	do{
		ifile >> chrname >> chroffset;
		chr_names.push_back(chrname);
		chr_offsets.push_back(chroffset);
	}while(!ifile.eof());
	ifile.close();

#ifdef WINDOWS
	ifstream chrbounds("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.chrLengths.txt");
#endif
#ifdef UNIX
s.clear();
s.append(dirname);
s.append("MM9.chrLengths.txt");
ifstream chrbounds(s.c_str());
#endif
	int st, end;
	chrstarts.resize(chr_names.size());
	chrends.resize(chr_names.size());
	do{
		chrbounds >> chrname >> st >> end;
		for (int i = 0; i < chr_names.size();++i){
			if(chrname == chr_names[i]){
				chrstarts[i] = st;
				chrends[i] = end;
				break;
			}
		}
	}while(!chrbounds.eof());
	chrbounds.close();

	cout << "RE site class initialised " << endl;


}
void  RESitesClass::FindBinaryPos(int chroffset, int coord, int &fileoffset, int &bitcount){

#ifdef WINDOWS
ifstream file2 ("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCPos.Index.binary",ios::binary);
#endif

#ifdef UNIX
string s;
s.append(dirname);
s.append("MM9.GATCPos.Index.binary");
ifstream file2 (s.c_str(),ios::binary);
// GCtext_MM9.binary file contains start, end count and offset information
#endif

int start, end, count, offset;

    file2.seekg(chroffset, ios::beg); // Get to the correct chromosome position 
	do{
		file2.read((char *)(&start),sizeof(start));
		file2.read((char *)(&end),sizeof(end));
		file2.read((char *)(&offset),sizeof(offset));
		file2.read((char *)(&count),sizeof(count));
		
		if(( start <= coord ) && ( end >= coord )){
				fileoffset = offset; // this is the offset
				bitcount = count;
				if (((coord + HalfClusterDist) > end)){
					file2.read((char *)(&start),sizeof(start));
					file2.read((char *)(&end),sizeof(end));
					file2.read((char *)(&offset),sizeof(offset));
					file2.read((char *)(&count),sizeof(count));
					
					bitcount += count;
				}
				break;
		}
	}while (!(file2.eof()));

	file2.close();

}
int RESitesClass::GetRESitesCount(string chr, int st, int end){

#ifdef WINDOWS
ifstream file ("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCPos.bin", ios::in | ios::binary | ios::ate);
#endif

#ifdef UNIX
string s;
s.append(dirname);
s.append("MM9.GATCPos.bin");
ifstream file(s.c_str(),ios::binary);
#endif

ifstream::pos_type fileSize;

int bitcount = 0;
int REcount = 0;
if(file.is_open())
{
	int offset;
	int chroffset,chrend;
	int REpos = 0;


	for(int i=0;i<chr_names.size();++i){
		if(chr_names[i] == chr){
			if ((end <= chrstarts[i] ) || st >= chrends[i] )
				return 0;
			if (st <= chrstarts[i] && (end >= chrstarts[i]))
				st = chrstarts[i];
			chroffset = chr_offsets[i];
			chrend = chrends[i];
		}
	}	
	FindBinaryPos(chroffset, st, offset, bitcount);
	file.seekg(offset, ios::beg);
	for (int i = 0; i < bitcount ; ++i){
		if (file.eof())
			break;
		file.read((char *)(&REpos),sizeof(REpos));

		while ((REpos >= st && REpos <= end)){
			if (file.eof())
				break;
			file.read((char *)(&REpos),sizeof(REpos));
			++REcount;
		}
		if(REcount > 0)
			break;
	}
	file.close();
}
//cout << chr << ":" << st << "-" << end << "   " << REcount << endl;
	return REcount;

}

bool RESitesClass::GettheREPositions(std::string chr, int pos, int* renums){
	
#ifdef WINDOWS
	ifstream file ("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCPos.bin", ios::in | ios::binary | ios::ate);
#endif
	
#ifdef UNIX
string s;
s.append(dirname);
s.append("MM9.GATCPos.bin");
ifstream file(s.c_str(),ios::binary);
#endif
	ifstream::pos_type fileSize;

	if(file.is_open())
	{
		int chroffset,chrend;
		
		// Find the right chromosome
		for(int i=0;i<chr_names.size();++i){
			if(chr_names[i] == chr){
				if ((pos <= chrstarts[i] ) || pos >= chrends[i] )// if outside the chromosome limits, return 0
					return 0;
				if (pos <= chrstarts[i] && (pos >= chrstarts[i]))
					pos = chrstarts[i];
				
				chroffset = chr_offsets[i];
				chrend = chrends[i];
			}
		}	
		// Get to the right offset
#ifdef WINDOWS
		ifstream file2 ("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCPos.Index.binary",ios::binary);
		// GATC_MM9.binary file contains start, end count and offset information
#endif
#ifdef UNIX
string s;
s.append(dirname);
s.append("MM9.GATCPos.Index.binary");
ifstream file2(s.c_str(),ios::binary);
#endif
		file2.seekg(chroffset, ios::beg); // Get to the correct chromosome position 
		
		int start, end;
		int fileoffset, offsetprev, offsetcurrent, offsetnext;
		int bitcount = 0, countprev = 0, countcurrent, countnext;
		
		file2.read((char *)(&start),sizeof(start));
		file2.read((char *)(&end),sizeof(end));
		file2.read((char *)(&offsetprev),sizeof(offsetprev));
		file2.read((char *)(&countprev),sizeof(countprev));
		do{
			file2.read((char *)(&start),sizeof(start));
			file2.read((char *)(&end),sizeof(end));
			file2.read((char *)(&offsetcurrent),sizeof(offsetcurrent));
			file2.read((char *)(&countcurrent),sizeof(countcurrent));
			
			if(( start <= pos ) && ( end >= pos)){
				fileoffset = offsetcurrent; // this is the offset
				bitcount = countprev;
				bitcount += countcurrent;
			
				file2.read((char *)(&start),sizeof(start));
				file2.read((char *)(&end),sizeof(end));
				file2.read((char *)(&offsetnext),sizeof(offsetnext));
				file2.read((char *)(&countnext),sizeof(countnext));
				
				bitcount += countnext;
				break;
			}
		offsetcurrent = offsetprev;
		countcurrent = countprev;
		}while (!(file2.eof()));
		file2.close(); // Offsets are read now read the restriction enzyme sites

	int REposprev, REposnext;
	file.seekg(offsetprev, ios::beg); // Get to the right position
	for (int i = 0; i < bitcount; ++i){
		file.read((char *)(&REposprev),sizeof(REposprev));
		file.read((char *)(&REposnext),sizeof(REposnext));
		while (REposnext <=  pos && pos < chrend){
			REposprev = REposnext;
			file.read((char *)(&REposnext),sizeof(REposnext));
		}
		renums[0] = REposprev;
		renums[1] = REposnext;
		break;
	}
	file.close();
	}
	return 1;
}

