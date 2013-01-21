class MappabilityClass{

public :
public :
	int span;
	double lowerLimit;
	double dataRange;

	vector <string> chr_names;
	boost::unordered::unordered_map< string, int > chroffsets_indexfile;
	vector < REindexes > indexes;
	vector <int> chr_offsets; //offset bit to start at the right chr
	vector < char > mapp;

	void InitialiseVars(void);
	double GetMappability(string, int, int);
private:
	void FindBinaryPos(int, int, int&, int&, int);
	int CalculateNumberofBitsToRead(int, int&, int&, int, int, int&);


};

void MappabilityClass::InitialiseVars(void){

	int span = 200; // Window size
	int count = 500; // how many windows per chunk
	lowerLimit = 0.0;
	dataRange = 1.0;

#ifdef UNIX
	string s;
	s.append(dirname);
	s.append("crgMapabilityAlign36mer.bedgraph");
	ifstream mappf(s.c_str());
	s.clear();
#endif
#ifdef WINDOWS
	ifstream mappf("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\crgMapabilityAlign36mer.bedgraph");
#endif
	ofstream moutf("MM9.Mappability.200bp.txt");

	string chrname, chrp;
	int startp, endp, chrstart, icount = 0; 
	double mp, mappfull;

	vector < int > starts;
	vector <int > ends;
	vector < double > tm;

	cout << "Initialising Mappability Class..." << endl;

	mappf >> chrp >> startp >> endp >> mp; // Read the first line
	indexes.push_back(REindexes());
	chroffsets_indexfile[chrp] = (indexes.size()-1);
	chrname = chrp; // get the chr name outside the loop
	while(!(mappf.eof())){
		chrstart = startp; // chromosome start
		int binstart = chrstart;
		int binend = binstart + span;
		while(chrname == chrp &&(startp >= binstart && startp <=binend)){
			if (endp <= binend) // it belongs to the next bin also
				tm.push_back((mp*(endp-startp)));				
			else{
				tm.push_back(mp*(binend-startp));
				//check how many bins it contains
				int remainder = endp - binend;
				int rbins = remainder/span;
				double meanm = 0.0;
				if (rbins == 0){
					for(int it=0; it<tm.size();++it)
						meanm += tm[it];
					mappfull = meanm/(double(span));
					int compressed = 127*(mappfull);
					char ascii = char(compressed);
					mapp.push_back(ascii);
					++icount;
					starts.push_back(binstart);
					ends.push_back(binend);
					binstart += span;
					binend += span;
				}
				else{ // if it contains more bins
					for(int it = 0; it < tm.size(); ++it)
						meanm += tm[it];
					mappfull = meanm/(double(span));
					int compressed = 127*(mappfull);
					char ascii = char(compressed);
					mapp.push_back(ascii);
					++icount;
					starts.push_back(binstart);
					ends.push_back(binend);
					binstart += span;
					binend += span;
					for (int b = 1; b <= rbins; ++b){
						mappfull = mp;
						int compressed = 127*(mappfull);
						char ascii = char(compressed);
						mapp.push_back(ascii);
						++icount;
						starts.push_back(binstart);
						ends.push_back(binend);
						binstart += span;
						binend += span;
					}
				}
				tm.clear();
				tm.push_back(mp*(endp-binstart));
			}
			int prevend = endp;
			mappf >> chrp >> startp >> endp >> mp; // Read the first line
			if(mappf.eof())
				break;
			if(chrname == chrp && (startp != prevend)){ //there is a gap
				int gapsize = startp - prevend;
				if( binend >= startp) // If the next line is within the current bin
					tm.push_back(0);
				else{
					//First finish the current bin
					tm.push_back(0);
					double meanm = 0.0;
					for(int it = 0; it < tm.size(); ++it)
						meanm += tm[it];
					mappfull = (meanm/(double(span)));
					int compressed = 127*(mappfull);
					char ascii = char(compressed);
					mapp.push_back(ascii);
					++icount;
					starts.push_back(binstart);
					ends.push_back(binend);
					int gapbins = (startp - binend)/span;
					binstart += span;
					binend += span;
					// If there are more bins in the gap
					for (int g = 1; g <= gapbins; ++g){
						mappfull = 0.0;
						int compressed = 127*(mappfull);
						char ascii = char(compressed);
						mapp.push_back(ascii);
						++icount;
						starts.push_back(binstart);
						ends.push_back(binend);
						binstart += span;
						binend += span;
					}
					tm.clear();
					if(startp >= binstart)
						tm.push_back(0);
				}
			}
			if(icount >= count){
				indexes.back().binstart.push_back(starts[0]);
				indexes.back().binend.push_back(ends.back());
				indexes.back().count.push_back(icount);
				indexes.back().offset.push_back(mapp.size()-icount);

				for(int x = 0; x < starts.size(); ++x){
					moutf << chrname << '\t' << starts[x] << '\t' << ends[x] << '\t' << mapp[x] << endl;
				}
				icount = 0;
				starts.clear();
				ends.clear();

			}
			if( (chrname != chrp) || (mappf.eof()) ){
				indexes.back().binstart.push_back(starts[0]);
				indexes.back().binend.push_back(ends.back());
				indexes.back().count.push_back(icount);
				indexes.back().offset.push_back(mapp.size()-icount);
				icount = 0;
				for(int x = 0; x < starts.size(); ++x){
					moutf << chrname << '\t' << starts[x] << '\t' << ends[x] << '\t' << mapp[x] << endl;
				}
				starts.clear();
				ends.clear();

				indexes.push_back(REindexes());
				chroffsets_indexfile[chrp] = (indexes.size()-1);
				chrname = chrp;
				chrstart = (startp);
				chr_names.push_back(chrp);
			}
		}
	}
	cout << "Mappability Class Initialised " << endl;
}

void  MappabilityClass::FindBinaryPos(int chroffset, int coord, int& fileoffset, int &bitcount, int NumberofBits){

int start, end, count, offset, NumberofBitsRead, i = 0;
	do{
		start = indexes[chroffset].binstart[i];
		end = indexes[chroffset].binend[i];
		offset = indexes[chroffset].offset[i];
		count = indexes[chroffset].count[i];
		++i;
		if( (start <= coord ) && ( end >= coord ) ){
			if(offset != 0)
				fileoffset = (offset + ((coord - start)/span)); // this is the offset
			else
				fileoffset = offset;
             //Check if there is any gaps, determine the bitcount
			 NumberofBitsRead = ((coord - start)/span);
			 bitcount = CalculateNumberofBitsToRead(chroffset, fileoffset, NumberofBitsRead, NumberofBits, coord, end);
			 break;
		}
	}while (i < indexes[chroffset].binstart.size());

}
double MappabilityClass::GetMappability(string chr, int st, int end){

int window = end - st; // RE Fragment size
int NumberofBits = (window/span) + 1; // How many bits to read

double Mappmean = 0;
int bitcount = 0;
int offset = 0;

int rightchr;
boost::unordered::unordered_map< string, int >::iterator it = chroffsets_indexfile.find(chr);
rightchr = it->second; // Get the right index vector

FindBinaryPos(rightchr, st, offset, bitcount, NumberofBits); // find the offset and bitcount

	char* data;
	if(bitcount > 0){
		data = new char[bitcount];
		for (int i = 0; i < bitcount; ++i){
			data [i] = mapp[offset + i ];
			if ( data[i] <128 ){
				double value = lowerLimit+(dataRange*((double(data[i]))/127.0));
				Mappmean += value;
			}
		}
		Mappmean/=bitcount;
	}
	return Mappmean;
}

int MappabilityClass::CalculateNumberofBitsToRead(int chroffset, int& fileoffset, int &bitcount, int NumberofBits, int coord, int &end_prev){
	
	int start, end, count, offset;
	
	start = indexes[chroffset].binstart[fileoffset];
	end = indexes[chroffset].binend[fileoffset];
	offset = indexes[chroffset].offset[fileoffset];
	count = indexes[chroffset].count[fileoffset];
	++fileoffset;
	if ( end_prev !=  start ) // there is a gap, you cannot read more
		return bitcount; // Return to the number of bits that can be read!
	else{
		if ((NumberofBits - bitcount) < count ) // if the remaining bits are all contained in the next step
			return NumberofBits;
		else{
			end_prev = end; // You will read more lines ...
			bitcount += count;
			CalculateNumberofBitsToRead(chroffset, fileoffset, bitcount, NumberofBits, coord, end_prev);
		}
	}
}

