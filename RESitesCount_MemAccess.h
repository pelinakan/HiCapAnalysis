class RESitesClass{
	friend class PromoterClass;
public :
	int span;
	vector <string> chr_names;
	boost::unordered::unordered_map< string, int > chroffsets_indexfile;
	vector < REindexes > indexes;
	vector < int > posvector;
	void InitialiseVars(void);
	bool GettheREPositions(std::string, int, int*);
	int GetRESitesCount(string, int , int);
	void CleanClass(void);
private:
};
void RESitesClass::InitialiseVars(void){
#ifdef UNIX
	string s;
	s.append(dirname);
	s.append("mm9.GATCpos.txt");
	ifstream RESitesf(s.c_str());
	s.clear();
#endif
#ifdef WINDOWS
	ifstream RESitesf("C:\\WORK\\3c-SeqCap\\CODES\\3C_Analysis\\HiCapAnalysis\\SupplementaryFiles\\MM9.GATCpos.txt");
#endif

	string chrname, chrp, temp;
	int pos, chrstart; 
	int startpos;
	int span = 1000000; // Window size
//For indexing
	cout << "Initialising RE site Class..." << endl;

bool f=0;
RESitesf >> chrp >> pos >> temp >> temp; // Read the first line
chrname = chrp; // get the chr name outside the loop
chr_names.push_back(chrp);
chrstart = (pos ); // chromosome start
indexes.push_back(REindexes());
chroffsets_indexfile[chrp] = (indexes.size()-1);
while(!(RESitesf.eof())){
	int binstart = (pos);
	int binend = binstart + span;
	startpos = posvector.size();
	while( chrname == chrp && (pos >= binstart && pos <= binend)){
		posvector.push_back(pos + 1); //UCSC coords
		RESitesf >> chrp >> pos >> temp >> temp; 
		if(RESitesf.eof()){
			f = 1; //End of file
			break;
		}
	}
	if(indexes.back().binend.empty())
		indexes.back().binstart.push_back(1);
	else
		indexes.back().binstart.push_back((indexes.back().binend.back())+1);
	indexes.back().binend.push_back((binend + 1));
	indexes.back().offset.push_back(startpos);
	indexes.back().count.push_back((posvector.size()-startpos));
	if((chrname != chrp) || f ){
		indexes.push_back(REindexes());
		chroffsets_indexfile[chrp] = (indexes.size()-1);
		chrname = chrp;
		chrstart = (pos);
		chr_names.push_back(chrp);
	}
}
	cout << "RE site class initialised " << endl;

}
int RESitesClass::GetRESitesCount(string chr, int st, int end){

int bitcount = 0;
int REcount = 0;
int rightchr;
int starttosearch = 0;
int HalfClusterDist = (coreprom_upstream + coreprom_downstream)/2;

boost::unordered::unordered_map< string, int >::iterator it = chroffsets_indexfile.find(chr);
rightchr = it->second; // Get the right index vector
for(int i = 0; i < indexes[rightchr].binstart.size();++i){ // Iterate over the bins
	if (indexes[rightchr].binstart[i] <= st && indexes[rightchr].binend[i] >= st){ // If start is within a bin
		starttosearch = indexes[rightchr].offset[i]; // mark the index in the posvector to start to search
		bitcount = indexes[rightchr].count[i]; // this many elements of the posvector is contained within that bin
		if (((st + HalfClusterDist) > indexes[rightchr].binend[i])) // if end is included in the next bin
			bitcount += indexes[rightchr].count[i+1]; // mark the number of elements in the next bin
			break;
	}
}
for (int i = starttosearch; i < starttosearch + bitcount ; ++i){ //This is where search in the posvector starts
	while ((posvector[i] >= st && posvector[i] <= end)){ // If that RE is within the fragment
		++i; 
		++REcount; // count RE sites
	}
	if(REcount > 0)
		break;
}
//cout << chr << ":" << st << "-" << end << "   " << REcount << endl;
	return REcount;

}

bool RESitesClass::GettheREPositions(std::string chr, int pos, int* renums){ // Returns closest RE sites to a position
	
int HalfClusterDist = (coreprom_upstream + coreprom_downstream)/2;
int rightchr;
int starttosearch = 0;
int bitcount = 0;

boost::unordered::unordered_map< string, int >::iterator it = chroffsets_indexfile.find(chr);
rightchr = it->second; // Get the right index vector
for(int i = 0; i < indexes[rightchr].binstart.size();++i){ // Iterate over the bins
	if (indexes[rightchr].binstart[i] <= pos && indexes[rightchr].binend[i] >= pos){ // If start is within a bin
		if(i > 0){
			starttosearch = indexes[rightchr].offset[i-1]; // mark the index in the posvector to start to search
			bitcount = indexes[rightchr].count[i-1]; // this many elements of the posvector is contained within that bin
			bitcount += indexes[rightchr].count[i];
		}
		else{
			starttosearch = indexes[rightchr].offset[i]; // mark the index in the posvector to start to search
			bitcount = indexes[rightchr].count[i]; // this many elements of the posvector is contained within that bin
		}
		if (((pos + HalfClusterDist) > indexes[rightchr].binend[i]) && (i+1 < indexes[rightchr].binstart.size())) // if end is included in the next bin
			bitcount += indexes[rightchr].count[i+1]; // mark the number of elements in the next bin
		break;
	}
}
	int REposprev, REposnext;
	for (int i = starttosearch; i < starttosearch + bitcount; i+=2){
		REposprev = posvector[i];
		REposnext = posvector[i+1];
		while (REposnext <=  pos){
			REposprev = REposnext;
			++i;
			REposnext = posvector[i];
		}
		renums[0] = REposprev;
		renums[1] = REposnext;
		break;
	}
	return 1;
}


void RESitesClass::CleanClass(void){

	posvector.clear();
	chr_names.clear();
	indexes.clear();
	chroffsets_indexfile.clear();
}