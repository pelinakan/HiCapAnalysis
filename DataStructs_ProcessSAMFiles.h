struct InteractionStruct{
	string chr;
	int resites[2]; // RE fragment start end
	int pos; //Actual read start
	int distance; // Distance from TSS
	int supp_pairs[2]; // the number of reads / the number of RE sites within the core promoter
	char type; // U: upstream , D: downstream , X: inter-chromosomal
	bool peakoverlap; // if there is a peak in that RE fragment
	string peakprofile; // Peak profile in binary format
	double mappability;
	double p_val[4];
};

struct SignalStruct{ // This struct keeps the numbers of restriction fragments where there is at least one pair coming from a feature, one struct per chr
	boost::unordered::unordered_map<int, int* > signal_ups; // upstream of the promoter activation direction key:REsite, [0]:Read Position, [1..n]: number of times observed for experiments 1..n
	boost::unordered::unordered_map<int, int* > signal_down; // downstream of the promoter activation direction key:REsite, [0]:Read Position, [1..n]: number of times observed for experiments 1..n
//	boost::unordered::unordered_map<int, float> distance; // Distance of the interaction
};

struct SignalStruct_CTX{ // This struct keeps the numbers restriction fragments where there is at least one pair coming from a feature, one struct per chr
	boost::unordered::unordered_map<int, int* > signal; 
	string maptochrname;
};
struct PeakMap{
	boost::unordered::unordered_set<int> peaks;
}; // each struct is a chromosome, chromosome names are mapped to struct indexes

struct MetaPeakMap{ // Key is the closest REsite to the peak, if there is a particular peak 1 otherwise 0
	boost::unordered::unordered_map<int, string > metapeaks;
};
boost::unordered::unordered_map<string, int> MetaPeakChrMap;


struct FeattoFeatSignalStruct{
	int feat_index; // within proms struct
	double normsignal[4]; // The number of reads within the core promoter of the interactor promoter
};

struct FeatureStruct{
	vector <int> ProbeIDs;
	string FeatureType; //Promoter, NegCtrl, etc..
	string chr;
	int nofRESites; // how many restriction sites the core promoter contains, will be used to normalise the signal
	double featmappability;
	int *closestREsitenums; // offset of the closest RE sites on each side of the negctrl this will be used to index interactions 
	int start;
	int end;
	vector < SignalStruct_CTX > Signals_CTX; // Each element of this vector will represent a chromosome
	SignalStruct Signals; // Intra-chromosomal interactions 
	vector < InteractionStruct > interactions;
};

struct GeneStruct: public FeatureStruct {
//Used for Genes
	string RefSeqName; 
	vector <string> TranscriptName;
	vector <int> isoformpromotercoords; 
	int TSS; //Transcription Start Site
	string strand;
	double *expression;
};

struct PromoterStruct: FeatureStruct{
	vector < string > transcripts;
	vector < string > genes;
	double *expression;
	string strand;
	int TSS;
	bool sharedpromoter; 
	vector < string > genes_sharingproms;
	//restriction fragments will be numbered with respect to the closest site upstream [0] or downstream [1].
	vector < FeattoFeatSignalStruct > promPromSignals; // Each element represent a promoter
};
struct NegativeControlStruct: public FeatureStruct{ 
//Used for NegativeControls
	string type; //genic or intergenic
	int midpoint;
	//restriction fragments will be numbered with respect to the closest site upstream [0] or downstream [1].
	vector < FeattoFeatSignalStruct > NegcNegcSignals; // Each element represent a negative control
};
struct PairStruct{
	string chr_1;
	string chr_2;
	int startcoord;
	int endcoord;
	int flagged;
};

//For RE Sites
struct REindexes{
	vector <int> binstart;
	vector <int> binend;
	vector <int> offset;
	vector <int> count;
}; // each index struct will keep a chromosome

struct Mappindexes{
	vector <int> binstart;
	vector <int> binend;
	vector <int> offset;
	vector <int> count;
}; // each index struct will keep a chromosome


struct BG_signals{

	boost::unordered::unordered_map< int, double > mean_upstream;
	boost::unordered::unordered_map< int, double > stdev_upstream;
	boost::unordered::unordered_map< int, double > mean_downstream;
	boost::unordered::unordered_map< int, double > stdev_downstream;

	double a[2]; // a[0] is upstream a[1] is downstream (power law fit)
	double b[2]; // b[0] is upstream b[1] is downstream (power law fit)
};