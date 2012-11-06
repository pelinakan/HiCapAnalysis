struct InteractionStruct{
	string chr; 
	vector <int > REsites;
	vector < int > position;
	int distance; // Distance from TSS
	int len; // Length of the interaction
	int signal; // the number of reads
	bool type; // 0.promoter-enhancer (PE) or 1.promoter-promoter (PP)
};
struct CenterBinStruct{
			int nofpairs;
			int recount;
			double norm_signal;
};

struct SignalStruct{ // This struct keeps the numbers of restriction fragments where there is at least one pair coming from a feature, one struct per chr
	boost::unordered::unordered_map<int, int> signal_ups; // upstream of the promoter activation direction
	boost::unordered::unordered_map<int, int> signal_down; // downstream of the promoter activation direction
};

struct SignalStruct_CTX{ // This struct keeps the numbers restriction fragments where there is at least one pair coming from a feature, one struct per chr
	boost::unordered::unordered_map<int, int > signal; 
	string maptochrname;
};
struct PeakMap{
	boost::unordered::unordered_set<int> peaks;
	string maptochrname;
};

struct MetaPeakMap{ // Key is the closest REsite to the peak, if there is a particular peak 1 otherwise 0
	boost::unordered::unordered_map<int, string > metapeaks;
	string maptochrname;
};

struct promPromSignalStruct{
	int promoter_index; // within proms struct
	double normsignal; // The number of reads within the core promoter of the interactor promoter
};

struct FeatureStruct{
	vector <int> ProbeIDs;
	string FeatureType; //Promoter, NegCtrl, etc..
	string chr;
	vector < InteractionStruct > AllExperiments_IntBins;
	vector < CenterBinStruct > centerbin;
};

struct GeneStruct: public FeatureStruct {
//Used for Genes
	string RefSeqName; 
	vector <string> TranscriptName;
	vector <int> isoformpromotercoords; 
	int TSS; //Transcription Start Site
	int gene_end; 
	string strand;
	double *expression;
};

struct NegativeControlStruct: public FeatureStruct{ 
//Used for NegativeControls
	string type; //genic or intergenic
	int start; 
	int end; 
	int midpoint;
};

struct PromoterStruct: FeatureStruct{
	vector < string > transcripts;
	vector < string > genes;
	string strand;
	string chr;
	int start;
	int end;
	int TSS;
	int nofRESites; // how many restriction sites the core promoter contains, will be used to normalise the signal
	bool sharedpromoter; 
	vector < string > genes_sharingproms;
	int *closestREsitenums; // offset of the closest RE sites on each side of the promoter this will be used to index interactions 
	//restriction fragments will be numbered with respect to the closest site upstream [0] or downstream [1].
	vector < SignalStruct_CTX > Signals_CTX; // Each element of this vector will represent a chromosome
	SignalStruct Signals; // Intra-chromosomal interactions 
	vector < promPromSignalStruct > promPromSignals; // Each element represent a promoter
};

struct PairStruct{
	string chr_1;
	string chr_2;
	int startcoord;
	int endcoord;
};

