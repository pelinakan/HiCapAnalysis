/*
Background interaction frequency is calculated 
*/
class DetermineBackgroundLevels{
	friend class NegCtrlClass;
public:
	BG_signals bglevels;
	void CalculateMeanandStdRegress(NegCtrlClass&,string);

private:
	void LinearRegression(int,double*,double*,string);

};
void DetermineBackgroundLevels::CalculateMeanandStdRegress(NegCtrlClass& ngs, string BaseFileName){


	boost::unordered::unordered_map< int, int > nofentries_perBin; // required for calculating mean
	boost::unordered::unordered_map< int, int > signal_square; // required for calculating stdev

	//Upstream
	for(int i = 0; i < ngs.NofNegCtrls; i++){
		boost::unordered::unordered_map< int, int >::const_iterator iter;
		for (iter = ngs.negctrls[i].Signals.signal_ups.begin(); iter != ngs.negctrls[i].Signals.signal_ups.end(); ++iter){
			int distance = (abs((ngs.negctrls[i].closestREsitenums[0] - iter->first)));
			int bin = distance / BinSize; 
			if(bglevels.mean_upstream.find(bin) == bglevels.mean_upstream.end())
				bglevels.mean_upstream[bin] = iter->second;
			else
				bglevels.mean_upstream[bin] = bglevels.mean_upstream[bin] + iter->second;
			if(nofentries_perBin.find(bin) == nofentries_perBin.end()){
				nofentries_perBin[bin] = 1;
				signal_square[bin] = (iter->second)*(iter->second);
			}
			else{
				nofentries_perBin[bin] = nofentries_perBin[bin] + 1;
				signal_square[bin] = (signal_square[bin] + ((iter->second)*(iter->second)));
			}
		}
	}
	
// Calculate Mean and stdev
	boost::unordered::unordered_map< int, double>::iterator it; // iterator for bin signals
	for (it = bglevels.mean_upstream.begin(); it != bglevels.mean_upstream.end(); ++it){
		it->second = it->second / (double(nofentries_perBin[it->first])); // Mean for that bin
		double mean_square = ((it->second)*(it->second));
		double signalsquare_mean = double(signal_square[it->first] / double(nofentries_perBin[it->first]));
		bglevels.stdev_upstream[it->first] = sqrt((signalsquare_mean - mean_square));
	}

//Downstream
nofentries_perBin.clear();
signal_square.clear();

	for(int i = 0; i < ngs.NofNegCtrls; i++){
		boost::unordered::unordered_map< int, int >::const_iterator iter;
		for (iter = ngs.negctrls[i].Signals.signal_down.begin(); iter != ngs.negctrls[i].Signals.signal_down.end(); ++iter){
			int  bin = ((abs((iter->first) - ngs.negctrls[i].closestREsitenums[1]))) / BinSize;
			if(bglevels.mean_downstream.find(bin) == bglevels.mean_downstream.end())
				bglevels.mean_downstream[bin] = iter->second;
			else
				bglevels.mean_downstream[bin] = bglevels.mean_downstream[bin] + iter->second;
			if(nofentries_perBin.find(bin) == nofentries_perBin.end()){
				nofentries_perBin[bin] = 1;
				signal_square[bin] = (iter->second)*(iter->second);
			}
			else{
				nofentries_perBin[bin] = nofentries_perBin[bin] + 1;
				signal_square[bin] = signal_square[bin] + ((iter->second)*(iter->second));
			}
		}
	}
// Calculate Mean and stdev
	for (it = bglevels.mean_downstream.begin(); it != bglevels.mean_downstream.end(); ++it){
		it->second = it->second / (double(nofentries_perBin[it->first])); // Mean for that bin
		double mean_square = ((it->second)*(it->second));
		double signalsquare_mean = double(signal_square[it->first] / double(nofentries_perBin[it->first]));
		bglevels.stdev_upstream[it->first] = sqrt((signalsquare_mean - mean_square));
	}
	
	string FileName;
	FileName.append(BaseFileName);
	FileName.append("BackgroundLevels.txt");
	ofstream outf(FileName.c_str());
	
	outf << "UPSTREAM"  << endl;
	for (it = bglevels.mean_upstream.begin(); it != bglevels.mean_upstream.end(); ++it){
		outf << it->first << '\t' << it->second << '\t' << bglevels.stdev_upstream[it->first] << endl;
	}
	outf << "DOWNSTREAM"  << endl;
	for (it = bglevels.mean_downstream.begin(); it != bglevels.mean_downstream.end(); ++it){
		outf << it->first << '\t' << it->second << '\t' << bglevels.stdev_downstream[it->first] << endl;
	}

}
void DetermineBackgroundLevels::LinearRegression(int n, double* x, double* y,string whichstream){

double intercept; //Power law y = (a) * (x^b)
int direction;

Maths::Regression::Linear A(n, x, y);


if (whichstream == "upstream")
	direction = 1;
else
	direction = 0;

	bglevels.b[direction] = A.getSlope();
    cout << "    Slope = " << bglevels.b[direction] << endl;
	
	intercept = A.getIntercept();
	cout << "Intercept = " << intercept << endl << endl;
   
	bglevels.a[direction] = pow(2,intercept);
	
	cout << "Regression coefficient = " << A.getCoefficient() << endl;

}


