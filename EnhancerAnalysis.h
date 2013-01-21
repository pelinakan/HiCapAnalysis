void EnhancerPositions(void){ // counts the number of inter and intra-genic enhancers

	PromoterClass Promoters;
	RESitesClass dpnIIsites;
	MappabilityClass mapp;

//	mapp.InitialiseVars();
	dpnIIsites.InitialiseVars();
//	Promoters.InitialiseData();
//	Promoters.ReadPromoterAnnotation(dpnIIsites, mapp);
	cout << "Promoters Annotated" << endl;

	ifstream infile("chr2_genesperenhancer.txt");
	ifstream infile2("EnhancerCoords_All.txt");

	ofstream outfile("chr2_enhancerfeats.txt");
	ofstream outfile2("closestREsites_chr2enhancers.txt");
	int coord, gcount,pi; 
	string chr, peakoverlap;
    // chr = "chr2";
/*
	for(int k = 0;k < Promoters.ChrNames_refseq.size(); ++k){
		if(Promoters.ChrNames_refseq[k].compare("chr2") == 0){ // Find the right chromosome (genes are indexed acc. to which chromosome they are on
			pi = k;
			break;
		}
	}
*/
	while(!infile2.eof()){
		//infile >> coord >> gcount >> peakoverlap;
		infile2 >> chr >> coord;
		int *repos;
		repos = new int [2];
		dpnIIsites.GettheREPositions(chr,coord,repos);
		outfile2 << repos[0] << '\t' << repos[1] << endl;
/*
		bool found = 0;
		for(int i = 0; i < Promoters.refseq_indexes[pi].size(); ++i){ //Iterate over all refseq genes on that chromosome
			if(Promoters.refseq[Promoters.refseq_indexes[pi][i]].strand == "+"){
				if((Promoters.refseq[Promoters.refseq_indexes[pi][i]].isoformpromotercoords[0] <= coord && Promoters.refseq[Promoters.refseq_indexes[pi][i]].gene_end >= coord)){ // If the readstart is contained within the core promoter
					outfile << coord << '\t' << gcount << '\t' << peakoverlap << '\t' << Promoters.refseq[Promoters.refseq_indexes[pi][i]].RefSeqName << '\t' << Promoters.refseq[Promoters.refseq_indexes[pi][i]].expression[0] << endl;
					found = 1;
					break;
				}
			}
			else{
				if((Promoters.refseq[Promoters.refseq_indexes[pi][i]].gene_end <= coord && Promoters.refseq[Promoters.refseq_indexes[pi][i]].isoformpromotercoords[0] >= coord)){ // If the readstart is contained within the core promoter
					outfile << coord << '\t' << gcount << '\t' << peakoverlap << '\t' << Promoters.refseq[Promoters.refseq_indexes[pi][i]].RefSeqName << '\t' << Promoters.refseq[Promoters.refseq_indexes[pi][i]].expression[0] << endl;
					found = 1;
					break;
				}

			}
		}
		if (!found)
			outfile << coord << '\t' << gcount << '\t' << peakoverlap << '\t' << "Inter-genic" << '\t' << -1 << endl;
	*/
	}

}