
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: Chromosome Conformation Analysis \n");
  fprintf(stderr, "Contact: Pelin Akan <pelin.akan@scilifelab.se>\n\n");
  fprintf(stderr, "Usage:   HiCapAnalysis inputfile  OutputFileNameBase CellType \n");
  fprintf(stderr, "Option:  inputfile [string] (input SAM file \n"); 
  fprintf(stderr, "Option:  OutputFileNameBase [string] (alloutput files will start with this name \n");
  fprintf(stderr, "Option:  CellType for expression data [INT]: 0:mES, 1:XEN, 2:TS \n");
  return 0;
}
