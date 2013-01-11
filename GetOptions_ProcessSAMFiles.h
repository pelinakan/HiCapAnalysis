
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: Chromosome Conformation Analysis \n");
  fprintf(stderr, "Contact: Pelin Akan <pelin.akan@scilifelab.se>\n\n");
  fprintf(stderr, "Usage:   AnalyseHiCap MinNumberofPairs MinJunctionDistance \n");
  fprintf(stderr, "Option:  SAMFILE names and details are read from Experiments.txt file \n");
  fprintf(stderr, "Option:  Minimum Number of Pairs to call an interaction [INT] \n");
  fprintf(stderr, "Option:  Minimum Junction Distance [INT] \n");

  return 0;
}
