
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: Chromosome Conformation Analysis \n");
  fprintf(stderr, "Contact: Pelin Akan <pelin.akan@scilifelab.se>\n\n");
  fprintf(stderr, "Usage:   AnalyseHiCap ExperimentListFile MinNumberofPairs MinJunctionDistance BaseFileName \n");
  fprintf(stderr, "Parameter 1:  BAMFILE names and details are read from ExperimentListFile: enter the file name \n");
  fprintf(stderr, "Parameter 2:  Minimum Number of Pairs to call an interaction [INT] \n");
  fprintf(stderr, "Parameter 3:  Minimum Junction Distance [INT] \n");
  fprintf(stderr, "Parameter 4:  Base File Name to append output files [string] \n");
  fprintf(stderr, "Uses interactions of promoters to calculate background levels and prints p-values \n");

  return 0;
}
