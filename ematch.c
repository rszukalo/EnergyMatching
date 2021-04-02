/*
 *  Here's a template for analysis/manipulation of stuff
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "safe_mem.h"
#include "io_read.h"
#include "gromacs_topology.h"
#include "LDD.h"
#include "ematch_types.h"
#include "ematch_func.h"

int main(int argc, char * argv[])
{
  /* Generic indices */
  int i, j, k;

  tW_esystem sys; 
  setup_ematch_sys(&sys);
  summarize_esystem(&sys);


  /* loop over frames */
  if (sys.bCalcMode_Whole)
  {
    while (read_next_frame(sys.fr,TRUE))
    {
      if (sys.bfBonds) { do_ematch_bonds(&sys); } 
      if (sys.bfAngles) { do_ematch_angles(&sys); }
      if (sys.bfDihedrals) { do_ematch_dihedrals(&sys); }
      if (sys.bfPairs) { do_ematch_pairs(&sys); }
      if (sys.bfVolume) { do_ematch_volume(&sys); }     
    }
  }
  fclose(sys.fr->fp);

  
  char *fnm = (char *) ecalloc(50, sizeof(char));
  if (sys.bCalcMode_Whole)
  {
    sprintf(fnm,"A_matrix.dat");
    save_A_matrix(sys.A, sys.n_frames, sys.n_total_bins, fnm);
  }
  else if (sys.bCalcMode_SecondHalf)
  {
    read_A_matrix(sys.A, sys.n_frames, sys.n_total_bins);
  }

  trim_A(&sys);
  sprintf(fnm,"A_matrix_trimmed.dat");
  save_A_matrix(sys.A_trim, sys.n_frames, sys.n_total_bins_trim, fnm);

  get_answer(&sys);
  save_answer(&sys);  
  check_answer(&sys);


  return 0;
}
