#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cgff_types.h"
#include "gromacs_topology.h"
#include "safe_mem.h"
#include "io_read.h"
#include "wnoid_math.h"
#include "LDD.h"
#include "ematch_types.h"


/* This first group of functions sets up the system */
void initialize_esystem(tW_esystem *sys)
{
  sys->N_Site_Types = 0;
  sys->Site_Types = NULL;
  sys->n_particles = 0;
  sys->n_frames = sys->n_total_bins = 0;
  sys->bond_types = NULL;
  sys->angle_types = NULL; 
  sys->dihedral_types = NULL;
  sys->bond_atom_types = NULL; 
  sys->angle_atom_types = NULL; 
  sys->dihedral_atom_types = NULL;
  sys->n_pair_types = sys->n_bond_types = sys->n_angle_types = sys->n_dihedral_types = 0;
  sys->dr = NULL; 
  sys->db = NULL; 
  sys->da = NULL; 
  sys->dd = NULL;
  sys->dv = NULL;				// RJS 10/07/2020
  sys->r_min = NULL; 
  sys->r_max = NULL;
  sys->b_min = NULL;
  sys->b_max = NULL;
  sys->a_min = NULL;
  sys->a_max = NULL;
  sys->d_min = NULL;
  sys->d_max = NULL;
  sys->v_min = NULL;				// RJS 10/07/2020
  sys->v_max = NULL;				// RJS 10/07/2020
  sys->nb_basis = NULL;
  sys->bond_basis = NULL;
  sys->angle_basis = NULL;
  sys->dihedral_basis = NULL;
  sys->volume_basis = NULL;			// RJS 10/07/2020
  sys->n_nb_pad = NULL;
  sys->n_nb_bins = NULL;
  sys->n_bond_bins = NULL;
  sys->n_angle_bins = NULL;
  sys->n_dihedral_bins = NULL;
  sys->n_volume_bins = NULL;			// RJS 10/07/2020
  sys->nb_i0s = NULL;
  sys->bond_i0s = NULL;
  sys->angle_i0s = NULL;
  sys->dihedral_i0s = NULL;
  sys->volume_i0s = NULL;
  sys->n_bonds = sys->n_angles = sys->n_dihedrals = 0;
  sys->bond_list = NULL;
  sys->angle_list = NULL;
  sys->dihedral_list = NULL;
  sys->A = NULL;
  sys->A_lin = NULL;
  sys->Ev = NULL;
  sys->V = NULL;
  sys->E_times = NULL;
  sys->V_times = NULL;
  sys->bfTraj = sys->bfTop = sys->bfEfile = sys->bfVfile = sys->bfSiteTypes = sys->bfPairs = sys->bfBonds = sys->bfAngles = sys->bfDihedrals = sys->bfVolume = sys->bfRcond = sys->bfCalcmode = FALSE;
  sys->excl_list = NULL;
  sys->nexcla = NULL;

  sys->rcond = 0.00000001; // default
  strcpy(sys->calc_mode,eCM_WHOLE); // default
  sys->bCalcMode_Whole = TRUE;
  sys->bCalcMode_SecondHalf = FALSE;

  sys->fr = init_tW_gmx_trxframe();
  sys->top = init_tW_gmx_topology();
}

void read_etrajectory(tW_esystem *sys, tW_line line)
{
  tW_word directive, fnm;
  int test_sscanf;
    
  test_sscanf = sscanf(line,"%s %s",&directive,&fnm);
    
  if (test_sscanf != 2)
  {
    fprintf(stderr,"ERROR: unable to find trajectory filename following\n");
    fprintf(stderr,"       directive [%s]\n",eKEY_TRAJ);
    fprintf(stderr," line: %s\n",line);
    fprintf(stderr," test_sscanf: %d\n",test_sscanf);
    fprintf(stderr," directive: %s   fnm: %s \n",directive,fnm);
    exit(EXIT_FAILURE);
  }

  read_first_frame(sys->fr,(const char *) fnm);
}

void read_etopology(tW_esystem *sys, tW_line line)
{
  tW_word directive, fnm;
  int test_sscanf;

  test_sscanf = sscanf(line,"%s %s",&directive,&fnm);

  if (test_sscanf != 2)
  {
    fprintf(stderr,"ERROR: unable to find topology filename following\n");
    fprintf(stderr,"       directive [%s]\n",eKEY_TOP);
    fprintf(stderr," line: %s\n",line);
    exit(EXIT_FAILURE);
  }

  if (! read_topology(sys->top,(const char *) fnm))
  {
    fprintf(stderr,"ERROR: unable to read topology file: %s \n",fnm);
    exit(EXIT_FAILURE);
  }
  sys->n_particles = sys->top->contents->atoms.nr;
}

void read_eenergy(tW_esystem *sys, tW_line line)
{
  tW_word directive, fnm;
  int test_sscanf;

  test_sscanf = sscanf(line,"%s %s",&directive,&fnm);

  if (test_sscanf != 2)
  {
    fprintf(stderr,"ERROR: unable to find energy filename following\n");
    fprintf(stderr,"       directive [%s]\n",eKEY_EFILE);
    fprintf(stderr," line: %s\n",line);
    exit(EXIT_FAILURE);
  }

  strcpy(sys->efile_name, fnm);
  FILE *fp = open_file(fnm,'r');
  tW_line inp_line;
  int n_lines = 0;
  
  while (get_next_line(fp,inp_line) != -1)
  {
    ++n_lines;
  }
  rewind(fp);
  sys->n_frames = n_lines;
  sys->Ev = (double *) ecalloc(n_lines, sizeof(double));
  sys->E_times = (double *) ecalloc(n_lines, sizeof(double));

  fprintf(stderr,"n_lines: %d\n",n_lines);

  int i;
  tW_word time, energy;
  for (i = 0; i < n_lines; ++i)
  {
    if (get_next_line(fp,inp_line) == -1) 
    {  
      fprintf(stderr,"ERROR!! ran out of lines in energy file!\n");
      fprintf(stderr,"i: %d   n_lines: %d\n",i,n_lines);
      exit(EXIT_FAILURE);
    }
    test_sscanf = sscanf(inp_line," %s %s ",&time,&energy);
    if (test_sscanf != 2)
    {
      fprintf(stderr,"ERROR: unable to read two strings from line %d in energy file %s\n",i,fnm);
      fprintf(stderr," line: %s\n",inp_line);
      exit(EXIT_FAILURE);
    }    
    sys->E_times[i] = atof(time);
    sys->Ev[i] = atof(energy);
  }
  fclose(fp);   
  
}

// RJS 10/07/2020 
void read_evolume(tW_esystem *sys, tW_line line)
{
  tW_word directive, fnm;
  int test_sscanf;

  test_sscanf = sscanf(line,"%s %s",&directive,&fnm);

  if (test_sscanf != 2)
  {
    fprintf(stderr,"ERROR: unable to find energy filename following\n");
    fprintf(stderr,"       directive [%s]\n",eKEY_VFILE);
    fprintf(stderr," line: %s\n",line);
    exit(EXIT_FAILURE);
  }

  strcpy(sys->vfile_name, fnm);
  FILE *fp = open_file(fnm,'r');
  tW_line inp_line;
  int n_lines = 0;

  while (get_next_line(fp,inp_line) != -1)
  {
    ++n_lines;
  }
  rewind(fp);
  sys->n_frames = n_lines;

  sys->V = (double *) ecalloc(n_lines, sizeof(double));
  sys->V_times = (double *) ecalloc(n_lines, sizeof(double));
  sys->dv = (double *) ecalloc(1, sizeof(double));
  sys->v_min = (double *) ecalloc(1, sizeof(double));
  sys->v_max = (double *) ecalloc(1, sizeof(double));
  sys->volume_basis = (int *) ecalloc(1, sizeof(int));
  sys->n_volume_bins = (int *) ecalloc(1, sizeof(int));
  sys->volume_i0s = (int *) ecalloc(1, sizeof(int));
  fprintf(stderr,"n_lines_volume: %d\n",n_lines);

  int i;
  tW_word time, volume;
  
  sys->dv[i] = (double) 0.1;
  sys->v_min[i] = (double) 0.0;
  sys->v_max[i] = (double) 0.0;
  sys->volume_basis[i] = eEM_DELTA;
  
  for (i = 0; i < n_lines; ++i)
  {
    if (get_next_line(fp,inp_line) == -1)
    {
      fprintf(stderr,"ERROR!! ran out of lines in volume file!\n");
        (n_lines, sizeof(double));
      fprintf(stderr,"i: %d   n_lines: %d\n",i,n_lines);
      exit(EXIT_FAILURE);
    }
    test_sscanf = sscanf(inp_line," %s %s ",&time,&volume);
    if (test_sscanf != 2)
    {
      fprintf(stderr,"ERROR: unable to read two strings from line %d in volume file %s\n",i,fnm);
      fprintf(stderr," line: %s\n",inp_line);
      exit(EXIT_FAILURE);
    }
    sys->V_times[i] = atof(time);
    sys->V[i] = atof(volume);
    
    sys->v_min[0] = sys->V[0];

    if ( sys->V[i] < sys->v_min[0])  
    {      
      sys->v_min[0] = sys->V[i]; 
    }
    if (sys->V[i] > sys->v_max[0])
    {
      sys->v_max[0] = sys->V[i];
    }
  }
  
  sys->n_volume_bins[0] = (int) ((sys->v_max[0] - sys->v_min[0])/sys->dv[0] + 1.01);
     
  fclose(fp);

  fprintf(stderr,"\tvmin: %g   vmax: %g  dv: %g \n",sys->v_min[0],sys->v_max[0],sys->dv[0]);
  fprintf(stderr,"\tn_volume_bins: %d \n",sys->n_volume_bins[0]);
  fprintf(stderr,"vmax-vmin: %g   (vmax-vmin)/dv: %g   + 1.01: %g     (int): %d \n",sys->v_max[0] - sys->v_min[0], (sys->v_max[0]-sys->v_min[0])/sys->dv[0],(sys->v_max[0] - sys->v_min[0]) / sys->dv[0] + 1.01, (int) ((sys->v_max[0] - sys->v_min[0]) / sys->dv[0] + 1.01));
}
// RJS 10/07/2020 

void read_esites(tW_esystem *sys, tW_line inp_line, FILE *fp)
{
  int nst, test_sscanf;
  tW_word dir;

  test_sscanf = sscanf(inp_line," %s %d ",&dir,&nst);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"Expected to find n_site_types after directive [%s]\n",eKEY_SITE_TYPES);
    fprintf(stderr,"  line: %s\n",inp_line);
    exit(EXIT_FAILURE);
  }
  sys->N_Site_Types = nst;

  fprintf(stderr,"nst: %d\n",nst);
 
  sys->Site_Types = (tW_word *) ecalloc(sys->N_Site_Types,sizeof(tW_word));
  int i;
  for (i = 0; i < sys->N_Site_Types; ++i)
  {
    get_next_line(fp,inp_line);
    if ((strstr(inp_line,"End") != NULL) && (strstr(inp_line,eKEY_SITE_TYPES) != NULL))
    {
      fprintf(stderr,"ERROR: you said there are %d site types, but I found \n",sys->N_Site_Types);
      fprintf(stderr,"       [End %s] on line %d\n",eKEY_SITE_TYPES,i);
      exit(EXIT_FAILURE);
    }
    char *newline = strstr(inp_line,"\n");
    char *p = &(inp_line[0]);
    while (((*p) == ' ') || ((*p) == '\t')) { ++p; }
    if ((*newline) == '\n') { (*newline) = '\0'; }
    strcpy(sys->Site_Types[i],p);
    fprintf(stderr,"Site_Types[%d]: %s\n",i,sys->Site_Types[i]);
  }  
  get_next_line(fp,inp_line);
  if ((strstr(inp_line,"End") == NULL) || (strstr(inp_line,eKEY_SITE_TYPES) == NULL))
  {
    fprintf(stderr,"ERROR: expected to find [End %s] after reading %d lines under [%s]\n",eKEY_SITE_TYPES,sys->N_Site_Types,eKEY_SITE_TYPES);
    fprintf(stderr,"       Instead, I found %s\n",inp_line);
    exit(EXIT_FAILURE);
  }
}

void read_epairs(tW_esystem *sys, tW_line inp_line, FILE *fp)
{
  int i, test_sscanf, npt;
  tW_word t1, t2, basis, dir;
  float r0, rc, dr;
  
  test_sscanf = sscanf(inp_line," %s %d ",&dir,&npt);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"Expected to find n_pair_types after directive [%s]\n",eKEY_NB_INTS);
    fprintf(stderr,"  line: %s\n",inp_line);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr,"npt: %d\n",npt);

  sys->n_pair_types = npt;
  sys->dr = (double *) ecalloc(sys->n_pair_types, sizeof(double));
  sys->r_min = (double *) ecalloc(sys->n_pair_types, sizeof(double));
  sys->r_max = (double *) ecalloc(sys->n_pair_types, sizeof(double));
  sys->nb_basis = (int *) ecalloc(sys->n_pair_types, sizeof(int));
  sys->n_nb_pad = (int *) ecalloc(sys->n_pair_types, sizeof(int));
  sys->n_nb_bins = (int *) ecalloc(sys->n_pair_types, sizeof(int));
  sys->nb_i0s = (int *) ecalloc(sys->n_pair_types, sizeof(int));
  sys->pair_atom_types = (tW_word **) ecalloc(sys->n_pair_types, sizeof(tW_word *));
  for (i = 0; i < sys->n_pair_types; ++i)
  {
    sys->pair_atom_types[i] = (tW_word *) ecalloc(2,sizeof(tW_word));
  } 
 
 
  for (i = 0; i < sys->n_pair_types; ++i)
  {
    get_next_line(fp,inp_line);
    if ((strstr(inp_line,"End") != NULL) && (strstr(inp_line,eKEY_NB_INTS) != NULL))
    {
      fprintf(stderr,"ERROR: you said there are %d nb pair types, but I found \n",sys->n_pair_types);
      fprintf(stderr,"       [End %s] on line %d\n",eKEY_NB_INTS,i);
      exit(EXIT_FAILURE);
    }
    test_sscanf = sscanf(inp_line," %s %s %s %f %f %f ",&t1,&t2,&basis,&r0,&rc,&dr);
    if (test_sscanf != 6)
    {
      fprintf(stderr,"ERROR: expected to read 6 items from line under directive [%s]\n",eKEY_NB_INTS);
      fprintf(stderr,"    However, I only read %d items from line: %s\n",test_sscanf,inp_line);
      fprintf(stderr,"    t1: %s  t2: %s  basis: %s  r0: %g  rc: %g  dr: %g \n",t1,t2,basis,r0,rc,dr);
      exit(EXIT_FAILURE);
    } 

    sys->dr[i] = (double) dr;
    sys->r_min[i] = (double) r0;
    sys->r_max[i] = (double) rc;
    strcpy(sys->pair_atom_types[i][0],t1);
    strcpy(sys->pair_atom_types[i][1],t2);
    if (strcmp(basis,DELTA_NAME) == 0) 
    { 
      sys->nb_basis[i] = eEM_DELTA; 
      sys->n_nb_pad[i] = eDELTA_PAD;
    }
    else if (strcmp(basis,LINEAR_NAME) == 0)
    {
      sys->nb_basis[i] = eEM_LINEAR;
      sys->n_nb_pad[i] = eLINEAR_PAD;
    }
    else if (strcmp(basis,BSPLINE_NAME) == 0)
    {
      sys->nb_basis[i] = eEM_BSPLINE;
      sys->n_nb_pad[i] = eBSPLINE_PAD;
    }
    else
    {
      fprintf(stderr,"ERROR: expected third item in line to be the basis type\n");
      fprintf(stderr," line: %s",inp_line);
      fprintf(stderr,"accepted basis types: %s %s %s\n",DELTA_NAME,LINEAR_NAME,BSPLINE_NAME);
      fprintf(stderr,"please fix the line under directive [%s]\n",eKEY_NB_INTS);
      exit(EXIT_FAILURE);
    }
    sys->r_min[i] -= sys->n_nb_pad[i] * sys->dr[i];
    sys->r_max[i] += sys->n_nb_pad[i] * sys->dr[i];
    sys->n_nb_bins[i] = (int) ((sys->r_max[i] - sys->r_min[i])/sys->dr[i] + 1.01);
    fprintf(stderr,"pairidx: %d\n",i);
    fprintf(stderr,"\trmin: %g   rmax: %g  dr: %g \n",sys->r_min[i],sys->r_max[i],sys->dr[i]);
    fprintf(stderr,"\tn_nb_bins: %d   n_nb_pad: %d \n",sys->n_nb_bins[i],sys->n_nb_pad[i]);
    fprintf(stderr,"rmax-rmin: %g   (rmax-rmin)/dr: %g   + 1.01: %g     (int): %d \n",sys->r_max[i] - sys->r_min[i], (sys->r_max[i]-sys->r_min[i])/sys->dr[i],(sys->r_max[i] - sys->r_min[i]) / sys->dr[i] + 1.01, (int) ((sys->r_max[i] - sys->r_min[i]) / sys->dr[i] + 1.01));
  }   
  
  get_next_line(fp,inp_line);
  if ((strstr(inp_line,"End") == NULL) || (strstr(inp_line,eKEY_NB_INTS) == NULL))
  {
    fprintf(stderr,"ERROR: expected to find [End %s] after reading %d lines under [%s]\n",eKEY_NB_INTS,sys->n_pair_types,eKEY_NB_INTS);
    fprintf(stderr,"       Instead, I found %s\n",inp_line);
    exit(EXIT_FAILURE);
  }

}

void read_ebonds(tW_esystem *sys, tW_line inp_line, FILE *fp)
{
  int i, test_sscanf, nbt;
  tW_word t1, t2, basis, dir;
  float r0, rc, dr;

  test_sscanf = sscanf(inp_line," %s %d ",&dir,&nbt);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"Expected to find n_bond_types after directive [%s]\n",eKEY_BOND_INTS);
    fprintf(stderr,"  line: %s\n",inp_line);
    exit(EXIT_FAILURE);
  } 
  fprintf(stderr,"nbt: %d\n",nbt);
  sys->n_bond_types = nbt;

  sys->db = (double *) ecalloc(sys->n_bond_types, sizeof(double));
  sys->b_min = (double *) ecalloc(sys->n_bond_types, sizeof(double));
  sys->b_max = (double *) ecalloc(sys->n_bond_types, sizeof(double));
  sys->bond_basis = (int *) ecalloc(sys->n_bond_types, sizeof(int));
  sys->n_bond_bins = (int *) ecalloc(sys->n_bond_types, sizeof(int));
  sys->bond_i0s = (int *) ecalloc(sys->n_bond_types, sizeof(int));
  sys->bond_atom_types = (tW_word **) ecalloc(sys->n_bond_types, sizeof(tW_word *));


  for (i = 0; i < sys->n_bond_types; ++i)
  {
    sys->bond_atom_types[i] = (tW_word *) ecalloc(2,sizeof(tW_word));
  } 
 
  for (i = 0; i < sys->n_bond_types; ++i)
  {
    get_next_line(fp,inp_line);
    if ((strstr(inp_line,"End") != NULL) && (strstr(inp_line,eKEY_BOND_INTS) != NULL))
    {
      fprintf(stderr,"ERROR: you said there are %d bond types, but I found \n",sys->n_bond_types);
      fprintf(stderr,"       [End %s] on line %d\n",eKEY_BOND_INTS,i);
      exit(EXIT_FAILURE);
    }
    test_sscanf = sscanf(inp_line," %s %s %s %f %f %f ",&t1,&t2,&basis,&r0,&rc,&dr);
    if (test_sscanf != 6)
    {
      fprintf(stderr,"ERROR: expected to read 6 items from line under directive [%s]\n",eKEY_BOND_INTS);
      fprintf(stderr,"    However, I only read %d items from line: %s\n",test_sscanf,inp_line);
      fprintf(stderr,"    t1: %s  t2: %s  basis: %s  r0: %g  rc: %g  dr: %g \n",t1,t2,basis,r0,rc,dr);
      exit(EXIT_FAILURE);
    } 

    fprintf(stderr,"    t1: %s  t2: %s  basis: %s  r0: %g  rc: %g  dr: %g \n",t1,t2,basis,r0,rc,dr);


    sys->db[i] = (double) dr;
    sys->b_min[i] = (double) r0;
    sys->b_max[i] = (double) rc;
    strcpy(sys->bond_atom_types[i][0],t1);
    strcpy(sys->bond_atom_types[i][1],t2);
    if (strcmp(basis,DELTA_NAME) == 0) 
    { 
      sys->bond_basis[i] = eEM_DELTA; 
    }
    else if (strcmp(basis,LINEAR_NAME) == 0)
    {
      sys->bond_basis[i] = eEM_LINEAR;
    }
    else if (strcmp(basis,BSPLINE_NAME) == 0)
    {
      sys->bond_basis[i] = eEM_BSPLINE;
    }
    else
    {
      fprintf(stderr,"ERROR: expected third item in line to be the basis type\n");
      fprintf(stderr," line: %s",inp_line);
      fprintf(stderr,"accepted basis types: %s %s %s\n",DELTA_NAME,LINEAR_NAME,BSPLINE_NAME);
      fprintf(stderr,"please fix the line under directive [%s]\n",eKEY_BOND_INTS);
      exit(EXIT_FAILURE);
    }

    sys->n_bond_bins[i] = (int) ((sys->b_max[i] - sys->b_min[i])/sys->db[i] + 1);
  }   
  
  get_next_line(fp,inp_line);
  if ((strstr(inp_line,"End") == NULL) || (strstr(inp_line,eKEY_BOND_INTS) == NULL))
  {
    fprintf(stderr,"ERROR: expected to find [End %s] after reading %d lines under [%s]\n",eKEY_BOND_INTS,sys->n_bond_types,eKEY_BOND_INTS);
    fprintf(stderr,"       Instead, I found %s\n",inp_line);
    exit(EXIT_FAILURE);
  }

}

void allocate_angle_stuff(tW_esystem *sys, int nat)
{
  int i;
  sys->da = (double *) ecalloc(nat, sizeof(double));
  sys->a_min = (double *) ecalloc(nat, sizeof(double));
  sys->a_max = (double *) ecalloc(nat, sizeof(double));
  sys->angle_basis = (int *) ecalloc(nat, sizeof(int));
  sys->n_angle_bins = (int *) ecalloc(nat, sizeof(int));
  sys->angle_i0s = (int *) ecalloc(nat, sizeof(int));
  sys->angle_atom_types = (tW_word **) ecalloc(nat, sizeof(tW_word *));

  for (i = 0; i < sys->n_angle_types; ++i)
  {
    sys->angle_atom_types[i] = (tW_word *) ecalloc(3,sizeof(tW_word));
  }
  
}

void read_eangles(tW_esystem *sys, tW_line inp_line, FILE *fp)
{
  int i, test_sscanf, nat;
  tW_word t1, t2, t3, basis, dir;
  float r0, rc, dr;

  test_sscanf = sscanf(inp_line," %s %d ",&dir,&nat);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"Expected to find n_angle_types after directive [%s]\n",eKEY_ANGLE_INTS);
    fprintf(stderr,"  line: %s\n",inp_line);
    exit(EXIT_FAILURE);
  }
  


  sys->n_angle_types = nat;
  allocate_angle_stuff(sys,nat);
/*
  sys->da = (double *) ecalloc(nat, sizeof(double));
  sys->a_min = (double *) ecalloc(nat, sizeof(double));
  sys->a_max = (double *) ecalloc(nat, sizeof(double));
  sys->angle_basis = (int *) ecalloc(nat, sizeof(int));
  sys->n_angle_bins = (int *) ecalloc(nat, sizeof(int));
  sys->angle_i0s = (int *) ecalloc(nat, sizeof(int));
  sys->angle_atom_types = (tW_word **) ecalloc(nat, sizeof(tW_word *));

  for (i = 0; i < sys->n_angle_types; ++i)
  {
    sys->angle_atom_types[i] = (tW_word *) ecalloc(3,sizeof(tW_word));
  } 
*/ 

  for (i = 0; i < sys->n_angle_types; ++i)
  {
    get_next_line(fp,inp_line);
    if ((strstr(inp_line,"End") != NULL) && (strstr(inp_line,eKEY_ANGLE_INTS) != NULL))
    {
      fprintf(stderr,"ERROR: you said there are %d angle types, but I found \n",sys->n_angle_types);
      fprintf(stderr,"       [End %s] on line %d\n",eKEY_ANGLE_INTS,i);
      exit(EXIT_FAILURE);
    }
    test_sscanf = sscanf(inp_line," %s %s %s %s %f %f %f ",&t1,&t2,&t3,&basis,&r0,&rc,&dr);
    if (test_sscanf != 7)
    {
      fprintf(stderr,"ERROR: expected to read 7 items from line under directive [%s]\n",eKEY_ANGLE_INTS);
      fprintf(stderr,"    However, I only read %d items from line: %s\n",test_sscanf,inp_line);
      fprintf(stderr,"    t1: %s  t2: %s  t3: %s  basis: %s  a0: %g  ac: %g  da: %g \n",t1,t2,t3,basis,r0,rc,dr);
      exit(EXIT_FAILURE);
    } 

    sys->da[i] = (double) dr;
    sys->a_min[i] = (double) r0;
    sys->a_max[i] = (double) rc;
    strcpy(sys->angle_atom_types[i][0],t1);
    strcpy(sys->angle_atom_types[i][1],t2);
    strcpy(sys->angle_atom_types[i][2],t3);

    if (strcmp(basis,DELTA_NAME) == 0) 
    { 
      sys->angle_basis[i] = eEM_DELTA; 
    }
    else if (strcmp(basis,LINEAR_NAME) == 0)
    {
      sys->angle_basis[i] = eEM_LINEAR;
    }
    else if (strcmp(basis,BSPLINE_NAME) == 0)
    {
      sys->angle_basis[i] = eEM_BSPLINE;
    }
    else
    {
      fprintf(stderr,"ERROR: expected fourth in line to be the basis type\n");
      fprintf(stderr," line: %s",inp_line);
      fprintf(stderr,"accepted basis types: %s %s %s\n",DELTA_NAME,LINEAR_NAME,BSPLINE_NAME);
      fprintf(stderr,"please fix the line under directive [%s]\n",eKEY_ANGLE_INTS);
      exit(EXIT_FAILURE);
    }

    if (sys->a_min[i] < 0.0) { sys->a_min[i] = 0.0; }
    if (sys->a_max[i] > 180.0) { sys->a_max[i] = 180.0; }
    sys->n_angle_bins[i] = (int) ((sys->a_max[i] - sys->a_min[i])/sys->da[i] + 1);
  }   
  
  get_next_line(fp,inp_line);
  if ((strstr(inp_line,"End") == NULL) || (strstr(inp_line,eKEY_ANGLE_INTS) == NULL))
  {
    fprintf(stderr,"ERROR: expected to find [End %s] after reading %d lines under [%s]\n",eKEY_ANGLE_INTS,sys->n_angle_types,eKEY_ANGLE_INTS);
    fprintf(stderr,"       Instead, I found %s\n",inp_line);
    exit(EXIT_FAILURE);
  }

}

void read_edihedrals(tW_esystem *sys, tW_line inp_line, FILE *fp)
{
  int i, test_sscanf, ndt;
  tW_word t1, t2, t3, t4, basis, dir;
  float r0, rc, dr;

  test_sscanf = sscanf(inp_line," %s %d ",&dir,&ndt);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"Expected to find n_dihedral_types after directive [%s]\n",eKEY_DIH_INTS);
    fprintf(stderr,"  line: %s\n",inp_line);
    exit(EXIT_FAILURE);
  }


  sys->n_dihedral_types = ndt;
  sys->dd = (double *) ecalloc(sys->n_dihedral_types, sizeof(double));
  sys->d_min = (double *) ecalloc(sys->n_dihedral_types, sizeof(double));
  sys->d_max = (double *) ecalloc(sys->n_dihedral_types, sizeof(double));
  sys->dihedral_basis = (int *) ecalloc(sys->n_dihedral_types, sizeof(int));
  sys->n_dihedral_bins = (int *) ecalloc(sys->n_dihedral_types, sizeof(int));
  sys->dihedral_i0s = (int *) ecalloc(sys->n_dihedral_types, sizeof(int));
  sys->dihedral_atom_types = (tW_word **) ecalloc(sys->n_dihedral_types, sizeof(tW_word *));
 
  for (i = 0; i < sys->n_dihedral_types; ++i)
  {
    sys->dihedral_atom_types[i] = (tW_word *) ecalloc(4,sizeof(tW_word));
  } 
 
  for (i = 0; i < sys->n_dihedral_types; ++i)
  {
    get_next_line(fp,inp_line);
    if ((strstr(inp_line,"End") != NULL) && (strstr(inp_line,eKEY_DIH_INTS) != NULL))
    {
      fprintf(stderr,"ERROR: you said there are %d dihedral types, but I found \n",sys->n_dihedral_types); 
      fprintf(stderr,"       [End %s] on line %d\n",eKEY_DIH_INTS,i);
      exit(EXIT_FAILURE);
    }
    test_sscanf = sscanf(inp_line," %s %s %s %s %s %f %f %f ",&t1,&t2,&t3,&t4,&basis,&r0,&rc,&dr);
    if (test_sscanf != 8)
    {
      fprintf(stderr,"ERROR: expected to read 8 items from line under directive [%s]\n",eKEY_DIH_INTS);
      fprintf(stderr,"    However, I only read %d items from line: %s\n",test_sscanf,inp_line);
      fprintf(stderr,"    t1: %s  t2: %s  t3: %s  t4: %s  basis: %s  d0: %g  dc: %g  dd: %g \n",t1,t2,t3,t4,basis,r0,rc,dr);
      exit(EXIT_FAILURE);
    } 

    sys->dd[i] = (double) dr;
    sys->d_min[i] = (double) r0;
    sys->d_max[i] = (double) rc;
    strcpy(sys->dihedral_atom_types[i][0],t1);
    strcpy(sys->dihedral_atom_types[i][1],t2);
    strcpy(sys->dihedral_atom_types[i][2],t3);
    strcpy(sys->dihedral_atom_types[i][3],t4);

    if (strcmp(basis,DELTA_NAME) == 0) 
    { 
      sys->dihedral_basis[i] = eEM_DELTA; 
    }
    else if (strcmp(basis,LINEAR_NAME) == 0)
    {
      sys->dihedral_basis[i] = eEM_LINEAR;
    }
    else if (strcmp(basis,BSPLINE_NAME) == 0)
    {
      sys->dihedral_basis[i] = eEM_BSPLINE;
    }
    else
    {
      fprintf(stderr,"ERROR: expected fifth in line to be the basis type\n");
      fprintf(stderr," line: %s",inp_line);
      fprintf(stderr,"accepted basis types: %s %s %s\n",DELTA_NAME,LINEAR_NAME,BSPLINE_NAME);
      fprintf(stderr,"please fix the line under directive [%s]\n",eKEY_DIH_INTS);
      exit(EXIT_FAILURE);
    }

    if (sys->d_min[i] < -180.0) { sys->d_min[i] = -180.0; }
    if (sys->d_max[i] > 180.0) { sys->d_max[i] = 180.0; }
    sys->n_dihedral_bins[i] = (int) ((sys->d_max[i] - sys->d_min[i])/sys->dd[i] + 1);
  }   
  
  get_next_line(fp,inp_line);
  if ((strstr(inp_line,"End") == NULL) || (strstr(inp_line,eKEY_DIH_INTS) == NULL))
  {
    fprintf(stderr,"ERROR: expected to find [End %s] after reading %d lines under [%s]\n",eKEY_DIH_INTS,sys->n_dihedral_types,eKEY_DIH_INTS);
    fprintf(stderr,"       Instead, I found %s\n",inp_line);
    exit(EXIT_FAILURE);
  }

}

void read_ercond(tW_esystem *sys, tW_line inp_line)
{
  tW_word dir, val;
  int test_sscanf;
  
  test_sscanf = sscanf(inp_line,"%s %s",&dir,&val);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"ERROR: expected to find [%s] RCOND\n",eKEY_RCOND);
    fprintf(stderr,"instead, I found line: %s",inp_line);
    fprintf(stderr,"dir: %s   val: %s \n",dir,val);
    exit(EXIT_FAILURE);
  }
  sys->rcond = atof(val);

}

void read_ecalcmode(tW_esystem *sys, tW_line inp_line)
{
  tW_word dir, val;
  int test_sscanf;
  
  test_sscanf = sscanf(inp_line,"%s %s",&dir,&val);
  if (test_sscanf != 2)
  {
    fprintf(stderr,"ERROR: expected to find [%s] CALCMODE\n",eKEY_CALCMODE);
    fprintf(stderr,"instead, I found like: %s",inp_line);
    fprintf(stderr,"dir: %s   val: %s \n",dir,val);
    exit(EXIT_FAILURE);
  }
  strcpy(sys->calc_mode,val);
}


void read_epar(tW_esystem *sys)
{
  FILE *fp = open_file("epar.txt",'r');
  tW_line inp_line;
  
  while (get_next_line(fp,inp_line) != -1)
  {
    if (strstr(inp_line,eKEY_TRAJ) != NULL)
    { 
      if (sys->bfTraj)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_TRAJ);
        exit(EXIT_FAILURE);
      }
      read_etrajectory(sys,inp_line); 
      sys->bfTraj = TRUE;
    }
    
    else if (strstr(inp_line,eKEY_TOP) != NULL)
    {
      if (sys->bfTop)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_TOP);
        exit(EXIT_FAILURE);
      } 
      read_etopology(sys,inp_line); 
      sys->bfTop = TRUE;
    }

    else if (strstr(inp_line,eKEY_EFILE) != NULL)
    { 

      if (sys->bfEfile)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_EFILE);
        exit(EXIT_FAILURE);
      } 
      read_eenergy(sys,inp_line); 
      sys->bfEfile = TRUE;
    }
    
    /* RJS 10/07/2020 */ 
    else if (strstr(inp_line,eKEY_VFILE) != NULL)
    {

      if (sys->bfVfile)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_VFILE);
        exit(EXIT_FAILURE);
      }
      read_evolume(sys,inp_line);
      sys->bfVolume = TRUE;
    }
    /* RJS 10/07/2020 */

    else if (strstr(inp_line,eKEY_SITE_TYPES) != NULL)
    { 
      if (sys->bfSiteTypes)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_SITE_TYPES);
        exit(EXIT_FAILURE);
      } 
      read_esites(sys,inp_line, fp); 
      sys->bfSiteTypes = TRUE; 
    }

    else if (strstr(inp_line,eKEY_NB_INTS) != NULL)
    { 
      if (sys->bfPairs)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_NB_INTS);
        exit(EXIT_FAILURE);
      } 
      read_epairs(sys, inp_line, fp); 
      sys->bfPairs = TRUE;
    }

    else if (strstr(inp_line,eKEY_BOND_INTS) != NULL)
    { 
      if (sys->bfBonds)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_BOND_INTS);
        exit(EXIT_FAILURE);
      }
      double *dummy = (double *) ecalloc(3,sizeof(double)); 
      read_ebonds(sys, inp_line, fp); 
      sys->bfBonds = TRUE;
      double *dummy2 = (double *) ecalloc(3,sizeof(double));
      efree(dummy);
      efree(dummy2);
    }

    else if (strstr(inp_line,eKEY_ANGLE_INTS) != NULL)
    { 
      if (sys->bfAngles)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_ANGLE_INTS);
        exit(EXIT_FAILURE);
      } 
      read_eangles(sys, inp_line, fp); 
      sys->bfAngles = TRUE;
    }

    else if (strstr(inp_line,eKEY_DIH_INTS) != NULL)
    { 
      if (sys->bfDihedrals)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_DIH_INTS);
        exit(EXIT_FAILURE);
      } 
      read_edihedrals(sys, inp_line, fp); 
      sys->bfDihedrals = TRUE;
    }

    else if (strstr(inp_line,eKEY_RCOND) != NULL)
    {
      if (sys->bfRcond)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_RCOND);
        exit(EXIT_FAILURE);
      }
      read_ercond(sys,inp_line);
      sys->bfRcond = TRUE;
    }
    
    else if (strstr(inp_line,eKEY_CALCMODE) != NULL)
    {
      if (sys->bfCalcmode)
      {
        fprintf(stderr,"ERROR: reading directive [%s] twice!\n",eKEY_CALCMODE);
        exit(EXIT_FAILURE);
      }
      read_ecalcmode(sys,inp_line);
      sys->bfCalcmode = TRUE;
    }

    else
    {
      fprintf(stderr,"WARNING: Unrecognized line in epar.txt: %s \n",inp_line);
    }
  }
  fclose(fp);

}

void check_mandatory_stuff(tW_esystem *sys)
{
  if (! sys->bfTraj)
  {
    fprintf(stderr,"ERROR: never found [%s] directive!\n",eKEY_TRAJ);
    exit(EXIT_FAILURE);
  }
  if (! sys->bfTop)
  {
    fprintf(stderr,"ERROR: never found [%s] directive!\n",eKEY_TOP);
    exit(EXIT_FAILURE);
  }
  if (! sys->bfEfile)
  {
    fprintf(stderr,"ERROR: never found [%s] directive!\n",eKEY_EFILE);
    exit(EXIT_FAILURE);
  }
  if (! sys->bfSiteTypes)
  {
    fprintf(stderr,"ERROR: never found [%s] directive!\n",eKEY_SITE_TYPES);
    exit(EXIT_FAILURE);
  }
  if ((! sys->bfPairs) && (! sys->bfBonds) && (! sys->bfAngles) && (! sys->bfDihedrals))
  {
    fprintf(stderr,"ERROR: never found any interaction directives: [%s] [%s] [%s] [%s]\n",
                    eKEY_NB_INTS, eKEY_BOND_INTS, eKEY_ANGLE_INTS, eKEY_DIH_INTS);
    exit(EXIT_FAILURE);
  }
}

void do_stuff(tW_esystem *sys)
{
  int i, j, k;
  int ai, aj, ak, al;
  int nst = sys->N_Site_Types;
  sys->n_bonds = sys->top->get_bond_count(sys->top);
  sys->n_angles = get_nangles(sys->top) / 4; // MRD TODO getting n_bonds, n_angles, and n_dihedrals should be symmetric
  sys->n_dihedrals = get_ndihs(sys->top) / 5;
 
  sys->bond_types = (int *) ecalloc(sys->n_bonds, sizeof(int));
  sys->angle_types = (int *) ecalloc(sys->n_angles, sizeof(int));
  sys->dihedral_types = (int *) ecalloc(sys->n_dihedrals, sizeof(int));
 
  sys->bond_list = get_bond_list(sys->top);
  sys->angle_list = get_angle_list(sys->top);
  sys->dihedral_list = get_dih_list(sys->top);
 
  sys->pair_types = (int **) ecalloc(sys->N_Site_Types, sizeof(int *));
  for (i = 0; i < sys->N_Site_Types; ++i)
  {
    sys->pair_types[i] = (int *) ecalloc(sys->N_Site_Types, sizeof(int));
    for (j = 0; j < sys->N_Site_Types; ++j)
    {
      sys->pair_types[i][j] = -1;
    }
  }

  bool bBroken = FALSE;
// MRD TODO we need to check to make sure pair, bond, angle, and dihedral types arent repeated.
  for (i = 0; i < sys->N_Site_Types; ++i)
  {
    for (j = 0; j < sys->N_Site_Types; ++j)
    {
      for (k = 0; k < sys->n_pair_types; ++k)
      {
        if (((strcmp(sys->top->atom_type_names[i],sys->pair_atom_types[k][0]) == 0) &&
             (strcmp(sys->top->atom_type_names[j],sys->pair_atom_types[k][1]) == 0)) ||
            ((strcmp(sys->top->atom_type_names[i],sys->pair_atom_types[k][1]) == 0) &&
             (strcmp(sys->top->atom_type_names[j],sys->pair_atom_types[k][0]) == 0)))
        {
          sys->pair_types[i][j] = k;
          break;
        }
      } 
    }
  }

  if (sys->bfBonds)
  {
    for (i = 0; i < sys->n_bonds; ++i)
    {
      sys->bond_types[i] = -1;
      ai = sys->bond_list[3*i+1];
      aj = sys->bond_list[3*i+2];
      for (j = 0; j < sys->n_bond_types; ++j)
      {
        if (((strcmp(sys->top->get_atomtype(sys->top,ai),sys->bond_atom_types[j][0]) == 0) &&
             (strcmp(sys->top->get_atomtype(sys->top,aj),sys->bond_atom_types[j][1]) == 0)) ||
            ((strcmp(sys->top->get_atomtype(sys->top,ai),sys->bond_atom_types[j][1]) == 0) &&
             (strcmp(sys->top->get_atomtype(sys->top,aj),sys->bond_atom_types[j][0]) == 0)))
        {
          sys->bond_types[i] = j;
          break;
        }
      }
      if (sys->bond_types[i] == -1)
      {
        fprintf(stderr,"WARNING: unable to determine type of bond idx %d between atoms %d and %d\n",i,ai,aj);
        fprintf(stderr,"atoms are types %s and %s\n",
                        sys->top->get_atomtype(sys->top,ai),
                        sys->top->get_atomtype(sys->top,aj));
        bBroken = TRUE;
      }
    }
  }

  if (sys->bfAngles)
  {
  for (i = 0; i < sys->n_angles; ++i)
  {
    sys->angle_types[i] = -1;
    ai = sys->angle_list[4*i+1];
    aj = sys->angle_list[4*i+2];
    ak = sys->angle_list[4*i+3];
    for (j = 0; j < sys->n_angle_types; ++j)
    {
      if (((strcmp(sys->top->get_atomtype(sys->top,ai),sys->angle_atom_types[j][0]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,aj),sys->angle_atom_types[j][1]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,ak),sys->angle_atom_types[j][2]) == 0)) ||
          ((strcmp(sys->top->get_atomtype(sys->top,ai),sys->angle_atom_types[j][2]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,aj),sys->angle_atom_types[j][1]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,ak),sys->angle_atom_types[j][0]) == 0)))
      {
        sys->angle_types[i] = j;
        break;
      }
    }
    if (sys->angle_types[i] == -1)
    {
      fprintf(stderr,"WARNING: unable to determine type of angle idx %d between atoms %d and %d and %d\n",i,ai,aj,ak);
      fprintf(stderr,"atoms are types %s and %s and %s\n",
                      sys->top->get_atomtype(sys->top,ai),
                      sys->top->get_atomtype(sys->top,aj),
                      sys->top->get_atomtype(sys->top,ak));
      bBroken = TRUE;
    }
  }
  }

  if (sys->bfDihedrals)
  {
  for (i = 0; i < sys->n_dihedrals; ++i)
  {
    sys->dihedral_types[i] = -1;
    ai = sys->dihedral_list[5*i+1];
    aj = sys->dihedral_list[5*i+2];
    ak = sys->dihedral_list[5*i+3];
    al = sys->dihedral_list[5*i+4];
    for (j = 0; j < sys->n_dihedral_types; ++j)
    {
      if (((strcmp(sys->top->get_atomtype(sys->top,ai),sys->dihedral_atom_types[j][0]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,aj),sys->dihedral_atom_types[j][1]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,ak),sys->dihedral_atom_types[j][2]) == 0) && 
           (strcmp(sys->top->get_atomtype(sys->top,al),sys->dihedral_atom_types[j][3]) == 0)) ||
          ((strcmp(sys->top->get_atomtype(sys->top,ai),sys->dihedral_atom_types[j][3]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,aj),sys->dihedral_atom_types[j][2]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,ak),sys->dihedral_atom_types[j][1]) == 0) &&
           (strcmp(sys->top->get_atomtype(sys->top,al),sys->dihedral_atom_types[j][0]) == 0)))
      {
        sys->dihedral_types[i] = j;
        break;
      }
    }
    if (sys->dihedral_types[i] == -1)
    {
      fprintf(stderr,"WARNING: unable to determine type of dihedral idx %d between atoms %d and %d and %d and %d\n",i,ai,aj,ak,al);
      fprintf(stderr,"atoms are types %s and %s and %s and %s\n",
                      sys->top->get_atomtype(sys->top,ai),
                      sys->top->get_atomtype(sys->top,aj),
                      sys->top->get_atomtype(sys->top,ak),
                      sys->top->get_atomtype(sys->top,al));
      bBroken = TRUE;
    }
  }
  }

  if (strcmp(sys->calc_mode,eCM_WHOLE) == 0)
  {
    sys->bCalcMode_Whole = TRUE;
    sys->bCalcMode_SecondHalf = FALSE;
  }
  else if (strcmp(sys->calc_mode,eCM_SECONDHALF) == 0)
  {
    sys->bCalcMode_Whole = FALSE;
    sys->bCalcMode_SecondHalf = TRUE;
  }
  else
  {
    fprintf(stderr,"ERROR: recongized calc modes are: %s %s \n",eCM_WHOLE,eCM_SECONDHALF);
    fprintf(stderr,"You specified calc mode: %s\n",sys->calc_mode);
    bBroken = TRUE;
  }

  if (bBroken)
  {
    exit(EXIT_FAILURE);
  }



  // Checks
  if ((nst * nst + nst)/2 != sys->n_pair_types)
  {
    fprintf(stderr,"WARNING: given %d site types, there are theoretically %d pair interactions\n",nst,(nst*nst+nst)/2);
    fprintf(stderr,"However, epar.txt says there are %d pair interaciton types\n",sys->n_pair_types);
  }
  if (sys->top->get_n_bond_types(sys->top) != sys->n_bond_types)
  {
    fprintf(stderr,"WARNING: topology says there are %d bond types\n",sys->top->get_n_bond_types(sys->top));
    fprintf(stderr,"However, epar.txt says there are %d bond types\n",sys->n_bond_types);
  }
/*
  if (sys->top->get_n_angle_types(sys->top) != sys->n_angle_types)
  {
    fprintf(stderr,"WARNING: topology says there are %d angle types\n",sys->top->get_n_angle_types(sys->top));
    fprintf(stderr,"However, epar.txt says there are %d angle types\n",sys->n_angle_types);
  }
  if (sys->top->get_n_dihedral_types(sys->top) != sys->n_dihedral_types)
  {
    fprintf(stderr,"WARNING: topology says there are %d dihedral types\n",sys->top->get_n_dihedral_types(sys->top));
    fprintf(stderr,"However, epar.txt says there are %d dihedral types\n",sys->n_dihedral_types);
  }
*/
  if ((int)sys->fr->n_frames != sys->n_frames)
  {
    fprintf(stderr,"WARNING: energy file contains data for %d frames\n",sys->n_frames);
    fprintf(stderr,"However, based on the .trr file size, we calculate %d frames in the trajectory file\n",sys->fr->n_frames);
  }

  if (sys->N_Site_Types != sys->top->n_atomtypes)
  {
    fprintf(stderr,"WARNING: epar.txt specifies %d site types\n",sys->N_Site_Types);
    fprintf(stderr,"However, the topology file specifies %d site types\n",sys->top->n_atomtypes);
  }
  else
  {
    bool found;
    for (i = 0; i < sys->N_Site_Types; ++i)
    {
      found = FALSE;
      for (j = 0; j < sys->top->n_atomtypes; ++j)
      {
        if (strcmp(sys->Site_Types[i],sys->top->atom_type_names[j]) == 0)
        {
          found = TRUE;
        }
      }
      if (! found)
      {
        fprintf(stderr,"WARNING: epar.txt site type %s not found in topolgoy file\n",sys->Site_Types[i]);
      }
    }
  }

  // Set some numbers and allocate memory
  int i0 = 0; 
  for (i = 0; i < sys->n_pair_types; ++i)
  {
    sys->nb_i0s[i] = i0;
    i0 += sys->n_nb_bins[i];
    fprintf(stderr,"nbt: %d  i0: %d  n_nb_bins: %d  i0_track: %d \n",i,sys->nb_i0s[i],sys->n_nb_bins[i],i0);
  }
  for (i = 0; i < sys->n_bond_types; ++i)
  {
    sys->bond_i0s[i] = i0;
    i0 += sys->n_bond_bins[i];
    fprintf(stderr,"bt: %d  i0: %d  n_b_bins: %d  i0_track: %d \n",i,sys->bond_i0s[i],sys->n_bond_bins[i],i0);
  }
  for (i = 0; i < sys->n_angle_types; ++i)
  {
    sys->angle_i0s[i] = i0;
    i0 += sys->n_angle_bins[i];
    fprintf(stderr,"at: %d  i0: %d  n_a_bins: %d  i0_track: %d \n",i,sys->angle_i0s[i],sys->n_angle_bins[i],i0);
  }
  for (i = 0; i < sys->n_dihedral_types; ++i)
  {
    sys->dihedral_i0s[i] = i0;
    i0 += sys->n_dihedral_bins[i];
    fprintf(stderr,"dt: %d  i0: %d  n_d_bins: %d  i0_track: %d \n",i,sys->dihedral_i0s[i],sys->n_dihedral_bins[i],i0);
  }

  sys->volume_i0s[0] = i0;
  i0 += sys->n_volume_bins[0];
  fprintf(stderr,"vt: %d  i0: %d  n_vol_bins: %d  i0_track: %d \n",i,sys->volume_i0s[0],sys->n_volume_bins[0],i0);

  sys->n_total_bins = i0;
 
  sys->trim_cols = (bool *) ecalloc(sys->n_total_bins, sizeof(bool));

  fprintf(stderr,"n_total_bins: %d  sys->n_frames: %d   prod: %d \n",sys->n_total_bins, sys->n_frames, sys->n_total_bins * sys->n_frames); 
  sys->A_lin = (double *) ecalloc(sys->n_total_bins * sys->n_frames, sizeof(double));
  sys->A = (double **) ecalloc(sys->n_frames, sizeof(double *));
  for (i = 0; i < sys->n_frames; ++i)
  {
    sys->A[i] = &(sys->A_lin[i*sys->n_total_bins]);
  }

  sys->nexcla = (int *) ecalloc(sys->n_particles, sizeof(int));
  sys->excl_list = (int **) ecalloc(sys->n_particles, sizeof(int *));
  for (i = 0; i < sys->n_particles; ++i)
  {
    sys->nexcla[i] = sys->top->get_nexcl(sys->top,i);
    sys->excl_list[i] = sys->top->get_excl_list(sys->top,i); 
  }

}

void summarize_esystem(tW_esystem *sys)
{
  int i;
  FILE *fp = open_file("summarize_epar.txt",'w');
 
  fprintf(fp,"\n\nTrajectory file: %s\n\n",sys->fr->filename); 
  fprintf(fp,"Topology file: %s\n\n",sys->top->filename);
  fprintf(fp,"Energy file: %s\n\n",sys->efile_name);
  fprintf(fp,"Volume file: %s\n\n",sys->vfile_name);
  fprintf(fp,"n_frames from energy file: %d\n",sys->n_frames);
  fprintf(fp,"n_frames from trajectory file: %d\n\n",(int) sys->fr->n_frames);
  fprintf(fp,"n_total_bins: %d\n\n",sys->n_total_bins);

  fprintf(fp,"n_site_types: %d\n",sys->N_Site_Types);
  for (i = 0; i < sys->N_Site_Types; ++i) { fprintf(fp,"site %d: %s\n",i,sys->Site_Types[i]); }
  fprintf(fp,"\n");

  fprintf(fp,"n_pair_interactions: %d\n",sys->n_pair_types);
  for (i = 0; i < sys->n_pair_types; ++i) 
  { 
    fprintf(fp,"pair int %d: %s %s \n",i,sys->pair_atom_types[i][0],sys->pair_atom_types[i][1]);
    fprintf(fp,"\ti0: %d  n_bins: %d  rmin: %g  rmax: %g dr: %g basis: %s\n",
                sys->nb_i0s[i],sys->n_nb_bins[i],sys->r_min[i],sys->r_max[i],sys->dr[i],
                (sys->nb_basis[i] == eEM_DELTA ? DELTA_NAME : (sys->nb_basis[i] == eEM_LINEAR ? LINEAR_NAME : BSPLINE_NAME)));
  }
  fprintf(fp,"\n");

  fprintf(fp,"n_bond_interactions: %d\n",sys->n_bond_types);
  for (i = 0; i < sys->n_bond_types; ++i) 
  { 
    fprintf(fp,"bond int %d: %s %s \n",i,sys->bond_atom_types[i][0],sys->bond_atom_types[i][1]); 
    fprintf(fp,"\ti0: %d  n_bins: %d  bmin: %g  bmax: %g db: %g basis: %s\n",
                sys->bond_i0s[i],sys->n_bond_bins[i],sys->b_min[i],sys->b_max[i],sys->db[i],
                (sys->bond_basis[i] == eEM_DELTA ? DELTA_NAME : (sys->bond_basis[i] == eEM_LINEAR ? LINEAR_NAME : BSPLINE_NAME)));
  }
  fprintf(fp,"\n");

  fprintf(fp,"n_angle_interactions: %d\n",sys->n_angle_types);
  for (i = 0; i < sys->n_angle_types; ++i) 
  {
    fprintf(fp,"angle int %d: %s %s %s \n",i,sys->angle_atom_types[i][0],sys->angle_atom_types[i][1],sys->angle_atom_types[i][2]); 
    fprintf(fp,"\ti0: %d  n_bins: %d  amin: %g  amax: %g da: %g basis: %s\n",
                sys->angle_i0s[i],sys->n_angle_bins[i],sys->a_min[i],sys->a_max[i],sys->da[i],
                (sys->angle_basis[i] == eEM_DELTA ? DELTA_NAME : (sys->angle_basis[i] == eEM_LINEAR ? LINEAR_NAME : BSPLINE_NAME)));
  }
  fprintf(fp,"\n");

  fprintf(fp,"n_dihedral_interactions: %d\n",sys->n_dihedral_types);
  for (i = 0; i < sys->n_dihedral_types; ++i) 
  { 
    fprintf(fp,"dihedral int %d: %s %s %s %s \n",i,sys->dihedral_atom_types[i][0],sys->dihedral_atom_types[i][1],
                                                   sys->dihedral_atom_types[i][2],sys->dihedral_atom_types[i][3]); 
    fprintf(fp,"\ti0: %d  n_bins: %d  dmin: %g  dmax: %g dd: %g basis: %s\n",
                sys->dihedral_i0s[i],sys->n_dihedral_bins[i],sys->d_min[i],sys->d_max[i],sys->dd[i],
                (sys->dihedral_basis[i] == eEM_DELTA ? DELTA_NAME : (sys->dihedral_basis[i] == eEM_LINEAR ? LINEAR_NAME : BSPLINE_NAME)));
  }
  fprintf(fp,"\n");
  fprintf(fp,"r_cond: %g \n\n",sys->rcond);
  fprintf(fp,"calc_mode: %s\n\n",sys->calc_mode);

  fclose(fp);
}

void setup_ematch_sys(tW_esystem *sys)
{
  initialize_esystem(sys);
  read_epar(sys);
  check_mandatory_stuff(sys);
  do_stuff(sys);
//  summarize_esystem(sys);
}

/* Populate the A matrix */

bool check_excl(int nr_excl, int *excl_list, int j)
{
    int i;

    for (i = 0; i < nr_excl; i++) {
        if (j == excl_list[i]) {
            return TRUE;
        }
    }

    return FALSE;
}


int get_grid_index_for_delta_basis_em(double r, int i_0, double dr, double R_0)
{
    int index, coeff;
    double resid;

    coeff = (int) floor((r - R_0) / dr);        /* ID last grid point before r.                            */
    resid = r - (R_0 + coeff * dr);     /* Distance from prev. grid point.                         */
    if (resid > 0.5 * dr) {
        coeff += 1;
    }                           /* If r closer to next grid_point update the next grid pt. */
    index = i_0 + coeff;

    return index;
}


int get_grid_index_for_linear_basis_em(double r, int i_0, double dr, double R_0)
{
    int index, coeff;

    coeff = (int) floor((r - R_0) / dr + FLOAT_EPS); /* ID last grid point before r. */
    index = i_0 + coeff;        /* index = lower gridpoint */

    if ((r - R_0)/dr + FLOAT_EPS - coeff >= 1.0)
    {
      index = index + 1;
    }
    return index;
}

double calc_linear_spline_A_em(double r, double r0, double dr, int idx, int i0)
{
  return (1.0 - (r - (r0 + dr * (idx - i0)))/dr);
}

void calc_cubic_Bspline_em(double dr, double R, double R_0, int i, double *B)
{

  double s = ((R-R_0) / dr) - (double) i;

  B[3] = s*s*s / 6.0;
  B[2] = (1.0 + 3.0 * (s + s * s - s * s * s))/6.0;
  B[1] = 2.0/3.0 - s * s + 0.5 * s * s * s;
  B[0] = (1.0 - 3.0 * s + 3.0 * s * s - s * s * s) / (6.0);
}


double get_pbc_distance(dvec xi, dvec xj, matrix box)
{
  double dx, dy, dz;
  if (fabs(xi[0] - xj[0]) > box[0][0]/2.0) 
  { 
    dx = box[0][0] - fabs(xi[0] - xj[0]); 
  }
  else 
  { 
    dx = fabs(xi[0] - xj[0]); 
  }
  if (fabs(xi[1] - xj[1]) > box[1][1]/2.0) 
  { 
    dy = box[1][1] - fabs(xi[1] - xj[1]); 
  }
  else 
  { 
    dy = fabs(xi[1] - xj[1]); 
  }
  if (fabs(xi[2] - xj[2]) > box[2][2]/2.0) 
  { 
    dz = box[2][2] - fabs(xi[2] - xj[2]); 
  }
  else 
  { 
    dz = fabs(xi[2] - xj[2]); 
  }
  return sqrt(dx*dx+dy*dy+dz*dz);
}

double get_angle(dvec xi, dvec xj, dvec xk)
{
  double rij, rkj;
  dvec xij, xkj;
  vect_diff(xi,xj,xij);
  vect_diff(xk,xj,xkj);
  rij = calc_norm(xij);
  rkj = calc_norm(xkj);
  return (acos(dot_prod(xij,xkj)/(rij*rkj)) * 180.0 / M_PI);
}

// MRD TODO I'm 99.9% sure this dihedral function is broken
double get_dihedral(dvec xi, dvec xj, dvec xk, dvec xl)
{
  int d;
  double cpn1n2x, cpn1n2y, cpn1n2z, dotprod, len_n_one, len_n_two;
  dvec v_one, v_two, v_three, n_one, n_two;

  dvec n1, n2, b2, m1; // MRD
  double x, y, phi, theta; // MRD

  for (d = 0; d < 3; ++d)
  {
    v_one[d] = xj[d] - xi[d];
    v_two[d] = xk[d] - xj[d];
    v_three[d] = xl[d] - xk[d];
  }

  scal_times_vect(1.0/calc_norm(v_two),v_two,b2); // MRD

  cross_product(v_one, v_two, n_one);
  cross_product(v_two, v_three, n_two);

  scal_times_vect(1.0/calc_norm(n_one),n_one,n1); // MRD
  scal_times_vect(1.0/calc_norm(n_two),n_two,n2); // MRD
  cross_product(n1,b2,m1); // MRD
  x = dot_prod(n1,n2); // MRD
  y = dot_prod(m1,n2); // MRD
  phi = atan2(y,x) * (180.0/M_PI); // MRD


  dotprod = dot_prod(n_one, n_two);
  len_n_one = sqrt(dot_prod(n_one, n_one));
  len_n_two = sqrt(dot_prod(n_two, n_two));

  cpn1n2x = (n_one[1]*n_two[2]) - (n_one[2]*n_two[1]);
  cpn1n2y = (n_one[2]*n_two[0]) - (n_one[0]*n_two[2]);
  cpn1n2z = (n_one[0]*n_two[1]) - (n_one[1]*n_two[0]);

  int num_match = 0;
  if (signbit(cpn1n2x) == signbit(v_two[0]))
  {
    num_match += 1;
  }
  if (signbit(cpn1n2y) == signbit(v_two[1]))
  {
    num_match += 1;
  }
  if (signbit(cpn1n2z) == signbit(v_two[2]))
  {
    num_match += 1;
  }

  if (num_match == 3)
  {
    theta = (acos(dotprod / (len_n_one * len_n_two)) * (180 / M_PI));
    return (acos(dotprod / (len_n_one * len_n_two)) * (180 / M_PI));
  }
  else if (num_match == 0)
  {
    theta = -(acos(dotprod / (len_n_one * len_n_two)) * (180 / M_PI));
    return -(acos(dotprod / (len_n_one * len_n_two)) * (180 / M_PI));
  }
  else
  {
    fprintf(stderr,"There is an error with dihedral.\n");
    fprintf(stderr,"xi: %g %g %g\nxj: %g %g %g\nxk: %g %g %g\nxl: %g %g %g\n",xi[0],xi[1],xi[2],xj[0],xj[1],xj[2],xk[0],xk[1],xk[2],xl[0],xl[1],xl[2]);
    fprintf(stderr,"v1: %g %g %g\nv2: %g %g %g\nv3: %g %g %g\n",v_one[0],v_one[1],v_one[2],v_two[0],v_two[1],v_two[2],v_three[0],v_three[1],v_three[2]);
    fprintf(stderr,"n1: %g %g %g\nn2: %g %g %g\n",n_one[0],n_one[1],n_one[2],n_two[0],n_two[1],n_two[2]);
    fprintf(stderr,"dp: %g ln1: %g ln2: %g \n",dotprod,len_n_one,len_n_two);
    fprintf(stderr,"cp: %g %g %g \n",cpn1n2x,cpn1n2y,cpn1n2z);
    fprintf(stderr,"x: %g  y: %g  phi: %g \n",x,y,phi);
  }
  return 0.0;



}

void do_ematch_bonds(tW_esystem *sys)
{
  int fr_idx = sys->fr->counter ;  
  int i;
  int ai, aj;
  int bond_type;
  double bond_coeffs[4];
  int bond_idxes[4];
  double r;
  for (i = 0; i < sys->n_bonds; ++i)
  {
    ai = sys->bond_list[3*i+1];
    aj = sys->bond_list[3*i+2];
    bond_type = sys->bond_types[i];
    r = get_pbc_distance(sys->fr->contents->x[ai],sys->fr->contents->x[aj],sys->fr->contents->box);
    if (r < sys->b_min[bond_type])
    {
      fprintf(stderr,"ERROR: in frame %d, bond between %d and %d is length %g < %g = b_min\n",
                      fr_idx, ai, aj, r, sys->b_min[bond_type]);
      exit(EXIT_FAILURE);
    }
    else if (r > sys->b_max[bond_type])
    {
      fprintf(stderr,"ERROR: in frame %d, bond between %d and %d is length %g > %g = b_max\n",
                      fr_idx, ai, aj, r, sys->b_max[bond_type]);
      exit(EXIT_FAILURE);
    }
    switch (sys->bond_basis[bond_type])
    {
      case (eEM_DELTA):
        bond_idxes[0] = get_grid_index_for_delta_basis_em(r, sys->bond_i0s[bond_type],
                                            sys->db[bond_type], sys->b_min[bond_type]);
        bond_coeffs[0] = 1.0;
        if ((bond_idxes[0] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[0] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[0]] += bond_coeffs[0]; }
        break;
      case (eEM_LINEAR):
        bond_idxes[0] = get_grid_index_for_linear_basis_em(r, sys->bond_i0s[bond_type], 
                                             sys->db[bond_type], sys->b_min[bond_type]);
        bond_idxes[1] = bond_idxes[0] + 1;
        bond_coeffs[0] = calc_linear_spline_A_em(r,sys->b_min[bond_type],sys->db[bond_type],
                                                    bond_idxes[0], sys->bond_i0s[bond_type]);
        bond_coeffs[1] = 1.0 - bond_coeffs[0];
        if ((bond_idxes[0] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[0] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[0]] += bond_coeffs[0]; }
        if ((bond_idxes[1] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[1] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[1]] += bond_coeffs[1]; }
        break;
      case (eEM_BSPLINE):
        bond_idxes[1] = get_grid_index_for_linear_basis_em(r, sys->bond_i0s[bond_type],
                                             sys->db[bond_type], sys->b_min[bond_type]);
        bond_idxes[0] = bond_idxes[1] - 1;
        bond_idxes[2] = bond_idxes[1] + 1;
        bond_idxes[3] = bond_idxes[1] + 2;
        calc_cubic_Bspline_em(sys->db[bond_type], r, sys->b_min[bond_type], bond_idxes[1] - sys->nb_i0s[bond_type], bond_coeffs);
        if ((bond_idxes[0] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[0] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[0]] += bond_coeffs[0]; }
        if ((bond_idxes[1] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[1] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[1]] += bond_coeffs[1]; }
        if ((bond_idxes[2] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[2] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[2]] += bond_coeffs[2]; }
        if ((bond_idxes[3] < sys->bond_i0s[bond_type] + sys->n_bond_bins[bond_type]) &&
            (bond_idxes[3] >= sys->bond_i0s[bond_type])) { sys->A[fr_idx][bond_idxes[3]] += bond_coeffs[3]; }
        break;
    }
  }
}

void do_ematch_angles(tW_esystem *sys)
{
  int fr_idx = sys->fr->counter ;  
  int i;
  int ai, aj, ak;
  int angle_type;
  double angle_coeffs[4];
  int angle_idxes[4];
  double theta;
  for (i = 0; i < sys->n_angles; ++i)
  {
    ai = sys->angle_list[4*i+1];
    aj = sys->angle_list[4*i+2];
    ak = sys->angle_list[4*i+3];
    angle_type = sys->angle_types[i];
    theta = get_angle(sys->fr->contents->x[ai],
                      sys->fr->contents->x[aj],
                      sys->fr->contents->x[ak]);
    if (theta < sys->a_min[angle_type])
    {
      fprintf(stderr,"ERROR: in frame %d, angle between %d and %d and %d is %g < %g = a_min\n",
                      fr_idx, ai, aj, ak, theta, sys->a_min[angle_type]);
      exit(EXIT_FAILURE);
    }
    else if (theta > sys->a_max[angle_type])
    {
      fprintf(stderr,"ERROR: in frame %d, angle between %d and %d and %d is %g > %g = a_max\n",
                      fr_idx, ai, aj, ak, theta, sys->a_max[angle_type]);
      exit(EXIT_FAILURE);
    }
    switch (sys->angle_basis[angle_type])
    {
      case (eEM_DELTA):
        angle_idxes[0] = get_grid_index_for_delta_basis_em(theta, sys->angle_i0s[angle_type],
                                                 sys->da[angle_type], sys->a_min[angle_type]);
        angle_coeffs[0] = 1.0;

        if ((angle_idxes[0] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[0] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[0]] += angle_coeffs[0]; }
        break;
      case (eEM_LINEAR):
        angle_idxes[0] = get_grid_index_for_linear_basis_em(theta, sys->angle_i0s[angle_type], 
                                                  sys->da[angle_type], sys->a_min[angle_type]);
        angle_idxes[1] = angle_idxes[0] + 1;
        angle_coeffs[0] = calc_linear_spline_A_em(theta, sys->a_min[angle_type],sys->da[angle_type],
                                                          angle_idxes[0],sys->angle_i0s[angle_type]);
        angle_coeffs[1] = 1.0 - angle_coeffs[0];
        if ((angle_idxes[0] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[0] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[0]] += angle_coeffs[0]; }
        if ((angle_idxes[1] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[1] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[1]] += angle_coeffs[1]; }
        break;
      case (eEM_BSPLINE):
        angle_idxes[1] = get_grid_index_for_linear_basis_em(theta, sys->angle_i0s[angle_type],
                                                  sys->da[angle_type], sys->a_min[angle_type]);
        angle_idxes[0] = angle_idxes[1] - 1;
        angle_idxes[2] = angle_idxes[1] + 1;
        angle_idxes[3] = angle_idxes[1] + 2;
        calc_cubic_Bspline_em(sys->da[angle_type], theta, sys->a_min[angle_type], angle_idxes[1] - sys->angle_i0s[angle_type], angle_coeffs);
        if ((angle_idxes[0] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[0] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[0]] += angle_coeffs[0]; }
        if ((angle_idxes[1] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[1] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[1]] += angle_coeffs[1]; }
        if ((angle_idxes[2] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[2] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[2]] += angle_coeffs[2]; }
        if ((angle_idxes[3] < sys->angle_i0s[angle_type] + sys->n_angle_bins[angle_type]) &&
            (angle_idxes[3] >= sys->angle_i0s[angle_type])) { sys->A[fr_idx][angle_idxes[3]] += angle_coeffs[3]; }
        break;
    }
  }
}


void do_ematch_dihedrals(tW_esystem *sys)
{
  int fr_idx = sys->fr->counter ;  
  int i;
  int ai, aj, ak, al;
  int dihedral_type;
  double dihedral_coeffs[4];
  int dihedral_idxes[4];
  double phi;
  for (i = 0; i < sys->n_dihedrals; ++i)
  {
    ai = sys->angle_list[5*i+1];
    aj = sys->angle_list[5*i+2];
    ak = sys->angle_list[5*i+3];
    al = sys->angle_list[5*i+4];
    dihedral_type = sys->dihedral_types[i];
    phi = get_dihedral(sys->fr->contents->x[ai],sys->fr->contents->x[aj],
                       sys->fr->contents->x[ak],sys->fr->contents->x[al]);
    if (phi < sys->d_min[dihedral_type])
    {
      fprintf(stderr,"ERROR: in frame %d, dihedral between %d and %d and %d and %d is %g < %g = d_min\n",
                      fr_idx, ai, aj, ak, al, phi, sys->d_min[dihedral_type]);
      exit(EXIT_FAILURE);
    }
    else if (phi > sys->d_max[dihedral_type])
    {
      fprintf(stderr,"ERROR: in frame %d, dihedral between %d and %d and %d and %d is %g > %g = d_max\n",
                      fr_idx, ai, aj, ak, al, phi, sys->d_max[dihedral_type]);
      exit(EXIT_FAILURE);
    }
    switch (sys->dihedral_basis[dihedral_type])
    {
      case (eEM_DELTA):
        dihedral_idxes[0] = get_grid_index_for_delta_basis_em(phi, sys->dihedral_i0s[dihedral_type],
                                             sys->dd[dihedral_type], sys->d_min[dihedral_type]);
        dihedral_coeffs[0] = 1.0;
        if ((dihedral_idxes[0] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[0] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[0]] += dihedral_coeffs[0]; }
        break;
      case (eEM_LINEAR):
        dihedral_idxes[0] = get_grid_index_for_linear_basis_em(phi, sys->dihedral_i0s[dihedral_type], 
                                                   sys->dd[dihedral_type], sys->d_min[dihedral_type]);
        dihedral_idxes[1] = dihedral_idxes[0] + 1;
        dihedral_coeffs[0] = calc_linear_spline_A_em(phi, sys->d_min[dihedral_type],sys->dd[dihedral_type],
                                                          dihedral_idxes[0],sys->dihedral_i0s[dihedral_type]);
        dihedral_coeffs[1] = 1.0 - dihedral_coeffs[0];
        if ((dihedral_idxes[0] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[0] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[0]] += dihedral_coeffs[0]; }
        if ((dihedral_idxes[1] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[1] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[1]] += dihedral_coeffs[1]; }
        break;
      case (eEM_BSPLINE):
        dihedral_idxes[1] = get_grid_index_for_linear_basis_em(phi, sys->dihedral_i0s[dihedral_type],
                                                   sys->dd[dihedral_type], sys->d_min[dihedral_type]);
        dihedral_idxes[0] = dihedral_idxes[1] - 1;
        dihedral_idxes[2] = dihedral_idxes[1] + 1;
        dihedral_idxes[3] = dihedral_idxes[1] + 2;
        calc_cubic_Bspline_em(sys->dd[dihedral_type], phi, sys->d_min[dihedral_type], dihedral_idxes[1] - sys->dihedral_i0s[dihedral_type], dihedral_coeffs);
        if ((dihedral_idxes[0] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[0] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[0]] += dihedral_coeffs[0]; }
        if ((dihedral_idxes[1] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[1] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[1]] += dihedral_coeffs[1]; }
        if ((dihedral_idxes[2] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[2] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[2]] += dihedral_coeffs[2]; }
        if ((dihedral_idxes[3] < sys->dihedral_i0s[dihedral_type] + sys->n_dihedral_bins[dihedral_type]) &&
            (dihedral_idxes[3] >= sys->dihedral_i0s[dihedral_type])) { sys->A[fr_idx][dihedral_idxes[3]] += dihedral_coeffs[3]; }
        break;
    }
  }
}

/* RJS 10/07/2020 */
void do_ematch_volume(tW_esystem *sys)
{
  int fr_idx = sys->fr->counter ;
  double volume_coeffs[1];
  int volume_idxes[1];
    switch (sys->volume_basis[0])
    {
      case (eEM_DELTA):
        volume_idxes[0] = get_grid_index_for_delta_basis_em(sys->V[fr_idx], sys->volume_i0s[0],
                                            sys->dv[0], sys->v_min[0]);
        volume_coeffs[0] = 1.0;
        if ((volume_idxes[0] < sys->volume_i0s[0] + sys->n_volume_bins[0]) &&
            (volume_idxes[0] >= sys->volume_i0s[0])) { sys->A[fr_idx][volume_idxes[0]] += volume_coeffs[0]; }
        break;
   }
}
/* RJS 10/07/2020 */

void do_ematch_pairs(tW_esystem *sys)
{
  int fr_idx = sys->fr->counter ;  
  int i, j, itype, jtype, pairtype;
  double r;
  int nb_idxes[4];
  double nb_coeffs[4];
  FILE *fp0;
  if (fr_idx == 0) { fp0 = open_file("zero_frame_stuff.dat",'w'); }
  for (i = 0; i < sys->n_particles; ++i)
  {
    itype = get_type(*(sys->top->contents->atoms.atomtype[i]),sys->top);
    for (j = i + 1; j < sys->n_particles; ++j)
    {
      if (check_excl(sys->nexcla[i],sys->excl_list[i],j))
      {
        continue;
      }
      jtype = get_type(*(sys->top->contents->atoms.atomtype[j]),sys->top);
      pairtype = sys->pair_types[itype][jtype];
      if (pairtype != -1)
      {
        r = get_pbc_distance(sys->fr->contents->x[i],sys->fr->contents->x[j],sys->fr->contents->box);
        if ((r >= sys->r_min[pairtype]) && (r < sys->r_max[pairtype]))
        {
          switch (sys->nb_basis[pairtype])
          {
            case (eEM_DELTA):
              nb_idxes[0] = get_grid_index_for_delta_basis_em(r, sys->nb_i0s[pairtype],
                                                 sys->dr[pairtype], sys->r_min[pairtype]);
              nb_coeffs[0] = 1.0;
              if ((nb_idxes[0] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[0] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[0]] += nb_coeffs[0]; }
              break;
            case (eEM_LINEAR):
              nb_idxes[0] = get_grid_index_for_linear_basis_em(r, sys->nb_i0s[pairtype],
                                                         sys->dr[pairtype], sys->r_min[pairtype]);
              nb_idxes[1] = nb_idxes[0] + 1;
              nb_coeffs[0] = calc_linear_spline_A_em(r, sys->r_min[pairtype],sys->dr[pairtype],
                                                                nb_idxes[0],sys->nb_i0s[pairtype]);
              nb_coeffs[1] = 1.0 - nb_coeffs[0];
              if ((nb_idxes[0] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[0] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[0]] += nb_coeffs[0]; }
              if ((nb_idxes[1] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[1] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[1]] += nb_coeffs[1]; }
              break;
            case (eEM_BSPLINE):
              nb_idxes[1] = get_grid_index_for_linear_basis_em(r, sys->nb_i0s[pairtype],
                                                         sys->dr[pairtype], sys->r_min[pairtype]);
              nb_idxes[0] = nb_idxes[1] - 1;
              nb_idxes[2] = nb_idxes[1] + 1;
              nb_idxes[3] = nb_idxes[1] + 2;
              calc_cubic_Bspline_em(sys->dr[pairtype], r, sys->r_min[pairtype], nb_idxes[1] - sys->nb_i0s[pairtype], nb_coeffs);
              if ((nb_idxes[0] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[0] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[0]] += nb_coeffs[0]; }
              if ((nb_idxes[1] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[1] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[1]] += nb_coeffs[1]; }
              if ((nb_idxes[2] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[2] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[2]] += nb_coeffs[2]; }
              if ((nb_idxes[3] < sys->nb_i0s[pairtype] + sys->n_nb_bins[pairtype]) &&
                  (nb_idxes[3] >= sys->nb_i0s[pairtype])) { sys->A[fr_idx][nb_idxes[3]] += nb_coeffs[3]; }
if (fr_idx == 0) { fprintf(fp0,"i: %d j: %d r: %g  nb_idxes: %d %d %d %d nb_coeffs: %g %g %g %g \n",i,j,r,
                   nb_idxes[0],nb_idxes[1],nb_idxes[2],nb_idxes[3],
                   nb_coeffs[0],nb_coeffs[1],nb_coeffs[2],nb_coeffs[3]); }
              break;
          }
        }
      }
    }
  }
  if (fr_idx == 0) { fclose(fp0); }
}

/* Save stuff */
void save_A_matrix(double **A, int nrows, int ncols, char * fnm)
{
  FILE *fp = open_file(fnm,'w');
  char *p = strstr(fnm,".dat");
  if (p != NULL) { strcpy(p,"_hm.dat") ; }
  else { sprintf(fnm,"A_matrix_hm.dat"); }
  FILE *fphm = open_file(fnm,'w');
  int i, j;
  for (i = 0; i < nrows; ++i)
  {
    for (j = 0; j < ncols; ++j)
    {
      fprintf(fphm,"%d %d %g \n",i,j,A[i][j]);
      fprintf(fp,"%g ",A[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fphm);
  fclose(fp);
}

void read_A_matrix(double **A, int nrows, int ncols)
{
  FILE *fp = open_file("A_matrix_hm.dat",'r');
  int i, row, col, iVal, test_sscanf;
  double dVal;
  tW_word val;
  tW_line inp_line;
  for (i = 0; i < nrows * ncols; ++i)
  {
    get_next_line(fp,inp_line);
    test_sscanf = sscanf(inp_line,"%d %d %s",&row,&col,&val);
    if (test_sscanf != 3)
    {
      fprintf(stderr,"ERROR: unable to read 3 values A_matrix_hm.dat line idx %d: %s\n",i,inp_line);
      fprintf(stderr,"row: %d  col: %d  val: %s\n",row,col,val);
      fclose(fp);
      exit(EXIT_FAILURE);
    }
    if (row >= nrows)
    {
      fprintf(stderr,"ERROR: row from A_matrix_hm.dat: %d   from line idx %d: %s\n",row,i,inp_line);
      fprintf(stderr,"nrows based on epar.txt: %d \n",nrows);
      fclose(fp);
      exit(EXIT_FAILURE);
    }
    if (col >= ncols)
    {
      fprintf(stderr,"ERROR: col from A_matrix_hm.dat: %d   from line idx %d: %s\n",col,i,inp_line);
      fprintf(stderr,"nrows based on epar.txt: %d \n",ncols);
      fclose(fp);
      exit(EXIT_FAILURE);
    }
    A[row][col] = atof(val);
  } 
  fclose(fp);
}

void get_zero_cols(double **A, bool *trim_cols, int nrows, int ncols, int *n_zero_cols, tW_esystem *sys)
{
  int i, j;
  (*n_zero_cols) = 0;
  for (i = 0; i < ncols; ++i)
  {
    bool nonzero = FALSE;
    trim_cols[i] = FALSE;
    for (j = 0; j < nrows; ++j)
    {
      if (A[j][i] != 0.0)
      {
        nonzero = TRUE;
        break;
      }
    }
    for (j = 0; j < sys->n_pair_types; ++j)
    {
      if ((sys->nb_i0s[j] <= i) && (sys->nb_i0s[j] + sys->n_nb_pad[j] > i)) { nonzero = FALSE; }
      if ((sys->nb_i0s[j] + sys->n_nb_bins[j] - sys->n_nb_pad[j] <= i) &&
          (sys->nb_i0s[j] + sys->n_nb_bins[j] > i)) {nonzero = FALSE; } 
    }
    if (! nonzero)
    {
      trim_cols[i] = TRUE;
      (*n_zero_cols) += 1;
    }
  }
}

void print_zero_cols(bool *tc, int nt)
{
  FILE *fp = fopen("trimmed_columns.dat","w");
  int i;
  for (i = 0; i < nt; ++i)
  {
    fprintf(fp,"%d %d \n",i,tc[i]);
  }
  fclose(fp);
}

void elim_zero_cols(tW_esystem *sys)
{
  int i, j;

  sys->n_total_bins_trim = sys->n_total_bins - sys->n_zero_cols;
  sys->A_trim = (double **) ecalloc(sys->n_frames, sizeof(double *));
  sys->A_trim_lin = (double *) ecalloc(sys->n_frames * sys->n_total_bins_trim, sizeof(double));
  sys->A_trim_lin_colmajor = (double *) ecalloc(sys->n_frames * sys->n_total_bins_trim, sizeof(double));
   
  for (i = 0; i < sys->n_frames; ++i)
  {
    sys->A_trim[i] = &(sys->A_trim_lin[i*sys->n_total_bins_trim]);
  }

  int tj = 0;
  for (j = 0; j < sys->n_total_bins; ++j)
  {
    if (sys->trim_cols[j] == FALSE)
    {
      for (i = 0; i < sys->n_frames; ++i)
      {
        sys->A_trim[i][tj] = sys->A[i][j];
        sys->A_trim_lin_colmajor[tj * sys->n_frames + i] = sys->A[i][j];
      }
      ++tj;
    }
  }
}

void trim_A(tW_esystem *sys)
{
  get_zero_cols(sys->A,sys->trim_cols,sys->n_frames, sys->n_total_bins, &(sys->n_zero_cols), sys);
  print_zero_cols(sys->trim_cols, sys->n_total_bins);
  elim_zero_cols(sys);
}

void get_answer(tW_esystem *sys)
{
  int i; 
  
 
  int M = sys->n_frames;
  int N = sys->n_total_bins_trim;
  int NRHS = 1;
  double *A = (double *) ecalloc(M*N,sizeof(double));
  int LDA = M;
  double *B = (double *) ecalloc(M, sizeof(double));
  int LDB = M; 
  double *S = (double *) ecalloc(N,sizeof(double));
  double RCOND = sys->rcond;
  int RANK;
  int LWORK = 3 * N + 5 * M;
  double *WORK = (double *) ecalloc(LWORK,sizeof(double));
  int INFO;

  for (i = 0; i < M*N; ++i) { A[i] = sys->A_trim_lin_colmajor[i]; }
  for (i = 0; i < M; ++i) { B[i] = sys->Ev[i]; }
  

  dgelss_(&M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);

  sys->phi = (double *) ecalloc(sys->n_total_bins_trim, sizeof(double));
  for (i = 0; i < N; ++i)
  {
    sys->phi[i] = B[i];
  }
 
  efree(A);
  efree(B);
  efree(S);
  efree(WORK); 
}

void save_answer(tW_esystem *sys)
{
  int i, j, phi_idx;
  char *fnm = (char *) ecalloc(100,sizeof(char)); 
  FILE *fp;
  double r, v;
  phi_idx = 0;
  for (i = 0; i < sys->n_pair_types; ++i)
  {
    sprintf(fnm,"pair_energy_%s_%s.dat",sys->pair_atom_types[i][0],sys->pair_atom_types[i][1]);
    fp = open_file(fnm,'w');
    r = sys->r_min[i];
    
    for (j = sys->nb_i0s[i]; j < sys->nb_i0s[i] + sys->n_nb_bins[i]; ++j)
    {
      if (! sys->trim_cols[j])
      {
        fprintf(fp,"%g %g \n",r,sys->phi[phi_idx]);
        ++phi_idx;
      }
      r += sys->dr[i];
    }
    fclose(fp);
  }

  for (i = 0; i < sys->n_bond_types; ++i)
  {
    sprintf(fnm,"bond_energy_%s_%s.dat",sys->bond_atom_types[i][0],sys->bond_atom_types[i][1]);
    fp = open_file(fnm,'w');
    r = sys->b_min[i];
    for (j = sys->bond_i0s[i]; j < sys->bond_i0s[i] + sys->n_bond_bins[i]; ++j)
    {
      if (! sys->trim_cols[j])
      {
        fprintf(fp,"%g %g \n",r,sys->phi[phi_idx]);
        ++phi_idx;
      }
      r += sys->db[i];
    }
    fclose(fp);
  }

  for (i = 0; i < sys->n_angle_types; ++i)
  {
    sprintf(fnm,"angle_energy_%s_%s_%s.dat",sys->angle_atom_types[i][0],sys->angle_atom_types[i][1],sys->angle_atom_types[i][2]);
    fp = open_file(fnm,'w');
    r = sys->a_min[i];
    for (j = sys->angle_i0s[i]; j < sys->angle_i0s[i] + sys->n_angle_bins[i]; ++j)
    {
      if (! sys->trim_cols[j])
      {
        fprintf(fp,"%g %g \n",r,sys->phi[phi_idx]);
        ++phi_idx;
      }
      r += sys->da[i];
    }
    fclose(fp);
  }

  for (i = 0; i < sys->n_dihedral_types; ++i)
  {
    sprintf(fnm,"dihedral_energy_%s_%s_%s_%s.dat",sys->dihedral_atom_types[i][0],sys->dihedral_atom_types[i][1],
                                                  sys->dihedral_atom_types[i][2],sys->dihedral_atom_types[i][3]);
    fp = open_file(fnm,'w');
    r = sys->d_min[i];
    for (j = sys->dihedral_i0s[i]; j < sys->dihedral_i0s[i] + sys->n_dihedral_bins[i]; ++j)
    {
      if (! sys->trim_cols[j])
      {
        fprintf(fp,"%g %g \n",r,sys->phi[phi_idx]);
        ++phi_idx;
      }
      r += sys->dd[i];
    }
    fclose(fp);
  }

  /* RJS 10/07/2020 */
  if (sys->bfVolume) 
  {
  sprintf(fnm,"volume_energy.dat");
    fp = open_file(fnm,'w');
    v = sys->v_min[0];
    for (j = sys->volume_i0s[0]; j < sys->volume_i0s[0] + sys->n_volume_bins[0]; ++j)
    {
      if (! sys->trim_cols[j])
      {
        fprintf(fp,"%g %g \n",v,sys->phi[phi_idx]);
        ++phi_idx;
      }
      v += sys->dv[0];
    }
    fclose(fp);

  } 
  /* RJS 10/07/2020 */

  efree(fnm);
}

void check_answer(tW_esystem *sys)
{
  int i, j;
  FILE *fp = open_file("Ev_check.dat",'w');
  double Ev;
  double chisq = 0.0;
  for (i = 0; i < sys->n_frames; ++i)
  {
    Ev = 0;
    for (j = 0; j < sys->n_total_bins_trim; ++j)
    {
      Ev += sys->A_trim[i][j] * sys->phi[j];
    }
    chisq += (sys->Ev[i] - Ev) * (sys->Ev[i] - Ev);
    fprintf(fp,"%d %g\n",i,Ev);
  }
  fclose(fp);
  chisq /= sys->n_frames;
  fp = open_file("chisq.dat",'w');
  fprintf(fp,"%g",chisq);
  fclose(fp);
}

void test_stuff()
{
  double dr = 0.01, r0 = 0.0, r = 0.5327438, B[4], A;
  int i = get_grid_index_for_linear_basis_em(r,0,dr,r0);
  A = calc_linear_spline_A_em(r,r0,dr,i,0);
  calc_cubic_Bspline_em(dr,r,r0,i,B);
  
  fprintf(stderr,"i: %d  r0: %g  dr: %g  r: %g  \n",i,r0,dr,r);
  fprintf(stderr,"A: %g  B: %g   A*0.53+B*0.54: %g \n",A,1.0-A,A*0.53+(1-A)*0.54);
  fprintf(stderr,"B0: %g   B1: %g   B2: %g   B3: %g \n",B[0],B[1],B[2],B[3]);
  fprintf(stderr,"B0*0.52+B1*0.53+B2*0.54+B3*0.55: %g \n",B[0]*0.52 + B[1]*0.53 + B[2] * 0.54 + B[3] * 0.55);
  fprintf(stderr,"B3*0.52+B2*0.53+B1*0.54+B0*0.55: %g \n",B[3]*0.52 + B[2]*0.53 + B[1] * 0.54 + B[0] * 0.55);
 



}
