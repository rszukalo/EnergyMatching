#ifndef NOID_ETYPES
#define NOID_ETYPES

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdint.h>
#include "cgff_types.h"

//Basis sets
enum {eEM_DELTA, eEM_LINEAR, eEM_BSPLINE};

#define eKEY_TRAJ "Trajectory"
#define eKEY_TOP "Topology"
#define eKEY_EFILE "Energy_File"
#define eKEY_VFILE "Volume_File"			// RJS 10/07/2020
#define eKEY_SITE_TYPES "Site_Types"
#define eKEY_NB_INTS "Pair_Interactions"
#define eKEY_BOND_INTS "Bond_Interactions"
#define eKEY_ANGLE_INTS "Angle_Interactions"
#define eKEY_DIH_INTS "Dihedral_Interactions"
#define eKEY_RCOND "Condition_Number"
#define eKEY_CALCMODE "Calculation_Mode"
#define eCM_WHOLE "Whole"
#define eCM_SECONDHALF "Second_Half"

#define eDELTA_PAD 1
#define eLINEAR_PAD 2
#define eBSPLINE_PAD 4

typedef struct {
    int N_Site_Types;		// no. of site types
    tW_word *Site_Types;	// list of names for site types
    tW_word efile_name;
    tW_word vfile_name;					// RJS 10/07/2020
    tW_word calc_mode;
    bool bCalcMode_Whole;
    bool bCalcMode_SecondHalf;

    bool bfTraj, bfTop, bfEfile, bfVfile, bfSiteTypes, bfPairs, bfBonds, bfAngles, bfDihedrals,bfRcond, bfCalcmode, bfVolume;		//added bfVfile, bfVolume RJS 10/07/2020

    int n_particles;
    int n_frames;
    int n_total_bins;
    tW_gmx_trxframe *fr;
    tW_gmx_topology *top;

    int *nexcla;
    int **excl_list;

    int n_pair_types;
    tW_word **pair_atom_types;
    int **pair_types;
    double *dr;
    double *r_min;
    double *r_max;
    int *nb_basis;
    int *n_nb_pad;
    bool *nb_start_pad;
    bool *nb_end_pad;
    int *n_nb_bins;
    int *nb_i0s;

    int n_bond_types;
    tW_word **bond_atom_types;
    int *bond_types;
    double *db;
    double *b_min;
    double *b_max;
    int *bond_basis;
    int *n_bond_bins;
    int n_bonds;
    int *bond_list;
    int *bond_i0s;

    int n_angle_types;
    tW_word **angle_atom_types;
    int *angle_types;
    double *da;
    double *a_min;
    double *a_max;
    int *angle_basis;
    int *n_angle_bins; 
    int n_angles;
    int *angle_list;
    int *angle_i0s;

    int n_dihedral_types;
    tW_word **dihedral_atom_types;
    int *dihedral_types;
    double *dd;
    double *d_min;
    double *d_max;
    int *dihedral_basis;
    int *n_dihedral_bins;
    int n_dihedrals;
    int *dihedral_list;
    int *dihedral_i0s;

    /* RJS 10/07/2020 */
    int *n_vol_types;
    double *dv;
    double *v_min;
    double *v_max;
    int *volume_basis;
    int *n_volume_bins;
    int *volume_i0s;
    /* RJS 10/02/2020 */

    double **A; // The A matrix
    double *A_lin; // Linear A Matrx
    double *Ev; // the energy vector
    double *E_times;
    double *V; //The volume vector		// RJS 10/07/2020
    double *V_times; 				// RJS 10/07/2020
    
    double **A_trim;
    double *A_trim_lin;
    double *A_trim_lin_colmajor;
    bool *trim_cols;
    int n_zero_cols;
    int n_total_bins_trim;
    double *phi;
    double rcond;    

} tW_esystem;

#ifdef __cplusplus
}
#endif

#endif
