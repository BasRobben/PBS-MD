#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters)
{
  double kT = BOLTZMANN * 298.0;    // thermal energy at 298K [J]

  // // Methane parameters (Task 1 & 2)
  // // physical units
  // double epsilon = 148.0 * BOLTZMANN;  // LJ interaction strength [J]
  // double sigma = 3.73e-10;             // LJ particle diameter [m]
  // double mass = 2.66e-26;              // mass of a CH4 particle [kg]
  // double density = 0.657;              // density of CH4 at 298K [kg/m3]
  // // reduced units
  // p_parameters->epsilon = 1.0;         // Set epsilon to 1 in reduced units
  // p_parameters->sigma = 1.0;           // Set sigma to 1 in reduced units
  // p_parameters->mass = 1.0;            // Set mass to 1 in reduced units
  // p_parameters->kT = kT / epsilon;     // Thermal energy in reduced units

  // Ethane parameters (Task 3)
  // physical units
  double epsilon = 98.0 * BOLTZMANN;          // LJ interaction strength [J]
  double sigma = 3.75e-10;                    // LJ particle diameter [m]
  double mass = 4.99e-26 / 2;                 // mass of a ethane particle [kg]
  double density = 1;                     // density of ethane at 298K [kg/m3]
  double kb = 2.5e5 * BOLTZMANN * (1e20);     // Bond force constant [J/m2]
  double r0 = 1.54e-10;                       // Equilibrium bond length [m]

  // reduced units
  p_parameters->epsilon = 1.0;      // Set epsilon to 1 in reduced units
  p_parameters->sigma = 1.0;        // Set sigma to 1 in reduced units
  p_parameters->mass = 1.0;         // Set mass to 1 in reduced units
  p_parameters->kT = kT / epsilon;  // Thermal energy in reduced units
  p_parameters->kb = kb * sigma * sigma / epsilon;  // Bond force constant in reduced units
  p_parameters->r0 = r0 / sigma;    // Equilibrium bond length in reduced units

  // Adjust the box size to achieve the desired density
  int num_part = 100;                          // number of particles
  double volume = (num_part * mass) / density;  // volume of the box
  double box_length = pow(volume, 1.0 / 3.0);   // length of the box
  double L_reduced = box_length / sigma;        // box length in reduced units

  // The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = num_part;                        //number of particles
  p_parameters->num_dt_steps = 20000;                        //number of time steps
  p_parameters->exclude_12_nb = 1;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 1;                          // 1-3 connected atoms exluded from non-bonded interactions
  p_parameters->dt = 0.01;                                  //integration time step
  p_parameters->L = (struct Vec3D){L_reduced, L_reduced, L_reduced}; //box size adjusted for density
  p_parameters->r_cut = 2.68;                               //cut-off distance used for neigbor list
  p_parameters->r_shell = 0.5;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 500;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file

  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
