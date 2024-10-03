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
  double epsilon_m = 148.0 * BOLTZMANN;  // LJ interaction strength [J]
  double sigma_m = 3.73e-10;             // LJ particle diameter [m]
  double mass_m = 2.66e-26;              // mass of a CH4 particle [kg]

  // Ethane parameters (Task 3)
  // physical units
  double epsilon_e = 98.0 * BOLTZMANN;          // LJ interaction strength [J]
  double sigma_e = 3.75e-10;                    // LJ particle diameter [m]
  double mass_e = 2.50e-26;                     // mass of a CH3 unit [kg]

  double kb = 2.5e5 * BOLTZMANN * (1e20);     // Bond force constant [J/m2]
  double r0 = 1.54e-10;                       // Equilibrium bond length [m]
  double density = 200;                       // density of the system [kg/m3]

  // reduced units
  p_parameters->epsilon_m = 1.0;                // Set epsilon to 1 in reduced units
  p_parameters->sigma_m = 1.0;                  // Set sigma to 1 in reduced units
  p_parameters->mass_m = 1.0;                   // Set mass to 1 in reduced units

  p_parameters->epsilon_e = epsilon_e / epsilon_m;  // LJ interaction strength in reduced units
  p_parameters->sigma_e = sigma_e / sigma_m;        // LJ particle diameter in reduced units
  p_parameters->mass_e = mass_e / mass_m;           // mass of a CH3 unit in reduced units

  p_parameters->kT = kT / epsilon_m;            // Thermal energy in reduced units
  p_parameters->kb = kb * sigma_m * sigma_m / epsilon_m;  // Bond spring constant in reduced units
  p_parameters->r0 = r0 / sigma_m;              // Equilibrium bond length in reduced units

  // Adjust the box size to achieve the desired density
  int num_part = 300;                          // number of particles

  double mole_fraction_methane = 0.5;          // Mole fraction of CH4

  // Calculate the number of CH4 and CH3 particles
  size_t num_methane = (size_t)((mole_fraction_methane * num_part)/(2.0 - mole_fraction_methane));
  size_t num_ethane = (size_t)((num_part - num_methane));

  // Check if the total number of particles matches the specified number
  if (num_methane + num_ethane != num_part) {
      fprintf(stderr, "Error: Total number of particles does not match specified number.\n");
      exit(EXIT_FAILURE);
  }

  // Assign the number of CH4 and CH3 particles to the parameters struct
  p_parameters -> num_methane = num_methane;
  p_parameters -> num_ethane = num_ethane;

  // Calculate box dimensions based on the number of particles, average mass and density
  double volume = (num_part * (mass_m + 2.0 * mass_e)) / (3.0 * density);   // volume of the box
  double box_length = pow(volume, 1.0 / 3.0);                               // length of the box
  double L_reduced = box_length / sigma_m;                                  // box length in reduced units

  // The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = num_part;                        //number of particles
  p_parameters->num_dt_steps = 20000;                        //number of time steps
  p_parameters->exclude_12_nb = 1;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 1;                          // 1-3 connected atoms exluded from non-bonded interactions
  p_parameters->dt = 0.000001;                                  //integration time step
  p_parameters->L = (struct Vec3D){L_reduced, L_reduced, L_reduced}; //box size adjusted for density
  p_parameters->r_cut = 10e-10 / sigma_m;                   //cut-off distance used for neigbor list
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