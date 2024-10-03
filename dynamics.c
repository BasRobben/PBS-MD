#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// This function updates particle positions using their velocities.
// The positions are advanced by one full time step (dt), and displacement vectors
// (dr) for one time step are updated. The displacement since the last neighbor 
// list creation (stored in p_nbrlist->dr) is also updated for each particle.
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr_loc;
    struct Vec3D *r = p_vectors->r;     // Particle positions
    struct Vec3D *dr = p_vectors->dr;   // Displacement in one timestep
    struct Vec3D *v = p_vectors->v;     // Particle velocities
    struct DeltaR *dr_nbrlist = p_nbrlist->dr;  // Displacement since last neighbor list creation
    size_t num_part = p_parameters->num_part;
    double dt = p_parameters->dt;

    // Loop over all particles to update their positions
    for (size_t i = 0; i < num_part; i++)
    {
        dr_loc.x = v[i].x * dt;  // Compute displacement in x-direction for one timestep
        dr_loc.y = v[i].y * dt;  // Compute displacement in y-direction for one timestep
        dr_loc.z = v[i].z * dt;  // Compute displacement in z-direction for one timestep

        dr[i] = dr_loc;          // Store the displacement for this timestep
        r[i].x += dr_loc.x;      // Update position in x-direction
        r[i].y += dr_loc.y;      // Update position in y-direction
        r[i].z += dr_loc.z;      // Update position in z-direction

        // Update the displacement since last neighbor list creation
        dr_nbrlist[i].x += dr_loc.x;
        dr_nbrlist[i].y += dr_loc.y;
        dr_nbrlist[i].z += dr_loc.z;
        dr_nbrlist[i].sq = (dr_nbrlist[i].x) * (dr_nbrlist[i].x) + 
                           (dr_nbrlist[i].y) * (dr_nbrlist[i].y) + 
                           (dr_nbrlist[i].z) * (dr_nbrlist[i].z);  // Square of total displacement since last neighbor list creation
    }
}

// This function updates particle velocities by half a time step using the current forces.
// The updated velocities are used in the velocity-Verlet integration scheme.
// The function also calculates and returns the kinetic energy of the system.
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    double Ekin = 0.0;  // Initialize kinetic energy
    const double factor_m = 0.5 / p_parameters->mass_m * p_parameters->dt;  // Factor for velocity update
    const double factor_e = 0.5 / p_parameters->mass_e * p_parameters->dt;  // Factor for velocity update
    struct Vec3D *v = p_vectors->v;  // Particle velocities
    struct Vec3D *f = p_vectors->f;  // Forces acting on particles

    // Loop over all particles and update their velocities
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        // Check particle type and update velocity accordingly
        if (p_vectors -> type[i] == 1) // Methane
        {
            v[i].x += factor_m * f[i].x;  // Update velocity in x-direction
            v[i].y += factor_m * f[i].y;  // Update velocity in y-direction
            v[i].z += factor_m * f[i].z;  // Update velocity in z-direction
            Ekin += p_parameters->mass_m * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z); // Accumulate kinetic energy
        }
        else if (p_vectors -> type[i] == 0) // Ethane
        {
            v[i].x += factor_e * f[i].x;  // Update velocity in x-direction
            v[i].y += factor_e * f[i].y;  // Update velocity in y-direction
            v[i].z += factor_e * f[i].z;  // Update velocity in z-direction
            Ekin += p_parameters->mass_e * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z); // Accumulate kinetic energy
        }

    }

    // Final kinetic energy calculation
    Ekin = 0.5 * Ekin;
    return Ekin;  // Return the system's kinetic energy
}

// This function applies periodic boundary conditions to ensure particles stay inside the simulation box.
// If a particle moves beyond the box, it is wrapped around to the opposite side.
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D invL;  // Inverse of the box size
    struct Vec3D *r = p_vectors->r;  // Particle positions
    struct Vec3D L = p_parameters->L;  // Box dimensions
    size_t num_part = p_parameters->num_part;  // Number of particles

    invL.x = 1.0 / L.x;
    invL.y = 1.0 / L.y;
    invL.z = 1.0 / L.z;

    // Loop over all particles and apply periodic boundary conditions
    for (size_t i = 0; i < num_part; i++)
    {
        r[i].x -= L.x * floor(r[i].x * invL.x);  // Apply periodic boundary in x-direction
        r[i].y -= L.y * floor(r[i].y * invL.y);  // Apply periodic boundary in y-direction
        r[i].z -= L.z * floor(r[i].z * invL.z);  // Apply periodic boundary in z-direction
    }
}

// This function applies a thermostat to maintain the system's temperature.
// The Berendsen thermostat is used to rescale the velocities of particles based on the system's kinetic energy.
double thermostat(struct Parameters *p_parameters, struct Vectors *p_vectors, double Ekin)
{
    struct Vec3D *v = p_vectors->v;             // Particle velocities
    size_t num_part = p_parameters->num_part;   // Number of particles
    double TmeasT0;                             // Measured temperature divided by T0
    double lambda;                              // Scaling factor
    double kT = p_parameters->kT;               // Thermal energy

    double Nfree = 3*(num_part - 1);            // Number of degrees of freedom
    double tau = 0.001;                          // Relaxation time for the thermostat
    Ekin = 0;                                   // Reset kinetic energy

    for (size_t i = 0; i < num_part; i++)
    {
        if (p_vectors -> type[i] == 1) // Methane
            Ekin += p_parameters->mass_m * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);  // Calculate kinetic energy
        else if (p_vectors -> type[i] == 0) // Ethane
            Ekin += p_parameters->mass_e * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);  // Calculate kinetic energy
    }

    TmeasT0 = Ekin / (Nfree * kT);                // Measured temperature divided by T0
    // lambda = sqrt(1.0 + p_parameters->dt / tau * (1 / TmeasT0 - 1.0));  // Scaling factor

    // for (size_t i = 0; i < num_part; i++)
    // {
    //     v[i].x *= lambda;  // Rescale velocity in x-direction
    //     v[i].y *= lambda;  // Rescale velocity in y-direction
    //     v[i].z *= lambda;  // Rescale velocity in z-direction
    // }
    
    return TmeasT0;  // Return the temperature ratio
}