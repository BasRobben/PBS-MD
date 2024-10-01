#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

// This function calculates all forces acting on the particles (bonded and non-bonded).
// It initializes the forces array, then calculates bond-stretch, angle-bend, dihedral-torsion,
// and non-bonded forces. The total potential energy is returned.
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    // Initialize the forces to zero for all particles
    for (size_t i = 0; i < num_part; i++)
        f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    // Calculate the forces and accumulate the potential energy from each type of interaction
    double Epot = calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_angle(p_parameters, p_vectors);
    Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);

    return Epot;
}

// This function calculates bond-stretch forces based on the current positions of the bonded particles.
// It applies the minimum image convention to calculate the distance between bonded pairs and then
// computes the force and potential energy due to the bond interaction.
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi = {0.0, 0.0, 0.0};

    double kb = p_parameters->kb;
    double r0 = p_parameters->r0;

    // Loop through each bond and calculate the forces
    for (size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        // Calculate distance and potential energy
        double r_ij = sqrt(rij.x * rij.x + rij.y * rij.y + rij.z * rij.z);
        Epot += 0.5 * kb * (r_ij - r0) * (r_ij - r0);

        // Calculate the force on particles i and j
        double inv_r_ij = 1.0 / r_ij;
        fi.x = kb * (r_ij - r0) * rij.x * inv_r_ij;
        fi.y = kb * (r_ij - r0) * rij.y * inv_r_ij;
        fi.z = kb * (r_ij - r0) * rij.z * inv_r_ij;

        // Update forces on particles i and j
        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;
    }

    return Epot;  // Return the potential energy due to bond-stretch interactions
}

// This function calculates angle-bend forces based on the current positions of the angle-defined particles.
// It uses the minimum image convention and computes forces due to angle interactions.
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Angle *angles = p_vectors->angles;
    size_t num_angles = p_vectors->num_angles;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0};

    // Loop through each angle and calculate the forces
    for (size_t q = 0; q < num_angles; ++q)
    {
        size_t i = angles[q].i;
        size_t j = angles[q].j;
        size_t k = angles[q].k;

        // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        rkj.x = r[k].x - r[j].x;
        rkj.x = rkj.x - L.x * floor(rkj.x / L.x + 0.5);
        rkj.y = r[k].y - r[j].y;
        rkj.y = rkj.y - L.y * floor(rkj.y / L.y + 0.5);
        rkj.z = r[k].z - r[j].z;
        rkj.z = rkj.z - L.z * floor(rkj.z / L.z + 0.5);

        // TODO: Provide the angle force calculation and assign forces to particles i, j, and k

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= (fi.x + fk.x);
        f[j].y -= (fi.y + fk.y);
        f[j].z -= (fi.z + fk.z);
        f[k].x += fk.x;
        f[k].y += fk.y;
        f[k].z += fk.z;
    }

    return Epot;  // Return the potential energy due to angle-bend interactions
}

// This function calculates dihedral-torsion forces based on the positions of four connected particles.
// It uses the minimum image convention and computes the forces resulting from the dihedral-torsion interaction.
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Dihedral *dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj, rkl;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0}, fl = {0.0, 0.0, 0.0};

    // Loop through each dihedral and calculate the forces
    for (size_t q = 0; q < num_dihedrals; ++q)
    {
        size_t i = dihedrals[q].i;
        size_t j = dihedrals[q].j;
        size_t k = dihedrals[q].k;
        size_t l = dihedrals[q].l;

        // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        rkj.x = r[k].x - r[j].x;
        rkj.x = rkj.x - L.x * floor(rkj.x / L.x + 0.5);
        rkj.y = r[k].y - r[j].y;
        rkj.y = rkj.y - L.y * floor(rkj.y / L.y + 0.5);
        rkj.z = r[k].z - r[j].z;
        rkj.z = rkj.z - L.z * floor(rkj.z / L.z + 0.5);

        rkl.x = r[l].x - r[k].x;
        rkl.x = rkl.x - L.x * floor(rkl.x / L.x + 0.5);
        rkl.y = r[l].y - r[k].y;
        rkl.y = rkl.y - L.y * floor(rkl.y / L.y + 0.5);
        rkl.z = r[l].z - r[k].z;
        rkl.z = rkl.z - L.z * floor(rkl.z / L.z + 0.5);

        // TODO: Provide the dihedral-torsion force calculation and assign forces to particles i, j, k, and l

    }

    return Epot;  // Return the potential energy due to dihedral-torsion interactions
}

// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the Lennard-Jones potential.
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D df;
    double r_cutsq, sigmasq, sr2, sr6, sr12, fr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;
    sigmasq = p_parameters->sigma * p_parameters->sigma;
    double epsilon = p_parameters->epsilon;

    double Epot = 0.0, Epot_cutoff;
    sr2 = sigmasq / r_cutsq;
    sr6 = sr2 * sr2 * sr2;
    sr12 = sr6 * sr6;
    Epot_cutoff = sr12 - sr6;

    // Loop through the neighbor list and calculate the forces for each particle pair
    for (size_t k = 0; k < num_nbrs; k++)
    {
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;

        // Compute forces if the distance is smaller than the cutoff distance
        if (rij.sq < r_cutsq)
        {
            sr2 = sigmasq / rij.sq;
            sr6 = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

            // Calculate the potential energy
            Epot += 4.0 * epsilon * (sr12 - sr6 - Epot_cutoff);

            // Compute the force and apply it to both particles
            fr = 24.0 * epsilon * (2.0 * sr12 - sr6) / rij.sq;  // Force divided by distance
            df.x = fr * rij.x;
            df.y = fr * rij.y;
            df.z = fr * rij.z;

            f[i].x += df.x;
            f[i].y += df.y;
            f[i].z += df.z;
            f[j].x -= df.x;
            f[j].y -= df.y;
            f[j].z -= df.z;
        }
    }

    return Epot;  // Return the potential energy due to non-bonded interactions
}
