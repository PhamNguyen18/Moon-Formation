/*
 * Moon accretion from a disk a planetesimals. 
 * Modified from the planetesimal disk migration
 * example
 */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

double roche_radius(double r_earth);
double particle_r(double m, double m_earth, double roche_radius);
void heartbeat(struct reb_simulation* r);
double disk_angular_momentum(struct reb_simulation* const r, double ang_sys);
double E0;

int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c);



int main(int argc, char* argv[]){
    
    bool do_integrate = true;
    struct reb_simulation* r = reb_create_simulation();
    char* filename = "/directory_to_save_simulation_archive.bin";
    
    // Simulation Setup
    //r->integrator    = REB_INTEGRATOR_MERCURIUS;
    r->integrator    = REB_INTEGRATOR_IAS15;
    r->heartbeat    = heartbeat;
    r-> G = 6.67430e-11; // MKS CODATA value
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge_pass_through;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Constants
    double m_earth = 5.9722e24; // Kg
    double r_earth = 6.371e6; // m
    double m_lunar = m_earth * 0.0123;
    double aroche = roche_radius(r_earth); // Roche radius
    double t_kepler = 25200; // Keplerian orbital period at roche (~7 h)
    double J_EM = 3.5e34; // Current Moon-Earth angular momentum
    
    // Boundaries
    // Note that particles that fall onto Earth will merge with it!
    r->boundary    = REB_BOUNDARY_OPEN;
    const double boxsize = 20*r_earth; // in 20 Earth radii
    reb_configure_box(r,boxsize,2,2,2);
    
    srand(12);
    r->dt = 0.002*t_kepler;

    // Earth 
    struct reb_particle earth = {0};
    earth.m         = m_earth;
    earth.r        = r_earth;        // m
    reb_add(r, earth);

    // Planetesimal disk parameters
    int N_planetesimals = 1500; // 1000
    double powerlaw = -2; // -1
    double amin = 0.35*aroche, amax = 1.5*aroche;
    double m_tot = 0.0;
    
    // Generate Planetesimal Disk
    while(r->N<N_planetesimals + 1){
        struct reb_particle pt = {0};
        double a = 0;
        double dist = 0.0;
        // Adjust to change initial disk mass
        double m = reb_random_powerlaw(m_lunar*5.7e-5, m_lunar*5e-2, -1.5); 
        m_tot += m;
        bool param = true;
        while (param) {
          a = reb_random_powerlaw(amin, amax, powerlaw);
          double e    = reb_random_rayleigh(0.24);
          double inc  = reb_random_rayleigh(0.24);
          double Omega = reb_random_uniform(0,2.*M_PI);
          double apsis = reb_random_uniform(0,2.*M_PI);
          double phi     = reb_random_uniform(0,2.*M_PI);
          pt = reb_tools_orbit_to_particle(r->G, earth, m, a, e, inc, Omega, apsis, phi);
          dist = pt.x*pt.x + pt.y*pt.y + pt.z*pt.z;

          // Make sure no particles are initialized inside the Earth!
          if (a > amin && a < amax) {
            if (dist > amin*amin) {
              param = false;
            }
          }
        }
        pt.r         = particle_r(m, m_earth, aroche);
        reb_add(r, pt);
    }
    // All particles experience gravity
    r->N_active = r->N;

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    printf("Disk mass total: %f \n", m_tot/m_lunar);
    printf("Disk angular momentum: %f \n", disk_angular_momentum(r, J_EM));

    if (do_integrate) {
      reb_simulationarchive_automate_interval(r, filename, 10*t_kepler);
      reb_integrate(r, 1550*t_kepler);
      reb_free_simulation(r);
    }

}

double disk_angular_momentum(struct reb_simulation* const r, double ang_sys) {
  struct reb_vec3d data = reb_tools_angular_momentum(r);
  double sq_data = data.x*data.x + data.y*data.y + data.z*data.z;
  double ang_tot = pow(sq_data, 1/2.0)/ ang_sys;
  return ang_tot;
}

double roche_radius(double r_earth) {
  return 2.456*pow((5.5/3.34), 1.0/3.0) * r_earth;
}

double particle_r(double m, double m_earth, double roche_radius) {
  // Ida 1997 implementation of particle radius based on mass
  return (1.0/2.456) * pow(m/m_earth, 1.0/3.0) * roche_radius;
}

int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c){
    // This function passes the collision to the default merging routine. 
    // If a merger occured, that routine will return a value other than 0.
    // This function then outputs some information about the merger.
    int result = reb_collision_resolve_merge(r, c);
    if (result!=0){
       printf("A merger occured! Particles involved: %d, %d.\n",c.p1,c.p2);
    }
    return result;
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10*25200)){
      double E = reb_tools_energy(r);
      double relE = fabs((E-E0)/E0);
      double t_kepler = 7*3600;
      
      //get orbital elements
      struct reb_particle p = r->particles[2];
      struct reb_particle star = r->particles[0];
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G,p,star);
      printf("t=%f, a2=%f,dE=%e,N=%d\n",r->t/t_kepler, o.a,relE,r->N);
      
    }
}
