/**
 * @file 	boundary.c
 * @brief 	Implementation of all boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * 
 * 
 * @section LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "rebound.h"
#include "boundary.h"
#include "tree.h"

// New function 
struct reb_vec3d angular_momentum(const struct reb_simulation* const r, int i, const double ang_scale){
	struct reb_particle pi = r->particles[i];
	//double Lx = 0, Ly = 0, Lz = 0;
	//double ang_p = 0;
	struct reb_vec3d L = {0};	

    L.x = pi.m*(pi.y*pi.vz - pi.z*pi.vy);
    L.y = pi.m*(pi.z*pi.vx - pi.x*pi.vz);
    L.z = pi.m*(pi.x*pi.vy - pi.y*pi.vx);
	//ang_p = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)/ang_scale;
	return L;
}

// Write angular momentum of particles that leave box to file
void write_to_file(struct reb_vec3d L) {
	FILE *fp;

    fp = fopen("/home/pham/Desktop/ang_out.txt", "a");
	fprintf(fp, "%f, %f, %f\n", L.x, L.y, L.z);
    fclose(fp);
}

// Open file to output data? 
void reb_boundary_check(struct reb_simulation* const r){
	struct reb_particle* const particles = r->particles;
	int N = r->N;
	const struct reb_vec3d boxsize = r->boxsize;
	switch(r->boundary){
		case REB_BOUNDARY_OPEN:
			for (int i=0;i<N;i++){ // run through loop backwards so we don't have to recheck same index after removing
				int removep = 0;
				if(particles[i].x>boxsize.x/2.){
					removep = 1;
				}
				if(particles[i].x<-boxsize.x/2.){
					removep = 1;
				}
				if(particles[i].y>boxsize.y/2.){
					removep = 1;
				}
				if(particles[i].y<-boxsize.y/2.){
					removep = 1;
				}
				if(particles[i].z>boxsize.z/2.){
					removep = 1;
				}
				if(particles[i].z<-boxsize.z/2.){
					removep = 1;
				}
				if (removep==1){
					struct reb_vec3d ang_mom = angular_momentum(r, i, 3.5e34);
                    if(r->track_energy_offset){
						write_to_file(ang_mom);
						printf("Particle left boundary!");
                        double Ei = reb_tools_energy(r);
                        reb_remove(r, i,1);
                        r->energy_offset += Ei - reb_tools_energy(r);
                    } else {
						write_to_file(ang_mom);
						printf("Particle left boundary!");
                    	reb_remove(r, i,0); // keepSorted=0 by default in C version
                    }
                    if (r->tree_root==NULL){
                        i--; // need to recheck the particle that replaced the removed one
                        N--; // This is the local N
                    }else{
                        // particle just marked, will be removed later
                        r->tree_needs_update= 1;
                    }
				}
			}
			break;
		case REB_BOUNDARY_SHEAR:
		{
			// The offset of ghostcell is time dependent.
			const double OMEGA = r->ri_sei.OMEGA;
			const double offsetp1 = -fmod(-1.5*OMEGA*boxsize.x*r->t+boxsize.y/2.,boxsize.y)-boxsize.y/2.; 
			const double offsetm1 = -fmod( 1.5*OMEGA*boxsize.x*r->t-boxsize.y/2.,boxsize.y)+boxsize.y/2.; 
			struct reb_particle* const particles = r->particles;
#pragma omp parallel for schedule(guided)
			for (int i=0;i<N;i++){
				// Radial
				while(particles[i].x>boxsize.x/2.){
					particles[i].x -= boxsize.x;
					particles[i].y += offsetp1;
					particles[i].vy += 3./2.*OMEGA*boxsize.x;
				}
				while(particles[i].x<-boxsize.x/2.){
					particles[i].x += boxsize.x;
					particles[i].y += offsetm1;
					particles[i].vy -= 3./2.*OMEGA*boxsize.x;
				}
				// Azimuthal
				while(particles[i].y>boxsize.y/2.){
					particles[i].y -= boxsize.y;
				}
				while(particles[i].y<-boxsize.y/2.){
					particles[i].y += boxsize.y;
				}
				// Vertical (there should be no boundary, but periodic makes life easier)
				while(particles[i].z>boxsize.z/2.){
					particles[i].z -= boxsize.z;
				}
				while(particles[i].z<-boxsize.z/2.){
					particles[i].z += boxsize.z;
				}
			}
		}
		break;
		case REB_BOUNDARY_PERIODIC:
#pragma omp parallel for schedule(guided)
			for (int i=0;i<N;i++){
				while(particles[i].x>boxsize.x/2.){
					particles[i].x -= boxsize.x;
				}
				while(particles[i].x<-boxsize.x/2.){
					particles[i].x += boxsize.x;
				}
				while(particles[i].y>boxsize.y/2.){
					particles[i].y -= boxsize.y;
				}
				while(particles[i].y<-boxsize.y/2.){
					particles[i].y += boxsize.y;
				}
				while(particles[i].z>boxsize.z/2.){
					particles[i].z -= boxsize.z;
				}
				while(particles[i].z<-boxsize.z/2.){
					particles[i].z += boxsize.z;
				}
			}
		break;
		default:
		break;
	}
}

const static struct reb_ghostbox nan_ghostbox = {.shiftx = 0, .shifty = 0, .shiftz = 0, .shiftvx = 0, .shiftvy = 0, .shiftvz = 0};

struct reb_ghostbox reb_boundary_get_ghostbox(struct reb_simulation* const r, int i, int j, int k){
	switch(r->boundary){
		case REB_BOUNDARY_OPEN:
		{
			struct reb_ghostbox gb;
			gb.shiftx = r->boxsize.x*(double)i;
			gb.shifty = r->boxsize.y*(double)j;
			gb.shiftz = r->boxsize.z*(double)k;
			gb.shiftvx = 0;
			gb.shiftvy = 0;
			gb.shiftvz = 0;
			return gb;
		}
		case REB_BOUNDARY_SHEAR:
		{
			const double OMEGA = r->ri_sei.OMEGA;
			struct reb_ghostbox gb;
			// Ghostboxes habe a finite velocity.
			gb.shiftvx = 0.;
			gb.shiftvy = -1.5*(double)i*OMEGA*r->boxsize.x;
			gb.shiftvz = 0.;
			// The shift in the y direction is time dependent. 
			double shift;
			if (i==0){
				shift = -fmod(gb.shiftvy*r->t,r->boxsize.y); 
			}else{
				if (i>0){
					shift = -fmod(gb.shiftvy*r->t-r->boxsize.y/2.,r->boxsize.y)-r->boxsize.y/2.; 
				}else{
					shift = -fmod(gb.shiftvy*r->t+r->boxsize.y/2.,r->boxsize.y)+r->boxsize.y/2.; 
				}	
			}
			gb.shiftx = r->boxsize.x*(double)i;
			gb.shifty = r->boxsize.y*(double)j-shift;
			gb.shiftz = r->boxsize.z*(double)k;
			return gb;
		}
		case REB_BOUNDARY_PERIODIC:
		{
			struct reb_ghostbox gb;
			gb.shiftx = r->boxsize.x*(double)i;
			gb.shifty = r->boxsize.y*(double)j;
			gb.shiftz = r->boxsize.z*(double)k;
			gb.shiftvx = 0;
			gb.shiftvy = 0;
			gb.shiftvz = 0;
			return gb;
		}
		default:
			return nan_ghostbox;
	}
}

int reb_boundary_particle_is_in_box(const struct reb_simulation* const r, struct reb_particle p){
	switch(r->boundary){
		case REB_BOUNDARY_OPEN:
		case REB_BOUNDARY_SHEAR:
			if(p.x>r->boxsize.x/2.){
				return 0;
			}
			if(p.x<-r->boxsize.x/2.){
				return 0;
			}
			if(p.y>r->boxsize.y/2.){
				return 0;
			}
			if(p.y<-r->boxsize.y/2.){
				return 0;
			}
			if(p.z>r->boxsize.z/2.){
				return 0;
			}
			if(p.z<-r->boxsize.z/2.){
				return 0;
			}
			return 1;
		default:
			return 1;
	}
}


