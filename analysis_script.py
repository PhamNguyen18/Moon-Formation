"""
Author: Pham Nguyen
Last Updated: 8-11-20
"""
import rebound
import argparse
import numpy as np
import pandas as pd
from IPython.display import display, clear_output
from matplotlib.patches import Ellipse, Circle
import matplotlib.animation as manimation
import matplotlib.pyplot as plt


class ReboundAnalysis:
    """
    Class to perform various data analysis of rebound simulation data. 
    Includes visualization tools to compare to Ida 1997.
    ------------------------------------------------------------------
    Params: file        Filename/ pathway to simulation archive
            time        Give time slice to analyze data
            t_kepler    Keplerian orbital period at Roche limit 
            dt          Time step for making movies
            r_earth     Radius of Earth (or primary body in general)
            m_earth     Mass of Earth (or primary body in general)
            m_lunar     Mass of Moon (or another primary moon)
            J_EM        Earth-Moon angular momentum
            rho         Density of Earth (or primary body)
            rho0        Density of Moon (or primary moon)
    """

    G_const = 6.67e-11  # MKS

    def __init__(self, file, time, t_kepler=25200,
                 dt=252, r_earth=6.371e6, m_earth=5.9722e24,
                 m_lunar=7.345806e+22, J_EM=3.5e34, rho=5.5,
                 rho0=3.3):

        self.file = file
        self.time = time
        self.sim = rebound.SimulationArchive(file)
        self.sim_time = self.sim.getSimulation(time)
        self.t_kepler = t_kepler if t_kepler is not None else t_kepler
        self.dt = dt if dt is not None else dt 

        # Planetary constants
        self.r_earth = r_earth if r_earth is not None else r_earth  # m
        self.m_earth = m_earth if m_earth is not None else m_earth  # Kg
        self.m_lunar = m_lunar if m_lunar is not None else m_lunar
        self.J_EM = J_EM if J_EM is not None else J_EM  # MKS 
        self.rho = rho if rho is not None else rho  # g 
        self.rho0 = rho0 if rho0 is not None else rho0  # g 
        self.roche = 2.456*(self.rho/self.rho0)**(1/3)*self.r_earth


    @staticmethod
    def mag(x, y, z):
        # General method to calculate the magnitude of a 3D vector.
        return np.sqrt(x**2 + y**2 + z**2)

    @staticmethod
    def particle_radius(m, m_earth, roche_radius):
        # Calculate particle radius from mass as in Ida 1997
        return (1/2.456) * (m/m_earth)**(1/3.0) * roche_radius

    def roche_radius(self):
        # Calculates the roche radius for a given density and object radius.
        self.roche = 2.456*(self.rho/self.rho0)**(1/3)*self.r_earth
        return None

    def angular_momentum(self, p, flag="particle", mi=1):
        """
        Calculates the angular momentum of a
        single particle or many particles.

        Scale is used to output the ang. momentum
        in different units. By default it is set
        to the current Earth-Moon ang. momentum.

        Parameters
        ----------

        particle : object
            Either a single particle object or
            a simulation class
        flag : str
            Used to calculate the angular momentum
            of a single particle or a larger collection
            of particles.

            disk : Angular momentum of the entire disk
                   excluding the Earth.
            bound : Angular momentum of all bound particles
                    outside the orbit of the largest moon
                    (includes the moon as well)
            p : Angular momentum of a single particle

        Returns
        -------
        Angular momentum of particle or collection of particles
        """

        Lx = 0
        Ly = 0
        Lz = 0
        scale = self.J_EM

        if flag == "disk":
            # Calculate angular momentum of ALL disk particles. Starting at index 1 excludes the Earth!
            for i in range(1, p.N):
                pi = p.particles[i]
                Lx += pi.m*(pi.y*pi.vz - pi.z*pi.vy) 
                Ly += pi.m*(pi.z*pi.vx - pi.x*pi.vz) 
                Lz += pi.m*(pi.x*pi.vy - pi.y*pi.vx) 
        elif flag == "bound":
            # Calculate angular momentum of largest mooon and all bound debris outside its orbit
            for i in range(mi, p.N):
                pi = p.particles[i]
                Lx += pi.m*(pi.y*pi.vz - pi.z*pi.vy) 
                Ly += pi.m*(pi.z*pi.vx - pi.x*pi.vz) 
                Lz += pi.m*(pi.x*pi.vy - pi.y*pi.vx) 
        elif flag == "particle":
            # Calculate angular momentum of a single particle
            Lx = p.m*(p.y*p.vz - p.z*p.vy) 
            Ly = p.m*(p.z*p.vx - p.x*p.vz) 
            Lz = p.m*(p.x*p.vy - p.y*p.vx) 

        return self.mag(Lx, Ly, Lz)/scale

    def data_summary(self, G_const):
        """
        Method to return simulation data related to moon and disk
        mass and specific angular momentum for comparison to Ida. 
        These values are calculated for the largest moon,
        bound particles with orbits outside the largest moon,
        and ejected material. 
        """

        sim_init_mass = 0
        sim_init = self.sim.tmin
        sim_end = self.sim.tmax
        sim_init_ang = self.angular_momentum(sim_init)
        # sim_end_ang = self.angular_momentum(sim_end)

        earth_accreted_mass = (sim_end.particles[0].m - sim_init.particles[0].m) / self.m_lunar
        # const = np.sqrt(G_const * self.m_earth * self.roche) 

        # End of simulation properties
        #angular_momentum_end = []
        masses_end = []
        semi_maj_end = []
        orbs_end = sim_end.calculate_orbits()

        for i in range(1, sim_init.N):
            sim_init_mass += sim_init.particles[i].m / self.m_lunar

        for i in range(1, sim_end.N):
            masses_end.append(sim_end.particles[i].m / self.m_lunar)
            semi_maj_end.append(orbs_end[i-1].a/self.roche)

        sim_end_mass = np.sum(masses_end)
        mi = int(np.argmax(masses_end)) + 1  # index of largest mass at end of simulation

        print("Ang. Mom. Initial:", sim_init_ang)
        print("Ang. Mom. (largest moon): ", self.angular_momentum(sim_end.particles[mi], "particle"))
        print("Ang. Mom. (disk bound): ", self.angular_momentum(sim_end, "bound", mi))
        #print("Ang. Mom. (ejected): ", ang_ejected) 
        print("Ang. Mom. Earth", self.angular_momentum(sim_end.particles[0], "particle"))
        print("Mass (largest moon): ", sim_end.particles[mi].m/self.m_lunar)
        print("Mass (disk bound): ", np.sum(masses_end[mi-1:]))
        print("Mass (ejected): ", sim_init_mass - earth_accreted_mass - sim_end_mass) 

        return None

    def plot_orbs(self, times, title, name, show=True, dpi=150, format='png', save=False):
        """
        Plots the orbits of planetesimals at three different
        times in the r-z plane. 

        TO DO:
        Should be a better way to replace "name" and "dpi"
        Could use kwargs perhaps to pass in other plotting
        parameters.
        """

        t0, t1, t2 = times * self.t_kepler 
        sims = [self.sim.getSimulation(t0),
                self.sim.getSimulation(t1),
                self.sim.getSimulation(t2)]

        # Get z and r coordinates and radii of all particles for each time
        z_coords = []
        r_coords = []
        radii = []

        for i in range(3):
            z_temp = []
            r_temp = []
            radii_temp = []
            for j in range(1, sims[i].N):
                # Switch to geocentric coordinates
                x = (sims[i].particles[j].x - sims[i].particles[0].x)/self.roche
                y = (sims[i].particles[j].y - sims[i].particles[0].y)/self.roche
                z = (sims[i].particles[j].z - sims[i].particles[0].z)/self.roche
                z_temp.append(np.abs(z))
                r_temp.append(self.mag(x, y, z))
                radii_temp.append(sims[i].particles[j].r/self.roche)
            z_coords.append(z_temp)
            r_coords.append(r_temp)
            radii.append(radii_temp)

        # Grab masses from last time
        masses = []
        for i in range(1, sims[2].N):
            masses.append(sims[2].particles[i].m/self.m_lunar)

        # Find a way to set global parameters for all subplots
        textsize = 16 
        fig, ax = plt.subplots(3, sharex=True, sharey=True, figsize=(7,8), gridspec_kw={"hspace": 0})

        # Top figure
        ax[0].set_title(title + " N = {}, M_Disk = {} M_L".format(1500, 2.67), size=textsize)
        ax[0].axvline(0.35, ls="--", color="gray")
        ax[0].axvline(0.35, ls="--", color="gray")
        ax[0].tick_params(direction="in", labelsize=textsize)

        # Scaling to keep points circular on plot
        x0, y0 = ax[0].transAxes.transform((0, 0))
        x1, y1 = ax[0].transAxes.transform((1.25, 3.0))
        dx = x1 - x0
        dy = y1 - y0
        maxd = max(dx, dy)

        for r, z, rad in zip(r_coords[0], z_coords[0], radii[0]):
            width = rad * maxd / dx
            height = rad * maxd / dy
            ax[0].add_artist(Ellipse((r, z), width, height, alpha=0.7))

        # Middle figure
        ax[1].set_ylabel("z $(a_R)$", size=textsize)
        ax[1].axvline(0.35, ls="--", color="gray")
        ax[1].tick_params(direction="in", labelsize=textsize)
        ax[1].text(0.4, 0.6, r'$R_\oplus$', size=textsize)

        for r, z, rad in zip(r_coords[1], z_coords[1], radii[1]):
            width = rad * maxd / dx
            height = rad * maxd / dy
            ax[1].add_artist(Ellipse((r, z), width, height, alpha=0.7))

        # Bottom figure
        ax[2].set_ylim(-0.05, 1.2)
        ax[2].set_xlim(0, 3)
        ax[2].tick_params(direction="in", labelsize=textsize)
        ax[2].set_xlabel("r $(a_R)$", size=textsize)
        ax[2].axvline(0.35, ls="--", color="gray")
        maxd = max(dx, dy)

        for r, z, rad in zip(r_coords[2], z_coords[2], radii[2]):
            width = rad * maxd / dx
            height = rad * maxd / dy
            ax[2].add_artist(Ellipse((r, z), width, height, alpha=0.7))

        # Annotate any masses that are greater than 30% a lunar mass as in Ida
        # This should actually be 30% of the largest luanr seed that is formed
        for i, m in enumerate(masses):
            if m > 0.3*masses[1]:
                m_label =  str(round(masses[i], 3)) + r"$M_L$"
                ax[2].annotate(m_label, xy=(r_coords[2][i]*1.05, z_coords[2][i]*1.05), size=14)
            else:
                continue

        if show:
            plt.tight_layout()
            plt.show()
        if save: 
            plt.tight_layout()
            plt.savefig(name, dpi=dpi, format=format)

        return fig, ax 

    def plot_ang_mom(self, times):
        """
        Plots the angular momentum of the planetesimal
        disk as a function of time. 
        """

        # Momentum components
        massl = np.zeros_like(times)
        tot_mom = np.zeros_like(times)
        for i, t in enumerate(times):
            px, py, pz = 0, 0, 0
            mass_tot = 0
            sim_t = self.sim.getSimulation(t)
            for j in range(sim_t.N):
                p = sim_t.particles[j]
                px += p.m * p.vx
                py += p.m * p.vy
                pz += p.m * p.vz
                mass_tot += p.m
            
            momentum = np.sqrt(px**2 + py**2 + pz**2)
            tot_mom[i] = momentum
            massl[i] = mass_tot
            
        # Angular momentum
        tot_ang_mom = np.zeros_like(times)
        axl = np.zeros_like(times)
        ayl = np.zeros_like(times)
        azl = np.zeros_like(times)
         
        #earth_ang_mom = np.zeros_like(times)
        #disk_ang_mom = np.zeros_like(times)

        for i, t in enumerate(times):
            sim_t = self.sim.getSimulation(t)
            ax, ay, az = sim_t.calculate_angular_momentum()
            axl[i] = ax
            ayl[i] = ay
            azl[i] = az
            tot_ang_mom[i] = self.mag(ax, ay, az)
            #earth_ang_mom[i] = self.angular_momentum(sim_t.particles[0], "p")
            #disk_ang_mom[i] = self.angular_momentum(sim_t)

        t_to_kepler = times / self.t_kepler
        tot_ang_mom = tot_ang_mom #/ self.J_EM
        #earth_ang_mom = earth_ang_mom
        #disk_ang_mom = disk_ang_mom

        fig, ax = plt.subplots()
        ax.set_title("Angular Momentum of Two Earth Mass Particles")
        ax.set_yscale('log')
        ax.plot(t_to_kepler, tot_ang_mom)
        ax.set_ylim(1e34, 1e36)
        ax.set_xlabel("Time (T_Kep)")
        ax.set_ylabel("Angular Momentum")
        ax.legend()
        plt.savefig("angular_momentum.png", dpi=150)
        plt.show()

        return fig, ax 



    def make_video(self, title, comment, vid_title, flag="disk",dpi=150):
        """
        Creates a movie of the planetesimal disk in the x-y plane.
        Also possible to create movies in the z-r plane as well. 

        Parameters
        ----------
        title : str
            Video title for metadata
        comment : str
            Comment for metadata
        vid_title : str
            Title appears on plot for video

        Returns
        -------
        Saves a video file of the simulation

        This code has not be tested yet!
        """

        if flag == "disk": # Top down view of planetary disk in the X-Y plane
            FFMpegWriter = manimation.writers["ffmpeg"]
            metadata = dict(title=title, comment=comment)
            writer = FFMpegWriter(fps=15, metadata=metadata)
            fig = plt.figure(figsize=(5, 5))

            with writer.saving(fig, vid_title+".mp4", dpi=dpi):
                for t in range(0, int(self.sim.tmax), int(0.05*self.t_kepler)):
                    sim = self.sim.getSimulation(t)
                    ax = plt.subplot(1, 1, 1)
                    x_e = sim.particles[0].x/self.roche
                    y_e = sim.particles[0].y/self.roche
                    circ = plt.Circle(
                        (x_e, y_e), radius=sim.particles[0].r/self.roche, fill=False, ls="--")
                    ax.set_title("Time = " + str(round(t/self.t_kepler, 0)) + " T_Kepler")
                    ax.add_artist(circ)
                    for i in range(1, sim.N):
                        x = sim.particles[i].x/self.roche  # - sim.particles[0].x)/roche
                        y = sim.particles[i].y/self.roche  # - sim.particles[0].y)/roche
                        circ = plt.Circle(
                            (x, y), radius=sim.particles[i].r/self.roche, color="grey", alpha=0.5)
                        ax.add_artist(circ)

                    ax.set_xlim(-4, 4)
                    ax.set_ylim(-4, 4)
                    ax.set_aspect("equal")

                    writer.grab_frame()
                    clear_output(wait=True)
                    fig.clf()
                    print("Time: ", t, "/", self.sim.tmax, end="\r")
        else: # View of the disk of R-Z plane (side view)
            FFMpegWriter = manimation.writers["ffmpeg"]
            metadata = dict(title=title, comment=comment)
            writer = FFMpegWriter(fps=15, metadata=metadata)
            fig = plt.figure(figsize=(10, 5))

            with writer.saving(fig, vid_title+".mp4", dpi=dpi):
                for t in range(0, int(self.sim.tmax), int(0.05*self.t_kepler)):
                    sim = self.sim.getSimulation(t)
                    z_coords = []
                    r_coords = []
                    radii = []

                    # Calculate z and r coordinates, also radius for each particle
                    z_temp = []
                    r_temp = []
                    radii_temp = []
                    for j in range(1, sim.N):
                        x = (sim.particles[j].x - sim.particles[0].x)/self.roche
                        y = (sim.particles[j].y - sim.particles[0].y)/self.roche
                        z = (sim.particles[j].z - sim.particles[0].z)/self.roche
                        z_temp.append(np.abs(z))
                        r_temp.append(self.mag(x, y, z))
                        radii_temp.append(sim.particles[j].r/self.roche)
                    z_coords.append(z_temp)
                    r_coords.append(r_temp)
                    radii.append(radii_temp)

                    fig, axs = plt.subplots()# 10, 12

                    # Top figure
                    axs.set_title("Test")
                    axs.axvline(0.35, ls="--", color="gray")
                    axs.axvline(0.35, ls="--", color="gray")
                    axs.tick_params(direction="in", labelsize=textsize)

                    x0, y0 = axs.transAxes.transform((0, 0))
                    x1, y1 = axs.transAxes.transform((1.25, 3.0))
                    dx = x1 - x0
                    dy = y1 - y0
                    maxd = max(dx, dy)
                    for r, z, rad in zip(r_coords[0], z_coords[0], radii[0]):
                        width = rad * maxd / dx
                        height = rad * maxd / dy
                        axs.add_artist(Ellipse((r, z), width, height, alpha=0.7))
                    axs.set_xlabel("r $(a_R)$")
                    axs.set_ylabel("z $(a_R)$")

                    writer.grab_frame()
                    clear_output(wait=True)
                    fig.clf()
                    print("Time: ", t, "/", self.sim.tmax, end="\r")


        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Options for rebound simulation archive analysis')
    parser.add_argument('file_name', nargs='?', default='simulationarchive_testparticles.bin', type=str,
                             help='Name of rebound simulation archive')
    parser.add_argument('time', nargs='?', default=0.0, type=float,
                             help='Time for plotting and analysis')
    args = parser.parse_args()

    t = np.array([0, 25, 40]) 
    sim = ReboundAnalysis(args.file_name, args.time)
    #sim.plot_orbs(t, "Ida9_test_newradius_363_Massive", 150, False)
    sim.make_video("test", "test", "Tidal merging 1500 particles")
