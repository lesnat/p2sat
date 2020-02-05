#coding:utf8
import numpy as np

class _SavePhaseSpace(object):
    r"""
    Export phase space into several file format.
    """
    def __init__(self,PhaseSpace):
        self._ds=PhaseSpace

    def txt(self,file_name,header=True,title="",sep=",",verbose=True):
        r"""
        Export particle phase space in a text file.

        Parameters
        ----------
        file_name : str
            name of the output file
        header : bool, optional
            True to put informations at the beginning of the file. Default is True
        title : str, optional
            short description of the file content
        sep : str, optional
            character to use to separate values. Default is ',' (csv file)
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is exported

        Notes
        -----
        The format in the output file is
        ::
            # title
            # legend
            w x y z px py pz t
            w x y z px py pz t
            . . . . .  .  .  .
            . . . . .  .  .  .
            . . . . .  .  .  .

        with 7 digits precision in scientific notation

        Some text can be written if the first character of the line is a '#'.
        """
        if verbose: print("Exporting %s phase space in %s ..."%(self._ds.metadata.specie["name"],file_name))

        d=self._ds.read

        # Opening the output file
        with open(file_name,'w') as f:
            if header:
                # Write title
                f.write("# Title : %s\n"%title)
                # Write legend
                f.write("# ")
                for legend in ('weight','x (um)','y (um)','z (um)','px (MeV/c)','py (MeV/c)','pz (MeV/c)','t (fs)'):
                    f.write("%-16s"%legend) # the chain is placed under 16 characters
                f.write("\n")

            # Write data
            for i in range(len(d.w)):
                for e in [d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],d.t[i]]:
                    tmp = "% .7E"%e # 7 digits precision with E notation
                    tmp+= sep       # separator
                    f.write("%-16s"%tmp) # the chain is placed under 16 characters
                f.write("\n")

        # np.savetxt(file_name,list(zip(d.w,d.x,d.y,d.z,d.px,d.py,d.pz,d.t)))

        if verbose: print('Done !')

    def gp3m2_input(self,file_name,title="",verbose=True):
        r"""
        Export particle phase space in a gp3m2 input file

        Parameters
        ----------
        file_name : str
            name of the file
        title : str, optional
            short description of the file content
        verbose : bool, optional
            verbosity
        """
        self.txt(file_name,header=True,title=title,sep=" ",verbose=verbose)

    def TrILEns_prop_ph(self,path,with_time=True,verbose=True):
        r"""
        Export particle phase space in a TrILEns input file.

        Parameters
        ----------
        path : str
            path to the output folder
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is exported
        """
        if verbose: print("Exporting %s phase space in %s ..."%(self._ds.metadata.specie["name"],path+"prop_ph.t"))

        d=self._ds.read

        # Opening the output file
        with open(path+'prop_ph.t','w') as f:
            if with_time:
                f.write('9 1.\n')
                f.write(' poi  phx  phy  phz  pdx  pdy  pdz  gph  tim\n')
            else:
                f.write('8 1.\n')
                f.write(' poi  phx  phy  phz  pdx  pdy  pdz  gph\n')

            # Write data
            for i in range(len(d.w)):
                if with_time:
                    data=[d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],d.ekin[i]/0.511,d.t[i]*1e-3]
                else:
                    data=[d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],d.ekin[i]/0.511]

                for e in data:
                    tmp="% .7E"%e # 7 digits precision with E notation
                    f.write("%-16s"%tmp) # the chain is placed under 16 characters
                f.write("\n")

        if verbose: print('Done !')

    def TrILEns_input(self,path,source_1,source_2,
                    lvl_diag=7, random_distrib=False, stat_col=True,
                    chi_to_time=True, normalized_nb=False,
                    pasdt=1.0,
                    max_lvl_kd_tree=18, max_lvl_octree=5, max_number_swarms=20,
                    threshold_merging=10, lvl_spatial=5, lvl_momentum=5,
                    maillage_spatial=None, nb_max_par_espece=2000000,
                    verbose=True):
        r"""
        Write default input file and export particle phase space in a TrILEns input file.

        Parameters
        ----------
        path : str
            path to the output folder
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is exported
        """
        if verbose: print("Writing input.txt and data in path %s ..."%(path))

        s1 = source_1
        s2 = source_2
        # Creation of the regular input file
        f = open(path+"input.txt", 'w')
        # Level of diagnistics : which lvl of diagnostics? 0 in case of no diagnostic. 7 Provides minimum output and should be used by default.
        # Zombies : Saving "lost"/annihilated particles as Zombies. Out-of-date. Should always be 0
        # Random : Photon distribution inside a MP randomized? Randomized positioning is somewhat slower. Only relevant with StatisticCollisions = false. The default value is 1.
        # Statistic : Only statistical handling of collisions. No "real" photons.
        f.write("--------Level of diagnostic (0 means no diagnostic), Zombies: yes/no, random distribution of photons, statistic collisions?\n")
        f.write("{0} {1} {2} {3}\n".format(lvl_diag                 , False        , random_distrib               , stat_col))
        # Flush : Defines, whether particles, that left the volume or annihilated are deleted (or saved as zombies) or kept until the end. Out-of-date. Should always be 1
        # Chi to t : If true, then in a_chi of the partic value, the collision time will be stored. Their position will be the position of the collision! The default value is 1
        # Normalized : Set to 1, if normalised particle numbers should be used, just as in CALDER. The default value is 0.
        f.write("--------flush dead particles, Chi-To-CollisionTime, normalised particle numbers?\n")
        f.write("{0} {1} {2}\n".format(True, chi_to_time, normalised_nb))
        # Timestep :
        #The time step will be the time it takes for the fastest particle to cover dtFactor times the length of the shortest bounding volume in each direction
        #in case of PIC code, dtFactor will be dt. A reasonable default is 1. or 2.
        f.write("--------temps: pasdt\n")
        f.write("{0}\n".format(pasdt))
        # Number of photon species
        f.write("--------radiations: nesp_ph\n")
        f.write("{0}\n".format(1))
        # Index of electron and positron
        f.write("--------paires: ipaire_e, ipaire_p\n")
        f.write("{0} {1}\n".format(1, 2))
        # Number of species ???
        f.write("--------particules: nesp\n")
        f.write("{0}\n".format(2))
        # ???
        f.write("--------particules: a_po, a_qx, a_qy, a_qz, a_px, a_py, a_pz, a_ga, a_ch, a_qxm1, a_qym1, a_qzm1\n")
        f.write("{0:d} {1:d} {2:d} {3:d} {4:d} {5:d} {6:d} {7:d} {8:d} {9:d} {10:d} {11:d}\n".format(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
        # Total number of input particles
        f.write("--------particules: nbpa\n")
        f.write("0 0 {0}\n".format(len(s1)+len(s2)))
        # ...
        # Species consist of density, charge, mass, number of cycles for all species. Number of cycles is not used at all (CALDER artifact)
        f.write("--------particules: density, charge, mass (one triplet per attribute)\n")
        f.write("1 1 1\n")
        f.write("1 -1 0\n")
        f.write("1 1 0\n")
        f.write("1 1 1\n") # ???
        # Simulation ...
        #max level of the BVtree. After destruction of nodes, the level can decrease. Lvl 18 corresponds to 131.072 leaf nodes per swarm.
        #Limit for the number of swarms before merging is done
        #max level of the octree. The octrees are always complete and full. Lvl 7 corresponds to 32.768 leaf nodes.
        f.write("--------BVtree(yes/no), Max Lvl of k-d tree, Max Lvl of octree, Max number of swarms\n")
        f.write("{0} {1} {2} {3}\n".format(True, max_lvl_kd_tree, max_lvl_octree, max_number_swarms))
        f.write("--------number of lvls for output, number of outputs, array with levels to be given out\n")
        f.write("{0} {1}\n\n".format(0, 0)) # ???
        # ...
        #Threshold for the number of particles per spatial cell in order to start merging. A number smaller 3 does not make any sense
        #The level of the octree, that is used in order to equidistantly divide momentum space into cells. Lvl 1 has 1 leaf, level 2 has 8.
        #Same, but for spatial cells.
        f.write("--------Threshold for particle merging, level of division of spatial- and momentum space\n")
        f.write("{0} {1} {2}\n".format(threshold_merging, lvl_spatial, lvl_momentum))
        # ...
        f.write("--------number of likely beams and their main momentum and gamma-factor:\n{0}\n".format(2))
        f.write("{0} {1} {2} {3}\n".format(s1.read.px.mean(), s1.read.py.mean(), s1.read.pz.mean(), 1. + s1.read.ekin.mean()/0.511))
        f.write("{0} {1} {2} {3}\n".format(s2.read.px.mean(), s2.read.py.mean(), s2.read.pz.mean(), 1. + s2.read.ekin.mean()/0.511))
        # Restart simulation
        #Important for Reruning simulations. MakeResumepoints makes the simulation create points from which a Simulation can be resumed. ResumepointPeriod defines the period of
        #those points in simulation steps. RestartRun makes a simulation use one of such points. Input.txt is still necessary!
        f.write("--------Re-start simulation, give-out resume points, formatted output, resume point period in time steps.\n")
        f.write("{0} {1} {2} {3}\n".format(False, False, False, 5))
        # ...
        f.write("--------Mesh parameters. Auto-mesh, lower corner, higher corner, cell size:\n")
        f.write("{0}\n{1} {2} {3}\n{4} {5} {6}\n{7} {8} {9}\n".format(False, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        f.close()

        # Writing data
        f = open(path+"data", 'w')
        # ...
        #Size of each MP (and therefor, at the moment, of each particle) all lengths are given in microns
        f.write("maillage_spatial\n******** taille des mailles en x, y, z (c/w0) -----------------------------\n")
        if maillage_spatial is None:
            f.write("{0} {1} {2}".format(1.0, 1.0, 1.0))
        else:
            f.write("{0} {1} {2}".format(maillage_spatial[0], maillage_spatial[1], maillage_spatial[2]))
        f.write("\n\n")
        # ...
        f.write("src_laser\n******** lambda, a_0, angle_y, angle_z, polar, chirp, c_e_p----------------\n")
        f.write("{0} 0 0 0 0 0 0 0".format(1))
        f.write("\n\n")
        # ...
        # Used as a threshold for particle merging, but also in order to conserve memory
        f.write("especes_plasma\n******** nombre d'especes et nb max par espece ----------------------------\n")
        f.write("0 0 {0} {1}".format(1, nb_max_par_espece))
        f.close()

        if verbose:print('Done !')
