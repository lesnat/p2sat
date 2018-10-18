#coding:utf8
import numpy as np
import matplotlib.pyplot as plt

class _Export(object):
    """
    Export phase space into several file format.
    """
    def __init__(self,PhaseSpace):
        self._ps=PhaseSpace

    def txt(self,file_name,header=True,title="",sep=",",verbose=True):
        """
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
        if verbose: print("Exporting %s phase space in %s ..."%(self._ps.particle["name"],file_name))

        d=self._ps.data.raw

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

        if verbose: print('Done !')

    def gp3m2_input(self,file_name,title="",verbose=True):
        """
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
        """
        Export particle phase space in a TrILEns input file.

        Parameters
        ----------
        path : str
            path to the output folder
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is exported
        """
        if verbose: print("Exporting %s phase space in %s ..."%(self._ps.particle["name"],path+"prop_ph.t"))

        d=self._ps.data

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
                    data=[d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],1.+d.ekin[i]/0.511,d.t[i]]
                else:
                    data=[d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],1.+d.ekin[i]/0.511]

                for e in data:
                    tmp="% .7E"%e # 7 digits precision with E notation
                    f.write("%-16s"%tmp) # the chain is placed under 16 characters
                f.write("\n")

        if verbose: print('Done !')

    # def TrILEns_prop_ph(self,path,with_time=True,other=None,R=None,T=None,htes=True,verbose=True):
    #     """
    #     Export particle phase space in a TrILEns input file.
    #
    #     Parameters
    #     ----------
    #     path : str
    #         path to the output folder
    #     with_time : bool, optional
    #         export phase-space with time coordinate. Default is True
    #     other : PhaseSpace, optional
    #         second gamma source. If None, only the phase space of current object is exported
    #     R : tuple of 3 float, optional
    #         rotation angles to pass to the `transformate` method. If None, no rotation
    #     T : tuple of 3 float, optional
    #         rotation angles to pass to the `transformate` method. If None, no translation
    #     htes : bool
    #         half tranform each source. Default is True
    #     verbose : bool, optional
    #         verbosity
    #
    #     Notes
    #     -----
    #
    #     """
    #     if other is None:
    #         sources = self._ps
    #     else:
    #         s1 = self._ps.copy()
    #         s2 = other.copy()
    #         if R is None: R=(0,0,0)
    #         if T is None: T=(0,0,0)
    #
    #         if htes:
    #             R1 = - np.array(R)/2.
    #             T1 = - np.array(T)/2.
    #             R2 = np.array(R)/2.
    #             T2 = np.array(T)/2.
    #         else:
    #             R1 = (0,0,0)
    #             T1 = (0,0,0)
    #             R2 = R
    #             T2 = T
    #         if verbose: print("Source 1 :")
    #         s1.data.transformate(R=R1,T=T1,verbose=verbose)
    #         if verbose: print("Source 2 :")
    #         s2.data.transformate(R=R2,T=T2,verbose=verbose)
    #
    #         sources = s1 + s2
    #
    #     d=sources.data
    #
    #     if verbose: print("Exporting %s phase space in %s ..."%(self._ps.particle["name"],path+"prop_th.t"))
    #
    #     # Opening the output file
    #     with open(path+'prop_ph.t','w') as f:
    #         if with_time:
    #             f.write('9 1.\n')
    #             f.write(' poi  phx  phy  phz  pdx  pdy  pdz  gph  tim\n')
    #         else:
    #             f.write('8 1.\n')
    #             f.write(' poi  phx  phy  phz  pdx  pdy  pdz  gph\n')
    #
    #         # Write data
    #         for i in range(len(d.w)):
    #             if with_time:
    #                 data=[d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],1.+d.ekin[i]/0.511,d.t[i]]
    #             else:
    #                 data=[d.w[i],d.x[i],d.y[i],d.z[i],d.px[i],d.py[i],d.pz[i],1.+d.ekin[i]/0.511]
    #
    #             for e in data:
    #                 tmp="% .7E"%e # 7 digits precision with E notation
    #                 f.write("%-16s"%tmp) # the chain is placed under 16 characters
    #             f.write("\n")
    #
    #     if verbose: print('Done !')

    def TrILEns_input(self,path,S1,S2,pasdt=1.0,maillage_spatial=None,with_time=True,verbose=True):
        """
        Write default input file and export particle phase space in a TrILEns input file.

        Parameters
        ----------
        path : str
            path to the output folder
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is exported
        """
        if verbose: print("Writing input.txt and data in path %s ..."%(path))

        d=self._ps.data

        # Creation of the regular input file
        f = open(path+"input.txt", 'w')
        # Level of diagnistics : which lvl of diagnostics? 0 in case of no diagnostic. 7 Provides minimum output and should be used by default.
        # Zombies : Saving "lost"/annihilated particles as Zombies. Out-of-date. Should always be 0
        # Random :
        # Statistic :
        f.write("--------Level of diagnostic (0 means no diagnostic), Zombies: yes/no, random distribution of photons, statistic collisions?\n")
        f.write("{0} {1} {2} {3}\n".format(7                        , False             , False                       , True))
        # Flush : Defines, whether particles, that left the volume or annihilated are deleted (or saved as zombies) or kept until the end. Out-of-date. Should always be 1
        # Chi to t :
        # Normalized :
        f.write("--------flush dead particles, Chi-To-CollisionTime, normalised particle numbers?\n")
        f.write("{0} {1} {2}\n".format(True, True, False))
        # Timestep :
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
        f.write("0 0 {0}\n".format(len(self._ps)))
        # ...
        f.write("--------particules: density, charge, mass (one triplet per attribute)\n")
        f.write("1 1 1\n")
        f.write("1 -1 0\n")
        f.write("1 1 0\n")
        f.write("1 1 1\n") # ???
        # Simulation ...
        f.write("--------BVtree(yes/no), Max Lvl of k-d tree, Max Lvl of octree, Max number of swarms\n")
        f.write("{0} {1} {2} {3}\n".format(True, 18, 5, 20))
        f.write("--------number of lvls for output, number of outputs, array with levels to be given out\n")
        f.write("{0} {1}\n\n".format(0, 0)) # ???
        # ...
        f.write("--------Threshold for particle merging, level of division of spatial- and momentum space\n")
        f.write("{0} {1} {2}\n".format(10, 5, 5))
        # ...
        f.write("--------number of likely beams and their main momentum and gamma-factor:\n{0}\n".format(2))
        f.write("{0} {1} {2} {3}\n".format(S1.data.px.mean(), S1.data.py.mean(), S1.data.pz.mean(), 1. + S1.data.ekin.mean()/0.511))
        f.write("{0} {1} {2} {3}\n".format(S2.data.px.mean(), S2.data.py.mean(), S2.data.pz.mean(), 1. + S2.data.ekin.mean()/0.511))
        # Restart simulation
        f.write("--------Re-start simulation, give-out resume points, formatted output, resume point period in time steps.\n")
        f.write("{0} {1} {2} {3}\n".format(False, False, False, 5))
        # ...
        f.write("--------Mesh parameters. Auto-mesh, lower corner, higher corner, cell size:\n")
        f.write("{0}\n{1} {2} {3}\n{4} {5} {6}\n{7} {8} {9}\n".format(False, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        f.close()

        # Writing data
        f = open(path+"data", 'w')
        # ...
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
        f.write("especes_plasma\n******** nombre d'especes et nb max par espece ----------------------------\n")
        f.write("0 0 {0} {1}".format(1, 2000000))
        f.close

        if verbose:print('Done !')

        self.TrILEns_prop_ph(path,with_time=with_time,verbose=verbose)
