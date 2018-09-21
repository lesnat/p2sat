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
        if verbose: print("Exporting %s phase space in %s ..."%(self._ps.specie["name"],file_name))

        d=self._ps.data

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

    def TrILEns_input(self,path,with_time=True,verbose=True):
        """
        Export particle phase space in a TrILEns input file.

        Parameters
        ----------
        path : str
            path to the output folder
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is exported
        """
        if verbose: print("Exporting %s phase space in %s ..."%(self._ps.specie["name"],path+"prop_th.t"))

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
