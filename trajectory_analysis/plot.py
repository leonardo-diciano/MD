import sys

track_e_tot, track_e_kin, track_e_pot, track_f_norm,track_temp,track_pressure, savefigs, helpmode, track_com_mom, track_bias, track_metaG = False,False,False,False,False,False,False,False,False,False, False

for flag in sys.argv:
    if flag == "-h":
        helpmode = True
        print("usage: python3 plot.py [-h] 'properties_file' [-energies] [-e_temp] [...] [-savefig]")
        print("\nyou may use different flags to track different properties:")
        print("     -e_tot          plot the total energy")
        print("     -e_pot          plot the potential energy")
        print("     -e_kin          plot the kinetic energy")
        print("     -energies       plot all three energy components in one window")
        print("     -f_norm         plot the norm of the force")
        print("     -temp           plot the temperature")
        print("     -pressure       plot the pressure")
        print("     -com_mom        plot the norm of the momentum of the center of mass (to check for translational motion)")
        print("     -bias           plot the total bias potential evolution")
        print("     -metaG          plot the reconstructed free energy surface from metadynamics")
        print("")
        print("     -savefig        instead of pop up windows, save a png of the plot")
    if flag == "-e_tot":
        track_e_tot = True
    if flag == "-e_kin":
        track_e_kin = True
    if flag == "-e_pot":
        track_e_pot = True
    if flag == "-f_norm":
        track_f_norm = True
    if flag == "-temp":
        track_temp = True
    if flag == "-pressure":
        track_pressure = True
    if flag == "-com_mom":
        track_com_mom = True
    if flag == "-bias":
        track_bias = True
    if flag == "-metaG":
        track_metaG = True
    if flag == "-energies":
        track_e_tot = True
        track_e_kin = True
        track_e_pot = True
    if flag == "-savefig":
        savefigs = True



import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def fes_reweight(cv, bias, temperature, gamma, nbins=500, discard=0.1):
    kB = 8.314462618e-3
    beta = 1.0 / (kB * temperature)

    # Discard initial part
    a=int(discard*len(cv))
    print(np.min(cv), np.argmin(cv),np.max(cv),np.argmax(cv))
    cv = np.degrees(np.asarray(cv[a:], dtype=float))
    bias = np.asarray(bias[a:], dtype=float)

    # Reweight the probability distribution and evaluate the histogram
    # with a Gaussian KDE for smoother FES
    weights = np.exp(beta * bias / gamma)
    kde = gaussian_kde(cv, weights=weights)
    s_grid = np.linspace(np.min(cv), np.max(cv), nbins)
    p_s = kde(s_grid)

    # Get the free energy and shift it to zero
    fes = -kB * temperature * np.log(p_s)
    fes -= np.nanmin(fes)

    return s_grid, fes



if not helpmode:
    properties_outfile = sys.argv[1]
    filename = ""
    print("reading from: ",properties_outfile)

    with open(properties_outfile, "r") as file:
        rows = []

        for row in file:
            rows.append(row.split())
        properties = rows[0]
        units = rows[1]

        if len(properties) != len(units):
            nprop = len(properties)
            nunits = len(units)
            print("WARNING: %i properties were read in, while %i units were given."%(nprop,nunits))
            print("properties: ",properties)
            print("units: ",units)
            print("check the input file if every property has a unit assigned to (or if it is unitless, it has a 'none')")
        if len(properties) != 11:
            print("WARNING: this analysis is for txt files containing the properties ['istep', 'E_tot', 'E_kin', 'E_pot', 'F_norm', 'Temp', 'Pressure', 'COM_momentum', 'CV_value', 'Inst_Bias', 'Tot_Bias'] only")

        istep,E_tot,E_kin,E_pot,F_norm,Temp,Pressure,com_mom,tot_bias,cv_value, inst_bias = [],[],[],[],[],[],[],[],[],[],[]
        for row in rows[2:]:
            istep.append(int(row[0]))
            if (track_e_tot):
                E_tot.append(float(row[1]))
            if (track_e_kin):
                E_kin.append(float(row[2]))
            if (track_e_pot):
                E_pot.append(float(row[3]))
            if (track_f_norm):
                F_norm.append(float(row[4]))
            if (track_temp):
                Temp.append(float(row[5]))
            if (track_pressure):
                Pressure.append(float(row[6]))
            if (track_com_mom):
                com_mom.append(float(row[7]))
            if (track_bias):
                tot_bias.append(float(row[10]))
            if (track_metaG):
                cv_value.append(float(row[8]))
                inst_bias.append(float(row[9]))

    if (track_e_tot and track_e_kin and track_e_pot):
        plt.figure()
        plt.plot(istep, E_tot, label = "E_tot")
        plt.plot(istep, E_pot, label = "E_pot")
        plt.plot(istep, E_kin, label = "E_kin")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("Energy throughout the simulation")
        plt.legend()
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"energy_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_e_tot and not(track_e_kin and track_e_pot)):
        plt.figure()
        plt.plot(istep, E_tot, label = "E_tot")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("Total energy throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"etot_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_e_pot and not(track_e_kin and track_e_tot)):
        plt.figure()
        plt.plot(istep, E_pot, label = "E_pot")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("potential energy throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"epot_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_e_kin and not(track_e_tot and track_e_pot)):
        plt.figure()
        plt.plot(istep, E_kin, label = "E_kin")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("kinetic energy throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"ekin_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_temp):
        plt.figure()
        plt.plot(istep, Temp, label = "Temp [K]")
        plt.ylabel("T [K]")
        plt.xlabel("istep")
        plt.title("Temperature throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"temperature_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_pressure):
        plt.figure()
        plt.plot(istep, Pressure, label = "Pressure [Pa]")
        plt.ylabel("P [Pa]")
        plt.xlabel("istep")
        plt.title("Pressure throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"pressure_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_f_norm):
        plt.figure()
        plt.plot(istep, F_norm, label = "F_norm [kJ/(Å mol)]")
        plt.ylabel("F [kJ/(Å mol)]")
        plt.xlabel("istep")
        plt.title("F_norm throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"forces_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()


    if (track_com_mom):
        plt.figure()
        plt.plot(istep, com_mom, label = "com_mom [gÅ/fs]")
        plt.ylabel("Momentum [gÅ/fs]")
        plt.xlabel("istep")
        plt.title("Total momentum (of the center of mass) throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"com_mom_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_bias):
        plt.figure()
        plt.plot(istep, tot_bias, label = "total bias potential [kJ/mol]")
        plt.ylabel("Total bias potential [kJ/mol]")
        plt.xlabel("istep")
        plt.title("Total bias potential applied throughout the simulation")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"tot_bias_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_metaG):
        bin,fes = fes_reweight(cv_value,inst_bias,298,10)
        plt.figure()
        plt.plot(bin, fes, label = "Free Energy[kJ/mol]")
        plt.ylabel("Free Energy [kJ/mol]")
        plt.xlabel("CV value")
        plt.title("Reconstructed Free Energy surface from Well-Tempered Metadynamics")
        plt.tight_layout()
        if (savefigs):
            plt.savefig(filename+"fes_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()



