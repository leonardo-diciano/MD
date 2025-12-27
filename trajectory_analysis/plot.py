import numpy as np
import matplotlib.pyplot as plt
import sys

track_e_tot, track_e_kin, track_e_pot, track_f_norm,track_temp,track_pressure, savefigs, helpmode = False,False,False,False,False,False,False,False

for flag in sys.argv:
    if flag == "-h":
        helpmode = True
        print("This script requires numpy and matplotlib")
        print("use it as:")
        print("python3 plot.py properties.txt")
        print("you may use different flags to track different properties:")
        print("-e_tot")
        print("-e_pot")
        print("-e_kin")
        print("-f_norm")
        print("-temp")
        print("-pressure")
        print("")
        print("To track all components of the energy at once, use:")
        print("-energies")
        print("")
        print("The most important flag, decide if you want the plots pop up as windows (default) or saved as figures (-savefig)")
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
    if flag == "-energies":
        track_e_tot = True
        track_e_kin = True
        track_e_pot = True
    if flag == "-savefig":
        savefigs = True

if not helpmode:
    properties_outfile = sys.argv[1]
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
        if len(properties) != 7:
            print("WARNING: this analysis is for txt files containing the properties ['istep', 'E_tot', 'E_kin', 'E_pot', 'F_norm', 'Temp', 'Pressure'] only")

        istep,E_tot,E_kin,E_pot,F_norm,Temp,Pressure = [],[],[],[],[],[],[]
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

    if (track_e_tot and track_e_kin and track_e_pot):
        plt.figure()
        plt.plot(istep, E_tot, label = "E_tot")
        plt.plot(istep, E_pot, label = "E_pot")
        plt.plot(istep, E_kin, label = "E_kin")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("Energy throughout the simulation")
        plt.legend()
        if (savefigs):
            plt.savefig(properties_outfile+".energy_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_e_tot and not(track_e_kin and track_e_pot)):
        plt.figure()
        plt.plot(istep, E_tot, label = "E_tot")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("Total energy throughout the simulation")
        if (savefigs):
            plt.savefig(properties_outfile+".etot_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_e_pot and not(track_e_kin and track_e_tot)):
        plt.figure()
        plt.plot(istep, E_pot, label = "E_pot")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("potential energy throughout the simulation")
        if (savefigs):
            plt.savefig(properties_outfile+".epot_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_e_kin and not(track_e_tot and track_e_pot)):
        plt.figure()
        plt.plot(istep, E_kin, label = "E_kin")
        plt.ylabel("E [kJ/mol]")
        plt.xlabel("istep")
        plt.title("kinetic energy throughout the simulation")
        if (savefigs):
            plt.savefig(properties_outfile+".ekin_plot.png", dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_temp):
        plt.figure()
        plt.plot(istep, Temp, label = "Temp [K]")
        plt.ylabel("T [K]")
        plt.xlabel("istep")
        plt.title("Temperature throughout the simulation")
        if (savefigs):
            plt.savefig(properties_outfile+".temperature_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_pressure):
        plt.figure()
        plt.plot(istep, Pressure, label = "Pressure [Pa]")
        plt.ylabel("P [Pa]")
        plt.xlabel("istep")
        plt.title("Pressure throughout the simulation")
        if (savefigs):
            plt.savefig(properties_outfile+".pressure_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()

    if (track_f_norm):
        plt.figure()
        plt.plot(istep, F_norm, label = "F_norm [kJ/(Å mol)]")
        plt.ylabel("F [kJ/(Å mol)]")
        plt.xlabel("istep")
        plt.title("F_norm throughout the simulation")
        if (savefigs):
            plt.savefig(properties_outfile+".forces_plot.png",dpi=300)
            plt.close()
        else:
            plt.show()
