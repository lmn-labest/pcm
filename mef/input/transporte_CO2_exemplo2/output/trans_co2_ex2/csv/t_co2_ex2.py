import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():

    data = pd.read_csv("t_co2_ex2_25_supg_alf_1.csv")
    x_mef_1, y_mef_1 = data["Points:0"].values, data["cCO2"].values

    data = pd.read_csv("t_co2_ex2_25_supg_alf_05.csv")
    x_mef_05, y_mef_05 = data["Points:0"].values, data["cCO2"].values

    data = pd.read_csv("t_co2_ex2_25_openfoam_25.csv")
    x_ofoam, y_ofoam = data["arc_length"].values, data["T"].values

    fig, ax = plt.subplots()

    ax.plot(x_mef_1 , y_mef_1 , color='red', ls='--', label=r'$\alpha=1.0$')
    ax.plot(x_mef_05, y_mef_05, color='blue', ls='--', label=r'$\alpha=0.5$')
    ax.plot(x_ofoam, y_ofoam, color='black', ls='-', label='OpenFOAM')


    ax.set_title(r'Resultado para $t= 0.25$')

    ax.set_xlim(0,1)
    ax.set_ylim(0,1.2)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\overline{C}$')
    ax.legend(fontsize=12)

    ax.grid(ls="--",color="gray")

    fig.savefig('t_co2_x2_25.pdf', dpi=300)
    fig.savefig('t_co2_x2_25.png', dpi=300)

    data = pd.read_csv("t_co2_ex2_50_supg_alf_1.csv")
    x_mef_1, y_mef_1 = data["Points:0"].values, data["cCO2"].values

    data = pd.read_csv("t_co2_ex2_50_supg_alf_05.csv")
    x_mef_05, y_mef_05 = data["Points:0"].values, data["cCO2"].values

    data = pd.read_csv("t_co2_ex2_25_openfoam_50.csv")
    x_ofoam, y_ofoam = data["arc_length"].values, data["T"].values

    fig, ax = plt.subplots()

    ax.plot(x_mef_1 , y_mef_1 , color='red', ls='--', label=r'$\alpha=1.0$')
    ax.plot(x_mef_05, y_mef_05, color='blue', ls='--', label=r'$\alpha=0.5$')
    ax.plot(x_ofoam, y_ofoam, color='black', ls='-', label='OpenFOAM')


    ax.set_title(r'Resultado para $t= 0.5$')

    ax.set_xlim(0,1)
    ax.set_ylim(0,1.2)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\overline{C}$')
    ax.legend(fontsize=12)

    ax.grid(ls="--",color="gray")

    fig.savefig('t_co2_x2_55.pdf', dpi=300)
    fig.savefig('t_co2_x2_55.png', dpi=300)


if __name__ == "__main__":
    main()