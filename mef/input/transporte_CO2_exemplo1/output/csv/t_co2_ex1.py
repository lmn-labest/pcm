import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def peclet_mesh(D, vel, h):

    return 0.5e0*vel*h/D


def func(x, D, vel, L):

    tmp1 = np.exp(vel*L/D) - 1
    tmp2 = np.exp(vel*x/D) - 1

    return 1.e0 - tmp2/tmp1


def main():


    for letra, vel in zip(('a', 'b'), (1.0,0.5)):

        data = pd.read_csv("t_co2_ex1_supg_" + letra + ".csv")
        x_s, y_s = data["Points:0"].values, data["cCO2"].values

        data = pd.read_csv("t_co2_ex1_garlekin_" + letra + ".csv")
        x_g, y_g = data["Points:0"].values, data["cCO2"].values

        x_a = np.linspace(0.0,1.0,200)
        y_a = func(x_a,0.04, vel, 1.0)

        Pe = peclet_mesh(0.04, vel, 0.1)

        fig, ax = plt.subplots()

        ax.plot(x_s, y_s, color='red', ls='--', label='SUPG')
        ax.plot(x_g, y_g, color='blue', ls='--', label='Garlekin')
        ax.plot(x_a, y_a, color='black', ls='-', label='Analítico')
        ax.set_title(r'Resultado para ${Pe_m} = $' + '{0}'.format(Pe))

        ax.set_xlim(0,1)
        ax.set_ylim(0,1.2)
        ax.set_xlabel('x')
        ax.set_ylabel(r'$\overline{C}$')
        ax.legend(fontsize=12)

        ax.grid(ls="--",color="gray")

        fig.savefig('t_co2_x1_'+letra+'.pdf', dpi=300)
        fig.savefig('t_co2_x1_'+letra+'.png', dpi=300)


if __name__ == "__main__":
    main()