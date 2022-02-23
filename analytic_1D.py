"""

From A novel lattice boltzmann model for the poisson equation

"""
import numpy as np 
from scipy.special import erfc
import matplotlib.pyplot as plt

def read_xml(f_name : str, fields : list) -> dict:
    retDict = {}
    with open(f_name,'r') as f:
        text = f.read().split("\n")
    for line in text:
        for field in fields:
            if "<{}>".format(field[0]) in line:
                temp = list(filter(None,line.split(" ")))
                retDict[field[0]] = field[1](temp[1])
    return retDict

def read_lbm_file(f_name : str, nx : int) -> np.array:
    """
    This function reads the *.dat file that Palabos writes on the
    """
    with open(f_name,'r') as f:
        text = f.read().split(" ")
    temp = np.array([float(j) for j in text[0:-1]]).reshape(1,len(text)-1)
    return temp

def u(x,k):
    return (np.exp(k)-1)/(np.exp(k)-np.exp(-k))*np.exp(-k*x)+(1-np.exp(-k))/(np.exp(k)-np.exp(-k))*np.exp(k*x)

if __name__ == "__main__":
    fields = [
                ('lx',int),
                ('ly',int),
                ('resolution',int),
                ('tau_phi',float),
                ('K_0',float),
            ]
    params = read_xml("params.xml",fields)

    resolution = params['resolution']
    nx = int(params['lx']*resolution)
    ny = int(params['ly']*resolution)
    K_0 = params['K_0']
    tau = params['tau_phi']

    k = 27.79
    x = np.linspace(0,1,1000)

    lbm_data=read_lbm_file("concentration_final.dat",ny).flatten()
    x2 = np.linspace(0,1,len(lbm_data))

    plt.figure()
    ax = plt.gca()
    ax.plot(x,u(x,k),linewidth = .5, label='analytic')
    ax.plot(x2,lbm_data,'kx',linewidth = .5, label ='lbm', markersize = 3)
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('u(x)')
    plt.savefig("figure1.png",dpi = 300)

    plt.show()
