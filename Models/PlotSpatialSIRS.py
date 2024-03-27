import matplotlib.pyplot as plt
import numpy as np

# Plots the Contour plot of the infection data for different combinations of p1 & p3
def ContourPlot():
    filename = 'Models/Data_Files/Infected_Data.txt'
    f = open(filename, 'r')
    lines = f.readlines() 
    
    p1s = np.zeros(len(lines))
    p3s = np.zeros(len(lines))
    Is = np.zeros(len(lines))

    i = 0
    for line in lines:
        line = line.strip('\n').split(',')
        p1s[i] = float(line[0])
        p3s[i] = float(line[1])
        Is[i] = float(line[2])
        i += 1

    f.close()

    fig,ax = plt.subplots(figsize=(8,8))    

    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlabel('P(S->I)')
    ax.set_ylabel('P(R->S)')
    ax.set_title("Fractional Average Infected Sites, <I>/N")
    ax.set_aspect('equal')
    tcf = ax.tricontour(p1s,p3s,Is)
    ax.tricontour(p1s,p3s,Is)
    fig.colorbar(tcf, ax=ax)

    plt.savefig('Documents/Images/ContourPlot.png')
    plt.show()

# Plots the Variance in the Average Infected fraction against p1
def VariancePlot():
    filename = 'Models/Data_Files/Variance_Data.txt'
    f = open(filename, 'r')
    lines = f.readlines() 
    
    p1s = np.zeros(len(lines))
    Vs = np.zeros(len(lines))
    Errs = np.zeros(len(lines))

    i = 0
    for line in lines:
        line = line.strip('\n').split(',') 
        p1s[i] = float(line[0])
        Vs[i] = float(line[1])
        Errs[i] = float(line[2])
        i += 1

    f.close()

    fig,ax = plt.subplots(figsize=(8,8))    

    ax.set_xlabel('P(S->I)')
    ax.set_ylabel('Variance')
    ax.set_title("Scaled Variance vs P1 (P2 = P3 = 0.5)")
    ax.errorbar(p1s,Vs, Errs, fmt='kx')
    
    plt.plot(p1s,Vs, linestyle='dashed')
    plt.savefig('Documents/Images/VariancePlot.png')
    plt.show()

# Plots the fraction of permanent immunity vs average infected fraction
def HerdPlot():
    filename = 'Models/Data_Files/Herd_Immunity_Data.txt'
    f = open(filename, 'r')
    lines = f.readlines() 
    
    hs = np.zeros(len(lines))
    Is = np.zeros(len(lines))

    i = 0
    for line in lines:
        line = line.strip('\n').split(',') 
        hs[i] = float(line[0])
        Is[i] = float(line[1])
        i += 1

    f.close()

    fig,ax = plt.subplots(figsize=(8,8))    

    ax.set_xlim(0,1)
    ax.set_ylim(0,0.5)
    ax.set_xlabel('Immunity Fraction (f_Im)')
    ax.set_ylabel('Average Infected Fraction (<I>/N)')
    ax.set_title("Average Infected Fraction vs Immunity Fraction (P1 = P2 = P3 = 0.5)")

    plt.plot(hs,Is)
    plt.savefig('Documents/Images/HerdPlot.png')
    plt.show()

def DataPlot():
    filename = 'Models/Data_Files/Wave_Data.txt'
    f = open(filename, 'r')
    lines = f.readlines() 
    
    Ss = np.zeros(len(lines))
    Is = np.zeros(len(lines))
    Rs = np.zeros(len(lines))
    ts = np.arange(0,len(lines))

    i = 0
    for line in lines:
        line = line.strip('\n').split(',') 
        Ss[i] = float(line[0]) * 100
        Is[i] = float(line[1]) * 100
        Rs[i] = float(line[2]) * 100
        i += 1

    f.close()

    fig,ax = plt.subplots(figsize=(18,8))  
        
    ax.set_xlim(100,900)
    ax.set_ylim(0,20)
    ax.set_xlabel('Time (Days)')
    ax.set_ylabel('Infected Percentage of Population')
    ax.set_title("Spatial SIRS (Beta = 1, Gamma = 0.1, Alpha = 0.01)", fontweight='bold')


    plt.rcParams.update({'font.size': 24})
    plt.rc('axes', labelsize=20) 
    for label in (ax.get_xticklabels() + ax.get_yticklabels()): 
        label.set_fontsize(16)

    #plt.plot(ts,Ss,label="Susceptible")
    plt.plot(ts,Is,label="Infected",color='#9fda9a')
    #plt.plot(ts,Rs,label="Resistant")

    #plt.legend()
    plt.tight_layout()
    plt.savefig('Documents/Images/InfectionPlot.png')
    plt.show()

# Allows command line inputs to determine choice of plot
if __name__ == "__main__":
    opt = input("Plot Type - Contour (C), Variance (V), Herd Immunity (H), Full SIR Data (S): ").upper()
    if opt == 'C':
        ContourPlot()
    elif opt == 'V':
        VariancePlot()
    elif opt == 'H':
        HerdPlot()
    elif opt == 'S':
        DataPlot()