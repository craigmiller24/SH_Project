import matplotlib.pyplot as plt
import numpy as np

# Plots the Contour plot of the infection data for different combinations of p1 & p3
def ContourPlot():
    filename = 'Models/Data_Files/Infected_Data.txt'

    # Read the data from the .txt file
    data = np.genfromtxt(filename, delimiter=',')

    # Separate the columns
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # Create an array of indices to sort the data based on x and y values
    sorted_indices = np.lexsort((x, y))

    # Sort the z values based on the sorted indices
    sorted_z = z[sorted_indices]


    # Reshape the sorted z values into a grid
    z_grid = sorted_z.reshape(36,36)

    fig,ax = plt.subplots(figsize=(8,8))    

    #ax.set_xlim(0,0.35)
    #ax.set_ylim(0,0.35)
    ax.set_xlabel('P(S->I)')
    ax.set_ylabel('P(R->S)')
    ax.set_title("Fractional Average Infected Sites, <I>/N")
    ax.set_aspect('equal')
    #ax.tricontour(x,y,z)

    # Plot the heatmap
    ht = ax.imshow(z_grid, cmap='gnuplot', origin='lower', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    fig.colorbar(ht,ax=ax,label="Value")  # Add a colorbar
    plt.savefig('Documents/Images/ContourPlot2.png')
    plt.show()

    #plt.show()

# Plots the Variance in the Average Infected fraction against p1
def VariancePlot():
    filename = 'Models/Data_Files/Variance_Data_updated.txt'
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
    plt.savefig('Documents/Images/VariancePlotNew.png')
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
    filename = 'Models/Data_Files/Wave_Data_2.txt'
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
    
    plt.figure(figsize=(15, 6))

    plt.xlim(100,900)
    plt.ylim(0,25)
    plt.xlabel('Time (Days)')
    plt.ylabel('Infected Percentage of Population')
    plt.title("Spatial SIRS Model")


    #plt.plot(ts,Ss,label="Susceptible")
    plt.plot(ts,Is,label="Infected",color='#9fda9a')
    #plt.plot(ts,Rs,label="Resistant")

    #plt.legend()
    plt.tight_layout()
    plt.savefig('Documents/Images/InfectionPlot2.png')
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