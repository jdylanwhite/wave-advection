# Import modules
import numpy as np
import matplotlib.pyplot as plt

# Load data
x = np.loadtxt("./data/x.dat")
y = np.loadtxt("./data/y.dat")
t = np.loadtxt("./data/t.dat")
psi = np.loadtxt("./data/psi.dat")

# Reshape psi
psi = psi.reshape(len(t),len(x),len(y))

# Create meshgrid
X, Y = np.meshgrid(x,y)

# Set up the figure
fig, ax = plt.subplots(figsize=(8,8))

# Loop through data and plot data
for n,time in enumerate(t):

    # Plot the data
    pc = ax.contourf(X,Y,psi[n,:,:],cmap="Blues")
    
    # Set titles
    ax.set_title(f"Time: {time}")
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    # Set axis limits
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])

    # Set the color bar
    if n==0:
        cbar = fig.colorbar(pc)
    else:
        cbar.update_normal(pc)
    
    # Pause the plot for a moment
    if n == 0:
        plt.pause(2)
    else:
        plt.pause(0.01)

    # Clear the plot to update the results
    if n != len(t)-1:
        ax.clear()
        
plt.show()

