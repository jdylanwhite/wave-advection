# Import modules
import numpy as np
import matplotlib.pyplot as plt

# Load data from output file
xData = np.loadtxt('./x.dat')
tData = np.loadtxt('./t.dat')
uData = np.loadtxt('./u.dat')

# Set up the figure
fig, ax = plt.subplots(figsize=(8,6))

# Loop through data and plot data
for n,t in enumerate(tData):

    # Plot the data
    ax.plot(xData,uData[n,:])

    # Set titles
    ax.set_title(f"Time: {t}")
    ax.set_xlabel('x')
    ax.set_ylabel('Wave amplitude')

    # Set axis limits
    ax.set_xlim(xData[0], xData[-1])
    ax.set_ylim(-1.5, 1.5)
    
    # Pause the plot for a moment
    if n == 0:
        plt.pause(2)
    else:
        plt.pause(0.01)

    # Clear the plot to update the results
    if n != len(tData)-1:
        ax.clear()
        
plt.show()
