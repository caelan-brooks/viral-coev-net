import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from coevolution_network_base import Population

# Close all open figures
plt.close('all')

# Initialize parameters
L = 40.0
dx = 0.3
x = np.arange(-L/2, L/2, dx)
r = 3
M = 15
beta = 2
alpha = 1
gamma = 0
D = 0.01
Nh = 2 * 10**6
viral_density = np.where(np.abs(x) <= 0.5 , 100.0, 0)
immune_density = np.zeros(len(x))  # Initializing immune_density as zeros

# Create a population object
population = Population( 
                        L, 
                        dx, 
                        r, 
                        M, 
                        beta, 
                        alpha, 
                        gamma, 
                        D, 
                        Nh, 
                        viral_density, 
                        immune_density, 
                        stochastic=True
                        )

#############################
## ANIMATION CODE
#############################

# Setting up the plotting
fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.set_xlabel('x')
ax1.set_ylabel('Viral Density', color=color)
ax1.set_xlim(min(population.xs), max(population.xs))
ax1.set_ylim(1, Nh)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_yscale('log')

ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylabel('Immune Density', color=color)
ax2.set_ylim(0, 1)
ax2.tick_params(axis='y', labelcolor=color)

line1, = ax1.plot([], [], color='tab:blue')
line2, = ax2.plot([], [], color='tab:red')

# Initializing lines
def init():
    line1.set_data(population.xs, [])
    line2.set_data(population.xs, [])
    return line1, line2

# Updating lines
def update(frame):
    population.single_step_evolve(dt=0.05)
    line1.set_ydata(population.viral_density)
    line2.set_ydata(population.immune_density)

    # Adding text to show the timestep
    ax1.text(0.1, 0.9, f'Timestep: {frame}', transform=ax1.transAxes, color='green')
   
    return line1, line2

ani = FuncAnimation(fig, update, frames=range(1600), init_func=init, blit=True, interval=5)

plt.title('Viral and Immune Density Evolution')
plt.show()
