import numpy as np


class Simulation:
    def __init__(self, T, dt, L, dx, r, M, beta, alpha, gamma, D, migration_rates, Nh, initial_infected, initial_width):
        """Initialize a simulation."""

        # Basic parameters
        self.T = T  # Total simulation time
        self.dt = dt  # Time step
        self.L = L  # Length of antigenic space
        self.dx = dx  # Discretization of antigenic space
        self.r = r  # Cross-reactivity
        self.M = M  # Number of immune pressures per host
        self.beta = beta  # Transmission rate
        self.alpha = alpha  # Death rate
        self.gamma = gamma  # Recover rate
        self.D = D  # Mutation rate (diffusion rate)
        self.migration_rates = migration_rates  # Migration rates
        self.Nh = Nh  # Number of hosts per deme
        self.initial_infected = initial_infected  # Initial number of infected individuals
        self.initial_width = initial_width  # Initial width of outbreak

        # Derived parameters
        self.ts = np.arange(0, T, dt)  # Vector of time points
        self.num_time_points = np.size(self.ts)  # Number of time points
        self.xs = np.arange(-L/2, L/2, dx)  # Vector of antigenic points
        self.num_antigen_points = np.size(self.xs)  # Number of antigenic points
        self.num_demes = np.shape(migration_rates)[0]  # Number of demes

        # Viral density "n_i(x,t)"
        self.viral_density = np.zeros([
            self.num_demes,
            self.num_antigen_points,
            self.num_time_points
        ])

        # Immune density "h_i(x,t)"
        self.immune_density = np.zeros([
            self.num_demes,
            self.num_antigen_points,
            self.num_time_points
        ])
        
