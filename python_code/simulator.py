class Simulation():
    def __init__(self, T, dt, L, dx, r, M, beta, alpha, gamma, D, migration_rate, Nh, initial_infected, initial_width):
        """Initialize a simulation"""
        self.T = T 
        self.dt = dt
        self.L = L
        self.dx = dx 
        self.r = r
        self.M = M
        self.beta = beta
        self.alpha = alpha
        self.gamma = gamma
        self.D = D
        self.migration_rate = migration_rate
        self.Nh = Nh
        self.initial_infected = initial_infected
        self.initial_width = initial_width

    
