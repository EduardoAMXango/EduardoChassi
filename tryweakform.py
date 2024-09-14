import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class BodyAndFrameFEM:
    def __init__(self, A, E, L, I, G, J, num_elements):
        self.A = A
        self.E = E
        self.L = L
        self.I = I
        self.G = G
        self.J = J
        self.num_elements = num_elements

    def shape_functions(self, xi):
        N_bar = np.array([0.5*(1-xi), 0.5*(1+xi)])
        N_beam = np.array([(1-xi)/2, (1+xi)/2])
        return N_bar, N_beam
    
    def weak_form_bar(self, u, v, f):
        du_dx = np.gradient(u, self.L/self.num_elements)
        dv_dx = np.gradient(v, self.L/self.num_elements)
        
        weak_form = self.E * self.A * np.dot(dv_dx, du_dx) * (self.L/self.num_elements) - np.dot(v, f)
        return weak_form

    def weak_form_beam(self, w, v, q):
        d2w_dx2 = np.gradient(np.gradient(w, self.L/self.num_elements), self.L/self.num_elements)
        d2v_dx2 = np.gradient(np.gradient(v, self.L/self.num_elements), self.L/self.num_elements)
        
        weak_form = self.E * self.I * np.dot(d2v_dx2, d2w_dx2) * (self.L/self.num_elements) - np.dot(v, q)
        return weak_form
    
    def analyze(self, F_axial, F_flexao, F_torsao):
        return (
            self.deformacao_axial(F_axial),
            self.flexao(F_flexao),
            self.torsao(F_torsao)
        )

    def deformacao_axial(self, F_axial):
        u = np.linspace(0, self.L, self.num_elements + 1)
        v = np.zeros_like(u)
        return self.weak_form_bar(u, v, F_axial)
    
    def flexao(self, F_flexao):
        w = np.linspace(0, self.L, self.num_elements + 1)
        v = np.zeros_like(w)
        return self.weak_form_beam(w, v, F_flexao)
    
    def torsao(self, F_torsao):
        return np.full(self.num_elements + 1, (F_torsao * self.L) / (self.G * self.J))
    
A = 0.01    # m^2
E = 210e9   # Pa
I = 1.6667e-5 # m^4
G = 81.2e9  # Pa
K = 1000    # N/m
rho = 7850  # kg/m^3
g = 9.81    # m/s^2
J = 1e-6    # m^4 (momento polar de inércia)
F1 = np.array([1000, 2000, 3000, 4000, 5000])
F2 = np.array([1000, 2000, 3000, 4000, 5000])
T = np.array([1000, 2000, 3000, 4000, 5000])
