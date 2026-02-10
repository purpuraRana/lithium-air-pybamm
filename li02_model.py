import pybamm
import numpy as np

class LiO2_1D(pybamm.BaseModel):
    """
    1D Li-O2 model
    """
    def __init__(self, name="Li-O2 1D"):
        super().__init__(name)

        # --------------------
        # Parameters
        # --------------------
        # Constants
        F = pybamm.Parameter("F")
        R = pybamm.Parameter("R")
        T = pybamm.Parameter("T")
        Eeq = pybamm.Parameter("Eeq")   # ORR equilibrium potential

        # Separator
        Lsep = pybamm.Parameter("Lsep")          # 2.5e-5
        eps_sep = pybamm.Parameter("eps_sep")    # 0.5

        # Cathode
        a0 = pybamm.Parameter("a0")     # Specific interfacial area of cathode
        D_O2 = pybamm.Parameter("D_O2") # O2 diffusion coefficient
        b = pybamm.Parameter("Bruggeman exponent")
        eps_s0 = pybamm.Parameter('eps_s0') # Initial volume fraction of solid phase
        r0 = pybamm.Parameter("r0")     # Cathode particle radius
        lm = pybamm.Parameter("lm")     # Mean product film thickness
        cO2_ext = pybamm.Parameter("cO2_ext")
        cO2_0 = pybamm.Parameter("cO2_0")
        S_O2 = pybamm.Parameter("S_O2")
        M_Li2O2 = pybamm.Parameter("M_Li2O2") # Li2O2 molecular weight
        rho_Li2O2 = pybamm.Parameter("rho_Li2O2") # Li2O2 density
        p = pybamm.Parameter("p")       # Geometrical factor
        R_film = pybamm.Parameter("R_film")
        beta = pybamm.Parameter("beta")
        k_c = pybamm.Parameter('k_c')
        Lc = pybamm.Parameter("Lc")
        sigma = pybamm.Parameter("sigma")
        kappa = pybamm.Parameter("kappa")
        D_Li = pybamm.Parameter("D_Li")
        tP = pybamm.Parameter("tP")     # Li transference number
        cLi2O2_max = pybamm.Parameter("c_Li2O2_max")
        
        # Anode
        k_a = pybamm.Parameter('k_a')

        # Other
        l_sigma = pybamm.Parameter("l_sigma")

        # Input
        Japp = pybamm.InputParameter("J") # A/m^2 (galvanostatic discharge current density)
        J_floor = pybamm.Parameter("J_floor")
        tau = pybamm.Parameter("tau_ramp")
        J = Japp * (J_floor + (1 - J_floor) * (1 - pybamm.exp(-pybamm.t / tau)))

        x = pybamm.SpatialVariable("x", domain="cathode")

        # State Variables
        cO2 = pybamm.Variable("O2 concentration", domain="cathode")
        eps = pybamm.Variable("Porosity", domain="cathode")
        cLi = pybamm.Variable("Lithium concentration", domain="cathode")
        cLi2O2 = pybamm.Variable("Li2O2 concentration", domain="cathode")

        eps_floor = pybamm.Scalar(1e-6)
        eps_safe = pybamm.maximum(eps, eps_floor)
        eps_cap = pybamm.Scalar(1e-6)
        #eps_clip = pybamm.minimum(eps_safe, (pybamm.Scalar(1) - eps_s0 - eps_cap))
        eps_clip = eps_safe

        # Algebraic Variables
        phi_s = pybamm.Variable("Solid potential", domain="cathode")
        phi_l = pybamm.Variable("Electrolyte potential", domain="cathode")

        # Derived Variables
        D_O2_eff = D_O2 * eps_clip**b
        D_Li_eff = D_Li * eps_clip**b
        sigma_eff = sigma * (1-eps_clip)**b
        kappa_eff = kappa * eps_clip**b
        eps_Li2O2 = 1 - eps - eps_s0  # Li2O2 volume fraction (Eq 24)
        eps_Li2O2_safe = pybamm.maximum(eps_Li2O2, pybamm.Scalar(0)) # keep it positive
        l_Li2O2 = ((eps_Li2O2_safe + eps_s0) / eps_s0) ** (pybamm.Scalar(1)/3) * r0 - r0 # Eq 23
        a = a0 * (1 - pybamm.Scalar(0.5) * pybamm.erf((l_Li2O2 - lm) / l_sigma)) # modified version of eq 22
        #a = a0
        #a = a0 * (1 - pybamm.erf(l_Li2O2 - lm)) / 2                 # Eq 22 (tunneling effect model) 
        cLi2O2_s = pybamm.minimum(cLi2O2, cLi2O2_max) # cap peroxide concentration at max for use in Eq 15 and 19/20

        # Electrochemical Reaction Kinetics
        #phi_film = j_c * R_film * eps_Li2O2         # Eq 17
        phi_film = pybamm.Scalar(0)
        eta_c = phi_s - phi_l - phi_film - Eeq      # Eq 16

        arg_a = (1-beta) * 2*F*eta_c/(R*T)
        arg_c = -beta      * 2*F*eta_c/(R*T)

        arg_a = pybamm.maximum(pybamm.minimum(arg_a, 25), -25)
        arg_c = pybamm.maximum(pybamm.minimum(arg_c, 25), -25)

        j_c = 2*F*(k_a * cLi2O2_s * pybamm.exp(arg_a) - k_c * cLi**2 * cO2 * pybamm.exp(arg_c)) # Eq 15

        # Conservation of Charge
        i_l = -kappa_eff * pybamm.grad(phi_l)       # Simplified Eq 9 (add full version when lithium transport is modeled)
        i_s = -sigma_eff * pybamm.grad(phi_s)

        # Solve Eq 7 and 8 algebraically
        self.algebraic = {
            phi_s: pybamm.div(i_s) + a * j_c,
            phi_l: pybamm.div(i_l) - a * j_c,
        }

        # Li2O2 Product Concentration (Eq 19)
        delta = pybamm.Scalar(1e-3) * cLi2O2_max
        H = 0.5*(1 + pybamm.tanh((cLi2O2 - cLi2O2_max)/delta))
        dcLi2O2_dt = (1-H) * (a * j_c / (2*F))

        # 1D Electrode Clogging (Eq 20)
        deps_dt = H * (-a * j_c * M_Li2O2 / (2 * F * rho_Li2O2))

        # Lithium ion migration (Eq 4)
        N_Li = -D_Li_eff * pybamm.grad(cLi) + i_l * tP / F
        R_Li = -a*j_c/F
        dcLi_dt = ((-pybamm.div(N_Li) + R_Li) - cLi * deps_dt )/ eps_safe

        # 1D O2 Cathode Diffusion (Eq 5)
        N_O2 = -D_O2_eff * pybamm.grad(cO2)     # O2 flux
        dcO2_dt = ((-pybamm.div(N_O2) - a * j_c / (2*F)) - cO2 * deps_dt) / eps_safe

        self.rhs = {cLi: dcLi_dt, cLi2O2: dcLi2O2_dt, cO2: dcO2_dt, eps: deps_dt}

        eps0 = pybamm.Parameter("eps0")
        eps_Li2O2_init = 1 - eps0 - eps_s0
        l_Li2O2_init = ((eps_Li2O2_init + eps_s0) / eps_s0) ** (pybamm.Scalar(1)/3) * r0 - r0
        a_init = a0 * (1 - pybamm.Scalar(0.5) * pybamm.erf((l_Li2O2_init - lm) / l_sigma))
        a_init = pybamm.maximum(a_init, pybamm.Scalar(1e-12) * a0)

        cLi0 = pybamm.Parameter("c_Li_0")
        cO20 = cO2_ext * S_O2

        cP0 = pybamm.minimum(eps_Li2O2_init * rho_Li2O2 / M_Li2O2, cLi2O2_max)
        cP0_reg = cP0 + pybamm.Scalar(1e-12) * cLi2O2_max
        eta0_raw = (R*T/(2*F)) * pybamm.log((k_c * cLi0**2 * cO20) / (k_a * cP0_reg))
        gate = cP0 / (cP0 + pybamm.Scalar(1e-12) * cLi2O2_max)
        eta0 = gate * eta0_raw

        kappa_eff0 = kappa * pybamm.Parameter("eps0") ** b
        phi_l_guess = -J * (x - Lc) / kappa_eff0   # so phi_l(L)=0 and dphi/dx = -J/kappa_eff0
        phi_s_guess = phi_l_guess + Eeq 

        # ICs
        
        c_seed = 1e-10 * cLi2O2_max
        #c_seed = 0
        self.initial_conditions = {
            cLi: pybamm.Parameter("c_Li_0"),
            cLi2O2: c_seed,         # seed a small amount of dissolved cLi2O2
            cO2: cO2_ext * S_O2, 
            eps: pybamm.Parameter("eps0"),
            phi_s: phi_s_guess, 
            phi_l: phi_l_guess, 
            
        } 

        # BCs
        sigma_eff_L = pybamm.boundary_value(sigma_eff, "right")
        kappa_eff_0 = pybamm.boundary_value(kappa_eff, "left")

        self.boundary_conditions = {
            cLi: {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(0), "Neumann")
            },

            cO2: {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (cO2_ext * S_O2, "Dirichlet")
            },

            phi_s: {
                "left":  (pybamm.Scalar(0), "Neumann"),
                "right": (-J / sigma_eff_L, "Neumann") 
            },
            phi_l: {
                "left":  (-J / kappa_eff_0, "Neumann"),  
                "right": (pybamm.Scalar(0), "Neumann") 
            },
        
        }
        # set gauge condition
        self.boundary_conditions[phi_l]["right"] = (pybamm.Scalar(0), "Dirichlet")

        # Outputs
        self.variables = {
            "O2 concentration": cO2,
            "Porosity": eps,
            "Solid potential": phi_s,
            "Electrolyte potential": phi_l,
            "Cathode interfacial current density": j_c,
            "Cathode overpotential": eta_c,
            "Film drop [V]": phi_film,
            "Active area": a,
        }


