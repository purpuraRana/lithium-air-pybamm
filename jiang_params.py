import pybamm

def jiang_2020_params():
    return pybamm.ParameterValues({
        # Constants
        "F": 96485,  # C/mol
        "T": 300,    # K (Table 1)
        "R": 8.314,     # J K-1 mol-1
        "tau_ramp": 1,   # seconds
        "J_floor": 1e-8,
        "eps_Li2O2_seed": 0, #1e-8,

        # Separator
        "Lsep": 2.5e-5,         # m
        "eps_sep": 0.5,

        # 0D structural / deposition parameters
        "eps0": 0.75,         # cathode porosity ε0 (original .75)
        "eps_s0": 0.25,        # <-- NOT given explicitly in Table 1; assumed 0.25
        "a0": 3.67e7,         # m^2/m^3
        "r0": 25e-9,          # m
        "lm": 7e-9,           # m (tunneling mean thickness)
        "l_sigma": 2e-9,      # 2 nm (distribution width)
        "Bruggeman exponent": 1.5,

        "sigma": 10,       # S/m (positive electrode conductivity)
        "kappa": 0.5,       # S/m (Li conductivity in electrolyte)
        "D_Li": 2.11e-9,     # m2/s
        "tP": 0.2594,         # lithium transference number
        "c_Li_0": 1000,      # mol m-3 (initial electrolyte concentration)
        "c_Li2O2_max": 0.09,   # mol m-3 (peroxide solubility limit)
        
        "k_a": 1.11e-15,     # anode reaction rate
        "k_c": 3.4e-20,      # cathode reaction rate

        # O2 boundary / solubility
        "S_O2": 0.38,         # -
        "cO2_ext": 9.46,      # mol/m^3
        "cO2_0": 3.5948,         # mol/m3 (initial O2 conc)
        "Lc": 8e-4,           # m (thickness of porous cathode)
        "D_O2": 1e-9,            # m2 s-1

        # Electrochem / voltage pieces
        "Eeq": 2.96,          # V
        "R_film": 50.0,        # Ohm*m^2

        # Li2O2 properties
        "M_Li2O2": 45.88e-3,  # kg/mol
        "rho_Li2O2": 2140.0,  # kg/m^3
        "p": 0.5,              # geometrical factor typically set to 0.5
        "beta": 0.5,

        # Simulation controls
        "V_cut": 2.4,         # V (they use 2.4 V cutoff in sensitivity section)
        "k_in": 1e3,          # 1/s "clamp strength" for 0D relaxation (your choice)
        "eta_c_placeholder": 0.0,  # start as 0 if you haven't implemented BV yet
        "J_ref": 1.0,         # A/m^2 just to scale placeholder

        # Other things
        "eps_min": 1e-3,    # clamp eps_Li2O2
        "R_ohm": 0.02,       # Ohm*m^2 (add later)
        "k_film": 0.0,      # tunable film penalty for 0D
        "S_O2_vol": 1e-4,
    })
