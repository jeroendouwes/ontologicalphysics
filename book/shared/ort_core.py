"""
ORT Core — Python port of all Java physics classes.

This module contains the complete physics engine for the ORT,
ported from 15 Java source files into a single Python module.

Classes:
    ORT — Core velocity model (v²_r + v²_t = c²)
    NewtonModel       — Classical Newtonian mechanics
    EinsteinModel     — Einstein's Special Relativity
    GravityModel      — Gravity via c_local(r), Schwarzschild, RN, Kerr
    CosmologicalModel — Cosmological shell, Hubble radius

All constants are defined at module level for easy access.
"""

import math
import numpy as np

# =============================================================================
# Physical constants
# =============================================================================

C = 299_792_458.0           # Speed of light [m/s]
G = 6.67430e-11             # Gravitational constant [m³ kg⁻¹ s⁻²]
K_E = 8.9875517873681764e9  # Coulomb constant [N·m²/C²]

# Astronomical masses [kg]
M_EARTH = 5.972e24
M_SUN = 1.989e30
M_JUPITER = 1.898e27

# Astronomical radii [m]
R_EARTH = 6.371e6
R_SUN = 6.957e8
R_JUPITER = 6.9911e7

# Particle constants
Q_ELEMENTARY = 1.602176634e-19   # Elementary charge [C]
M_PROTON = 1.67262192369e-27     # Proton mass [kg]
M_ELECTRON = 9.1093837015e-31    # Electron mass [kg]

# Quantum constants
H_PLANCK = 6.62607015e-34        # Planck constant [J·s]
HBAR = H_PLANCK / (2 * math.pi)  # Reduced Planck constant [J·s]

# Planetary orbital parameters
A_MERCURY = 5.7909050e10         # Semi-major axis [m]
E_MERCURY = 0.20563              # Eccentricity
T_MERCURY = 87.969 * 86400       # Orbital period [s]
A_VENUS = 1.08208e11
E_VENUS = 0.006772
A_EARTH_ORBIT = 1.49598023e11
E_EARTH_ORBIT = 0.01671
A_MARS = 2.27939200e11
E_MARS = 0.0934
R_SATURN_ORBIT = 1.4335e12
PARSEC = 3.0857e16              # [m]

# Black holes and distances
M_SGR_A_STAR = 4.0e6 * M_SUN
D_SGR_A_STAR = 8.178e3 * PARSEC
M_M87_STAR = 6.5e9 * M_SUN
D_M87_STAR = 16.8e6 * PARSEC

# Cosmological constants
H_0 = 67.4e3 / (1e6 * PARSEC)   # Hubble constant [s⁻¹]
M_OBSERVABLE = 9.2e52            # Observable universe mass [kg]
R_HUBBLE = C / H_0              # Hubble radius [m]
RHO_CRIT = 3 * H_0**2 / (8 * math.pi * G)  # Critical density [kg/m³]
OMEGA_MATTER = 0.315
OMEGA_LAMBDA = 0.685
T_UNIVERSE = 13.8e9 * 365.25 * 86400  # Age of universe [s]

# Gravitational wave constants — GW150914
GW150914_M1 = 36.0 * M_SUN
GW150914_M2 = 29.0 * M_SUN
GW150914_DISTANCE = 1.3e25       # ~410 Mpc [m]
GW150914_STRAIN = 1.0e-21
GW150914_M_FINAL = 62.3 * M_SUN
GW150914_A_FINAL = 0.67
GW150914_E_RAD = 3.0 * M_SUN * C**2

# Hulse-Taylor pulsar
HT_M1 = 1.4408 * M_SUN
HT_M2 = 1.3886 * M_SUN
HT_PERIOD = 7.752 * 3600         # [s]
HT_PDOT_OBS = -2.4211e-12        # Observed period decay [s/s]
HT_PDOT_GRT = -2.4025e-12        # GRT prediction [s/s]
HT_ECCENTRICITY = 0.6171334

# NR fitting coefficients
NR_SPIN_COEFF_A = 2.0 * math.sqrt(3.0)  # Rezzolla 2008
NR_SPIN_COEFF_B = -3.871
NR_SPIN_COEFF_C = 4.028
NR_ERAD_COEFF = 0.194             # Healy 2014
NR_LPEAK_COEFF = 0.0164           # Healy/Hofmann

# QNM coefficients (Berti, Cardoso & Will 2006)
QNM_F1 = 1.5251
QNM_F2 = -1.1568
QNM_F3 = 0.1292
QNM_Q1 = 0.7000
QNM_Q2 = 1.4187
QNM_Q3 = 0.4990

# S2 star at Sgr A*
A_S2 = 1.534e14                  # Semi-major axis [m] (~1025 AU)
E_S2 = 0.8843                   # Eccentricity
P_S2 = 16.05 * 365.25 * 86400   # Orbital period [s]

# Kerr / Frame-dragging
J_EARTH = 5.86e33               # Angular momentum [kg·m²/s]
J_SUN = 1.92e41
A_STAR_EARTH = J_EARTH * C / (G * M_EARTH**2)
GPB_ORBIT_RADIUS = R_EARTH + 642e3  # GP-B orbit [m]
GPB_FD_MEASURED = 37.2           # GP-B frame-dragging [mas/yr]
GPB_FD_PREDICTED = 39.2
A_LAGEOS = 1.227e7              # LAGEOS semi-major axis [m]

# Orbital radii
R_ISS = R_EARTH + 408_000
R_GPS = R_EARTH + 20_200_000
R_GEO = 42_164_000
R_MOON = 384_400_000


# =============================================================================
# ORT — Core velocity model
# =============================================================================

class ORT:
    """
    Core ORT: every object moves at c through spacetime.
    v_space² + v_time² = c²

    The state is fully determined by the angle θ from the time axis.
    """

    def __init__(self, theta):
        """Create from angle θ (radians, 0 ≤ θ ≤ π/2)."""
        if not (0 <= theta <= math.pi / 2 + 1e-10):
            raise ValueError(f"θ must be in [0, π/2], got {theta}")
        self._theta = min(theta, math.pi / 2)

    @classmethod
    def from_v(cls, v_space):
        """Create from spatial velocity v [m/s]."""
        v_space = min(abs(v_space), C)
        return cls(math.asin(v_space / C))

    @classmethod
    def from_beta(cls, beta):
        """Create from β = v/c."""
        beta = min(abs(beta), 1.0)
        return cls(math.asin(beta))

    @classmethod
    def from_angle(cls, theta):
        """Create from angle θ (radians)."""
        return cls(theta)

    @property
    def theta(self):
        """Angle from time axis [rad]."""
        return self._theta

    @property
    def theta_degrees(self):
        """Angle from time axis [degrees]."""
        return math.degrees(self._theta)

    @property
    def v_space(self):
        """Spatial velocity [m/s] — formula (2a)."""
        return C * math.sin(self._theta)

    @property
    def v_time(self):
        """Time velocity [m/s] — formula (2b)."""
        return C * math.cos(self._theta)

    @property
    def beta(self):
        """Normalized velocity v/c — formula (3)."""
        return math.sin(self._theta)

    @property
    def gamma(self):
        """Lorentz factor — formula (6)."""
        ct = math.cos(self._theta)
        if ct < 1e-15:
            return float('inf')
        return 1.0 / ct

    @property
    def time_dilation(self):
        """Time dilation factor 1/γ = cos(θ) — formula (5)."""
        return math.cos(self._theta)

    def proper_time(self, t_rest):
        """Proper time elapsed for coordinate time t_rest."""
        return t_rest * math.cos(self._theta)

    def relativistic_mass(self, m0):
        """Relativistic mass γm₀ — formula (14)."""
        return self.gamma * m0

    def total_energy(self, m0):
        """Total energy E = γm₀c² — formula (7)."""
        return self.gamma * m0 * C**2

    def kinetic_energy(self, m0):
        """Kinetic energy E_k = (γ-1)m₀c² — formula (8)."""
        return (self.gamma - 1) * m0 * C**2

    def kinetic_energy_integral(self, m0, n_steps=1_000_000):
        """
        Kinetic energy via numerical work integral ∫v dp (§4.4).

        Integrates from θ=0 to θ=self.theta using the trapezoid rule.
        This is the non-circular derivation of E=mc²:
          v = c·sin(θ),  dp = m₀c/cos²(θ) dθ
          v·dp = m₀c² · sin(θ)/cos²(θ) dθ
        Result: E_kin = m₀c²(1/cos(θ) - 1) = (γ-1)m₀c²
        """
        if self._theta < 1e-15:
            return 0.0
        dtheta = self._theta / n_steps
        total = 0.0
        for i in range(n_steps):
            th = (i + 0.5) * dtheta  # midpoint rule
            cos_th = math.cos(th)
            integrand = math.sin(th) / (cos_th * cos_th)
            total += integrand * dtheta
        return m0 * C**2 * total

    def momentum(self, m0):
        """Momentum p = γm₀v — formula (9)."""
        return self.gamma * m0 * self.v_space

    def zijn_mass(self, m0):
        """Mass-Being S_m = m_rel × v_t = m₀c — formula (12)."""
        return m0 * C

    def zijn_energy(self, m0):
        """Energy-Being S_E = E × v_t = m₀c³ — derived consequence of §4.4."""
        return m0 * C**3

    def zijn_vector(self, m0):
        """
        Being-vector S⃗ = (S_space, S_time) = (p, m₀c) — formulas (16)-(17).
        Returns (S_space, S_time, |S|).
        """
        s_space = self.momentum(m0)       # = p
        s_time = m0 * C                   # = m₀c (invariant)
        s_mag = self.gamma * m0 * C       # = γm₀c = E/c — formula (18)
        return (s_space, s_time, s_mag)

    @staticmethod
    def velocity_addition(v1, v2):
        """
        Relativistic velocity addition — formula (10).
        Returns v_total = (v1 + v2) / (1 + v1·v2/c²).
        """
        return (v1 + v2) / (1 + v1 * v2 / C**2)

    @staticmethod
    def rapidity(v):
        """Rapidity φ = atanh(v/c) — formula (11)."""
        beta = v / C
        if abs(beta) >= 1.0:
            return float('inf') if beta > 0 else float('-inf')
        return math.atanh(beta)

    def lorentz_transform(self, x, t):
        """
        Lorentz transformation — formulas (20)-(22).
        Returns (x', t') for event (x, t).
        """
        v = self.v_space
        g = self.gamma
        x_prime = g * (x - v * t)
        t_prime = g * (t - v * x / C**2)
        return (x_prime, t_prime)

    def length_contraction(self, L0):
        """Contracted length L = L₀/γ = L₀·cos(θ)."""
        return L0 * math.cos(self._theta)

    # --- Quantum mechanical connections (§13) ---

    def de_broglie_wavelength(self, m0):
        """De Broglie wavelength λ = h/p = h/S_space — formula (108)."""
        p = self.momentum(m0)
        if p < 1e-50:
            return float('inf')
        return H_PLANCK / p

    @staticmethod
    def compton_wavelength(m0):
        """Compton wavelength λ_C = h/(m₀c) = h/S_time — formula (111)."""
        return H_PLANCK / (m0 * C)

    @staticmethod
    def compton_frequency(m0):
        """Compton frequency f₀ = m₀c²/h — formula (110)."""
        return m0 * C**2 / H_PLANCK

    def particle_frequency(self, m0):
        """Relativistic frequency f = E/h = γm₀c²/h."""
        return self.total_energy(m0) / H_PLANCK

    def __repr__(self):
        return (f"ORT(θ={math.degrees(self._theta):.2f}°, "
                f"β={self.beta:.6f}, γ={self.gamma:.4f})")


# =============================================================================
# NewtonModel — Classical mechanics
# =============================================================================

class NewtonModel:
    """Newtonian (classical) model for comparison."""

    @staticmethod
    def proper_time(t_rest, v):
        """Newton: absolute time, no dilation."""
        return t_rest

    @staticmethod
    def length(L0, v):
        """Newton: absolute length, no contraction."""
        return L0

    @staticmethod
    def velocity_addition(v1, v2):
        """Newton: simple addition v₁ + v₂."""
        return v1 + v2

    @staticmethod
    def momentum(m0, v):
        """Newton: p = m₀v."""
        return m0 * v

    @staticmethod
    def kinetic_energy(m0, v):
        """Newton: E_k = ½m₀v²."""
        return 0.5 * m0 * v**2


# =============================================================================
# EinsteinModel — Standard Special Relativity
# =============================================================================

class EinsteinModel:
    """Einstein's Special Relativity for comparison."""

    @staticmethod
    def gamma(v):
        """Lorentz factor γ = 1/√(1 - v²/c²)."""
        beta2 = (v / C)**2
        if beta2 >= 1.0:
            return float('inf')
        return 1.0 / math.sqrt(1 - beta2)

    @staticmethod
    def proper_time(t_rest, v):
        """Proper time τ = t/γ."""
        return t_rest / EinsteinModel.gamma(v)

    @staticmethod
    def length(L0, v):
        """Length contraction L = L₀/γ."""
        return L0 / EinsteinModel.gamma(v)

    @staticmethod
    def velocity_addition(v1, v2):
        """Relativistic velocity addition."""
        return (v1 + v2) / (1 + v1 * v2 / C**2)

    @staticmethod
    def momentum(m0, v):
        """Relativistic momentum p = γm₀v."""
        return EinsteinModel.gamma(v) * m0 * v

    @staticmethod
    def kinetic_energy(m0, v):
        """Kinetic energy E_k = (γ-1)m₀c²."""
        return (EinsteinModel.gamma(v) - 1) * m0 * C**2

    @staticmethod
    def total_energy(m0, v):
        """Total energy E = γm₀c²."""
        return EinsteinModel.gamma(v) * m0 * C**2

    @staticmethod
    def rest_energy(m0):
        """Rest energy E₀ = m₀c²."""
        return m0 * C**2

    @staticmethod
    def lorentz_x(x, t, v):
        """Lorentz transformation x' = γ(x - vt)."""
        g = EinsteinModel.gamma(v)
        return g * (x - v * t)

    @staticmethod
    def lorentz_t(x, t, v):
        """Lorentz transformation t' = γ(t - vx/c²)."""
        g = EinsteinModel.gamma(v)
        return g * (t - v * x / C**2)


# =============================================================================
# GravityModel — Gravity via c_local(r)
# =============================================================================

class GravityModel:
    """
    Gravity model based on the ORT principle: v_grav consumes part of c,
    leaving c_local² = c² - v_grav² for space and time.

    Derivation: v_grav = √(2GM/r) → c_local(r) = c·√(1 - r_s/r),
    with r_s = 2GM/c² as consequence (where v_grav = c).

    Supports Schwarzschild (mass only), Reissner-Nordström (mass + charge),
    and Kerr (mass + spin) metrics.
    """

    def __init__(self, mass, charge=0.0, spin=0.0):
        """
        Parameters:
            mass:   Mass [kg]
            charge: Electric charge [C] (for Reissner-Nordström)
            spin:   Dimensionless spin parameter a* ∈ [0, 1) (for Kerr)
        """
        self.mass = mass
        self.charge = charge
        self.spin = spin

        # Schwarzschild radius
        self.rs = 2 * G * mass / C**2

        # Reissner-Nordström charge radius²
        self.rQ2 = K_E * charge**2 * G / C**4

        # RN horizons
        discriminant = self.rs**2 - 4 * self.rQ2
        if discriminant >= 0 and charge != 0:
            sqrt_disc = math.sqrt(discriminant)
            self.r_plus = (self.rs + sqrt_disc) / 2
            self.r_minus = (self.rs - sqrt_disc) / 2
        else:
            self.r_plus = self.rs  # Schwarzschild case
            self.r_minus = 0.0

        # Kerr parameters
        self.kerr_a = spin * G * mass / C**2  # [m]
        kerr_disc = 1 - spin**2 if spin < 1 else 0
        self.kerr_r_plus = (self.rs / 2) * (1 + math.sqrt(max(0, kerr_disc)))
        self.kerr_r_minus = (self.rs / 2) * (1 - math.sqrt(max(0, kerr_disc)))

    # --- Static utility methods ---

    @staticmethod
    def schwarzschild_radius(mass):
        """r_s = 2GM/c²."""
        return 2 * G * mass / C**2

    @staticmethod
    def newton_potential(mass, r):
        """Newtonian gravitational potential Φ = -GM/r."""
        return -G * mass / r

    @staticmethod
    def hubble_recession(d):
        """Hubble recession velocity v = H₀d."""
        return H_0 * d

    @staticmethod
    def v_time_at_distance(d):
        """Time velocity at cosmological distance: v_t = √(c² - (H₀d)²)."""
        v_hub = H_0 * d
        if v_hub >= C:
            return 0.0
        return math.sqrt(C**2 - v_hub**2)

    @staticmethod
    def pericenter_distance(a, e):
        """Pericenter distance r_p = a(1-e)."""
        return a * (1 - e)

    @staticmethod
    def lens_magnification(u):
        """Gravitational lens magnification μ = (u²+2)/(u√(u²+4))."""
        return (u**2 + 2) / (u * math.sqrt(u**2 + 4))

    @staticmethod
    def chirp_mass(m1, m2):
        """Chirp mass M_c = (m1·m2)^(3/5) / (m1+m2)^(1/5)."""
        return (m1 * m2)**0.6 / (m1 + m2)**0.2

    @staticmethod
    def peters_enhancement_factor(e):
        """Peters enhancement factor f(e) for eccentric orbits."""
        return (1 + (73/24)*e**2 + (37/96)*e**4) / (1 - e**2)**3.5

    @staticmethod
    def gw_power(m1, m2, a, e=0.0):
        """Gravitational wave luminosity [W]."""
        base = (32/5) * G**4 / C**5 * (m1*m2)**2 * (m1+m2) / a**5
        if e > 0:
            base *= GravityModel.peters_enhancement_factor(e)
        return base

    @staticmethod
    def orbital_decay_rate(m1, m2, a, e=0.0):
        """Orbital decay rate da/dt [m/s]."""
        base = -(64/5) * G**3 * m1 * m2 * (m1+m2) / (C**5 * a**3)
        if e > 0:
            base *= GravityModel.peters_enhancement_factor(e)
        return base

    @staticmethod
    def time_to_merger(m1, m2, a0):
        """Time to merger for circular orbit [s]."""
        return (5/256) * C**5 * a0**4 / (G**3 * m1 * m2 * (m1+m2))

    @staticmethod
    def symmetric_mass_ratio(m1, m2):
        """Symmetric mass ratio η = m1·m2/(m1+m2)²."""
        return m1 * m2 / (m1 + m2)**2

    @staticmethod
    def final_spin(m1, m2):
        """Final spin a_f from NR fit (Rezzolla 2008)."""
        eta = GravityModel.symmetric_mass_ratio(m1, m2)
        return NR_SPIN_COEFF_A * eta + NR_SPIN_COEFF_B * eta**2 + NR_SPIN_COEFF_C * eta**3

    @staticmethod
    def radiated_energy_fraction(m1, m2):
        """Radiated energy fraction ε from NR fit (Healy 2014)."""
        eta = GravityModel.symmetric_mass_ratio(m1, m2)
        return NR_ERAD_COEFF * 4 * eta**2

    @staticmethod
    def radiated_energy(m1, m2):
        """Radiated energy E_rad = ε·(m1+m2)·c² [J]."""
        eps = GravityModel.radiated_energy_fraction(m1, m2)
        return eps * (m1 + m2) * C**2

    @staticmethod
    def final_mass(m1, m2):
        """Final mass M_f = (m1+m2)(1-ε) [kg]."""
        eps = GravityModel.radiated_energy_fraction(m1, m2)
        return (m1 + m2) * (1 - eps)

    @staticmethod
    def peak_gw_luminosity(m1, m2):
        """Peak GW luminosity L_peak [W]."""
        eta = GravityModel.symmetric_mass_ratio(m1, m2)
        return eta**2 * NR_LPEAK_COEFF * C**5 / G

    @staticmethod
    def qnm_frequency(m_final, a_final):
        """Quasi-normal mode frequency f_QNM [Hz] (Berti et al. 2006)."""
        return (C**3 / (2 * math.pi * G * m_final)) * (
            QNM_F1 + QNM_F2 * (1 - a_final)**QNM_F3)

    @staticmethod
    def qnm_quality_factor(a_final):
        """QNM quality factor Q (Berti et al. 2006)."""
        return QNM_Q1 + QNM_Q2 * (1 - a_final)**(-QNM_Q3)

    @staticmethod
    def ringdown_damping_time(m_final, a_final):
        """Ringdown damping time τ = Q/(π·f_QNM) [s]."""
        f = GravityModel.qnm_frequency(m_final, a_final)
        q = GravityModel.qnm_quality_factor(a_final)
        return q / (math.pi * f)

    @staticmethod
    def ringdown_amplitude(h_peak, m_final, a_final, t):
        """Ringdown waveform h(t) = h_peak·exp(-t/τ)·cos(2πft)."""
        f = GravityModel.qnm_frequency(m_final, a_final)
        tau = GravityModel.ringdown_damping_time(m_final, a_final)
        return h_peak * math.exp(-t / tau) * math.cos(2 * math.pi * f * t)

    @staticmethod
    def lense_thirring(J, r):
        """Lense-Thirring precession rate Ω_LT = 2GJ/(c²r³) [rad/s]."""
        return 2 * G * J / (C**2 * r**3)

    @staticmethod
    def lense_thirring_mas_per_year(J, r):
        """Lense-Thirring precession in mas/yr."""
        omega = GravityModel.lense_thirring(J, r)
        return omega * (180 / math.pi) * 3600 * 1000 * (365.25 * 86400)

    # --- Instance methods ---

    def metric_function(self, r):
        """Metric function f(r) = 1 - r_s/r + r_Q²/r² — formula (24)/(67)."""
        if r <= 0:
            return float('-inf')
        return 1 - self.rs / r + self.rQ2 / r**2

    def c_local(self, r):
        """Local speed of light c_local(r) = c·√f(r) — formula (23)."""
        f = self.metric_function(r)
        if f <= 0:
            return 0.0
        return C * math.sqrt(f)

    def time_dilation_factor(self, r):
        """Gravitational time dilation √f(r) — formula (25)."""
        f = self.metric_function(r)
        if f <= 0:
            return 0.0
        return math.sqrt(f)

    def proper_time(self, t_inf, r):
        """Proper time at radius r: τ = t·√f(r) — formula (26)."""
        return t_inf * self.time_dilation_factor(r)

    def gravitational_redshift(self, r_emit, r_obs):
        """Gravitational redshift z — formula (27)-(28)."""
        f_emit = self.metric_function(r_emit)
        f_obs = self.metric_function(r_obs)
        if f_emit <= 0 or f_obs <= 0:
            return float('inf')
        return math.sqrt(f_obs / f_emit) - 1

    def combined_time_dilation(self, r, v):
        """Combined gravitational + kinematic dilation — formula (29)-(30)."""
        cl = self.c_local(r)
        if cl <= 0:
            return 0.0
        return math.sqrt(cl**2 - v**2) / C

    def v_grav(self, r):
        """Gravitational velocity v_grav = c·√(1-f(r)) — formula (33)."""
        f = self.metric_function(r)
        if f >= 1:
            return 0.0
        return C * math.sqrt(1 - f)

    def dc_local_dr(self, r):
        """Gradient dc_local/dr = GM/(c·r²·√(1 - r_s/r)) — formula (23b)."""
        f = self.metric_function(r)
        if f <= 0:
            return float('inf')
        return G * self.mass / (C * r**2 * math.sqrt(f))

    def proper_acceleration(self, r):
        """Proper acceleration g = GM/(r²·√(1 - r_s/r)) [m/s²]."""
        f = self.metric_function(r)
        if f <= 0:
            return float('inf')
        return G * self.mass / (r**2 * math.sqrt(f))

    def spatial_stretching(self, r):
        """Spatial stretching factor 1/√f(r) — formula (38)."""
        f = self.metric_function(r)
        if f <= 0:
            return float('inf')
        return 1.0 / math.sqrt(f)

    def effective_refractive_index(self, r):
        """Effective refractive index n_eff = (c/c_local)² — formula (56)."""
        cl = self.c_local(r)
        if cl <= 0:
            return float('inf')
        return (C / cl)**2

    def light_deflection(self, b):
        """Light deflection angle [rad] — formula (41)."""
        return 2 * self.rs / b - 3 * math.pi * self.rQ2 / (4 * b**2)

    def light_deflection_arcsec(self, b):
        """Light deflection in arcseconds."""
        return self.light_deflection(b) * (180 / math.pi) * 3600

    def half_light_deflection(self, b):
        """Soldner half-deflection α = r_s/b [rad]."""
        return self.rs / b

    def shapiro_delay(self, r1, r2, b):
        """Shapiro delay [s] — formula (73)."""
        return (self.rs / C) * math.log(4 * r1 * r2 / b**2) - \
               math.pi * self.rQ2 / (C * b)

    def shapiro_delay_microseconds(self, r1, r2, b):
        """Shapiro delay in microseconds."""
        return self.shapiro_delay(r1, r2, b) * 1e6

    def shapiro_delay_roundtrip(self, r1, r2, b):
        """Shapiro round-trip delay [s]."""
        return 2 * self.shapiro_delay(r1, r2, b)

    def half_shapiro_delay(self, r1, r2, b):
        """Half Shapiro delay (temporal component only) [s]."""
        return (self.rs / (2 * C)) * math.log(4 * r1 * r2 / b**2)

    def orbital_precession(self, a, e):
        """Orbital precession per orbit [rad] — formula (49)."""
        return 3 * math.pi * self.rs / (a * (1 - e**2)) - \
               2 * math.pi * self.rQ2 / (a**2 * (1 - e**2)**2)

    def orbital_precession_arcsec(self, a, e):
        """Orbital precession in arcseconds per orbit."""
        return self.orbital_precession(a, e) * (180 / math.pi) * 3600

    def orbital_precession_arcsec_century(self, a, e, period):
        """Orbital precession in arcsec per century."""
        per_orbit = self.orbital_precession_arcsec(a, e)
        orbits_per_century = (100 * 365.25 * 86400) / period
        return per_orbit * orbits_per_century

    def orbital_period(self, r):
        """Keplerian orbital period T = 2π√(r³/(GM)) [s]."""
        return 2 * math.pi * math.sqrt(r**3 / (G * self.mass))

    def geodetic_precession(self, r):
        """Geodetic precession per orbit [rad] — formula (76)."""
        return 2 * math.pi * (1 - math.sqrt(max(0,
            1 - 1.5 * self.rs / r + 2 * self.rQ2 / r**2)))

    def geodetic_precession_per_year(self, r, period):
        """Geodetic precession in rad/yr."""
        per_orbit = self.geodetic_precession(r)
        orbits_per_year = (365.25 * 86400) / period
        return per_orbit * orbits_per_year

    def geodetic_precession_mas_per_year(self, r, period):
        """Geodetic precession in mas/yr."""
        return self.geodetic_precession_per_year(r, period) * \
               (180 / math.pi) * 3600 * 1000

    def pericenter_velocity(self, a, e):
        """Pericenter velocity [m/s]."""
        return math.sqrt(G * self.mass * (1 + e) / (a * (1 - e)))

    def einstein_ring_angle(self, d_L, d_S, d_LS):
        """Einstein ring angle θ_E [rad] — formula (81)."""
        return math.sqrt(2 * self.rs * d_LS / (d_L * d_S))

    def einstein_ring_angle_arcsec(self, d_L, d_S, d_LS):
        """Einstein ring angle in arcseconds."""
        return self.einstein_ring_angle(d_L, d_S, d_LS) * (180 / math.pi) * 3600

    def photon_sphere(self):
        """Photon sphere radius [m] — formula (79)."""
        disc = 9 * self.rs**2 - 32 * self.rQ2
        if disc < 0:
            return None
        return (3 * self.rs + math.sqrt(disc)) / 4

    def shadow_radius(self):
        """Shadow radius (critical impact parameter) [m] — formula (80)."""
        r_ph = self.photon_sphere()
        if r_ph is None:
            return None
        f = self.metric_function(r_ph)
        if f <= 0:
            return float('inf')
        return r_ph / math.sqrt(f)

    def shadow_angular_diameter(self, distance):
        """Shadow angular diameter [μas]."""
        b_crit = self.shadow_radius()
        if b_crit is None:
            return None
        return 2 * b_crit / distance * (180 / math.pi) * 3600 * 1e6

    def extremal_charge(self):
        """Extremal charge Q_ext where r+ = r- [C]."""
        return C**2 * math.sqrt(self.rs**2 / (4 * K_E * G))

    def is_extremal(self):
        """Check if black hole is extremal (r+ = r-)."""
        return abs(self.rs**2 - 4 * self.rQ2) < 1e-10 * self.rs**2

    def is_naked_singularity(self):
        """Check if parameters give a naked singularity."""
        return self.rs**2 < 4 * self.rQ2

    # --- Interior / Cosmological ---

    def c_interior(self, r):
        """Interior c for r < r_s: c_int(r) = c·√(r_s/r - 1) — formula (86)."""
        if r <= 0:
            return float('inf')
        ratio = self.rs / r - 1
        if ratio <= 0:
            return 0.0
        return C * math.sqrt(ratio)

    def c_local_gw(self, r, h0, k, omega, t):
        """c_local with GW perturbation — formula (92)."""
        c_bg = self.c_local(r) if r > self.rs else C
        return c_bg * (1 + h0 / 2 * math.sin(k * r - omega * t))

    # --- Kerr methods ---

    def kerr_sigma(self, r, theta):
        """Kerr Σ = r² + a²cos²θ — formula (97)."""
        return r**2 + self.kerr_a**2 * math.cos(theta)**2

    def kerr_delta(self, r):
        """Kerr Δ = r² - r_s·r + a² — formula (98)."""
        return r**2 - self.rs * r + self.kerr_a**2

    def kerr_ergosphere(self, theta):
        """Kerr ergosphere radius r_ergo(θ) — formula (99)."""
        return (self.rs / 2) * (1 + math.sqrt(max(0,
            1 - self.spin**2 * math.cos(theta)**2)))

    def c_local_kerr(self, r, theta):
        """Kerr c_local(r,θ) = c·√(1 - r_s·r/Σ) — formula (100)."""
        sigma = self.kerr_sigma(r, theta)
        if sigma <= 0:
            return 0.0
        val = 1 - self.rs * r / sigma
        if val <= 0:
            return 0.0
        return C * math.sqrt(val)

    def omega_frame_drag(self, r, theta):
        """Frame-dragging angular velocity ω(r,θ) [rad/s] — formula (101)."""
        a = self.kerr_a
        sigma = self.kerr_sigma(r, theta)
        delta = self.kerr_delta(r)
        denom = (r**2 + a**2)**2 - a**2 * delta * math.sin(theta)**2
        if denom <= 0:
            return 0.0
        return 2 * G * self.mass * a * r / (C * denom)

    def kerr_isco(self, prograde=True):
        """
        Kerr ISCO radius (Bardeen 1972) — formula (102)-(103).
        Returns ISCO radius [m].
        """
        a_star = self.spin
        # Intermediate variables Z1, Z2
        z1 = 1 + (1 - a_star**2)**(1/3) * (
            (1 + a_star)**(1/3) + (1 - a_star)**(1/3))
        z2 = math.sqrt(3 * a_star**2 + z1**2)

        if prograde:
            r_isco = (self.rs / 2) * (3 + z2 - math.sqrt(
                (3 - z1) * (3 + z1 + 2 * z2)))
        else:
            r_isco = (self.rs / 2) * (3 + z2 + math.sqrt(
                (3 - z1) * (3 + z1 + 2 * z2)))
        return r_isco

    # --- Effective potential ---

    def effective_potential_newton(self, r, L):
        """Newtonian effective potential."""
        return -G * self.mass / r + L**2 / (2 * r**2)

    def effective_potential(self, r, L):
        """GR effective potential with corrections."""
        newton = self.effective_potential_newton(r, L)
        gr_correction = -G * self.mass * L**2 / (r**3 * C**2)
        rn_correction = self.rQ2 * C**2 / (2 * r**2) if self.rQ2 > 0 else 0
        return newton + gr_correction + rn_correction

    def combined_redshift(self, r, v):
        """Combined gravitational + kinematic redshift z."""
        factor = self.combined_time_dilation(r, v)
        if factor <= 0:
            return float('inf')
        return 1.0 / factor - 1

    def __repr__(self):
        parts = [f"M={self.mass:.3e}kg, r_s={self.rs:.3e}m"]
        if self.charge != 0:
            parts.append(f"Q={self.charge:.3e}C")
        if self.spin != 0:
            parts.append(f"a*={self.spin:.4f}")
        return f"GravityModel({', '.join(parts)})"


# =============================================================================
# CosmologicalModel — Universe as black hole interior
# =============================================================================

class CosmologicalModel:
    """
    Cosmological extension: observable universe as BH interior.

    Key insight: R_Hubble ≈ r_s(M_observable) — the Hubble radius
    approximately equals the Schwarzschild radius of the observable mass.
    """

    def __init__(self):
        self.r_hubble = R_HUBBLE
        self.r_s_observable = GravityModel.schwarzschild_radius(M_OBSERVABLE)

    @staticmethod
    def hubble_radius():
        """Hubble radius R_H = c/H₀ [m] — formula (83)."""
        return R_HUBBLE

    @staticmethod
    def critical_density():
        """Critical density ρ_crit = 3H₀²/(8πG) [kg/m³] — formula (84)."""
        return RHO_CRIT

    @staticmethod
    def schwarzschild_equals_hubble():
        """Critical mass where r_s = R_H — formula (85)."""
        return C**3 / (2 * G * H_0)

    @staticmethod
    def c_local_cosmological(r):
        """c_local for cosmological distances using Hubble flow."""
        v_hub = H_0 * r
        if v_hub >= C:
            return 0.0
        return math.sqrt(C**2 - v_hub**2)

    @staticmethod
    def scale_factor(t, t0=None):
        """FLRW scale factor a(t) = (t/t₀)^(2/3) for matter-dominated."""
        if t0 is None:
            t0 = T_UNIVERSE
        if t <= 0:
            return 0.0
        return (t / t0)**(2/3)

    def ratio_hubble_schwarzschild(self):
        """Ratio R_H / r_s(M_observable) ≈ 1."""
        return self.r_hubble / self.r_s_observable

    def lambda_effective(self, mass):
        """Effective cosmological constant for BH of given mass."""
        rs = GravityModel.schwarzschild_radius(mass)
        if rs <= 0:
            return 0.0
        return 3 * (C / rs)**2 * OMEGA_LAMBDA

    def dark_matter_acceleration(self, r):
        """Boundary acceleration a₀·(r/R_H) for dark matter analogy."""
        a0 = C**2 / R_HUBBLE
        return a0 * (r / R_HUBBLE)


# =============================================================================
# GW helper functions (standalone, for plotting convenience)
# =============================================================================

def inspiral_separation(t, m1, m2, a0):
    """
    Inspiral orbital separation a(t) for circular orbit.
    a(t)⁴ = a₀⁴ - (256/5)·G³·m1·m2·(m1+m2)/c⁵ · t
    """
    coeff = (256/5) * G**3 * m1 * m2 * (m1 + m2) / C**5
    a4 = a0**4 - coeff * t
    if a4 <= 0:
        return 0.0
    return a4**0.25


def inspiral_gw_frequency(m1, m2, a):
    """GW frequency from orbital separation: f_GW = (1/π)·√(G(m1+m2)/a³)."""
    if a <= 0:
        return float('inf')
    return (1 / math.pi) * math.sqrt(G * (m1 + m2) / a**3)


def gw_strain(r, m_chirp, f):
    """GW strain amplitude h = (4/r)·(G·M_c/c²)^(5/3)·(πf/c)^(2/3)."""
    return (4 / r) * (G * m_chirp / C**2)**(5/3) * (math.pi * f / C)**(2/3)


# =============================================================================
# Predefined model instances for common objects
# =============================================================================

EARTH = GravityModel(M_EARTH)
SUN = GravityModel(M_SUN)
JUPITER = GravityModel(M_JUPITER)
SGR_A_STAR = GravityModel(M_SGR_A_STAR)
M87_STAR = GravityModel(M_M87_STAR)

# Kerr models
EARTH_KERR = GravityModel(M_EARTH, spin=A_STAR_EARTH)
