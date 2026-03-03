"""
ORT Plots — Shared visualization functions for all notebooks.

All function names are in English. Plot titles and axis labels switch
between Dutch and English via the `lang` parameter ('nl' or 'en').

Dependencies: matplotlib, numpy, plotly, ipywidgets
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc, FancyArrowPatch
from IPython.display import display, HTML
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:
    HAS_WIDGETS = False

from ort_core import (
    C, G, ORT, NewtonModel, EinsteinModel,
    GravityModel, CosmologicalModel, R_HUBBLE, H_0,
    M_SUN, M_EARTH, R_EARTH, R_SUN, M_SGR_A_STAR, M_M87_STAR,
    D_SGR_A_STAR, D_M87_STAR, J_EARTH, J_SUN, GPB_ORBIT_RADIUS,
    A_MERCURY, E_MERCURY, T_MERCURY, GW150914_M1, GW150914_M2,
    GW150914_DISTANCE, GW150914_STRAIN, GW150914_M_FINAL, GW150914_A_FINAL,
    inspiral_separation, inspiral_gw_frequency, gw_strain,
    A_STAR_EARTH, A_LAGEOS, PARSEC,
    H_PLANCK, HBAR, M_ELECTRON, M_PROTON
)

# =============================================================================
# Translations
# =============================================================================

_T = {
    'velocity_circle_title': {'nl': 'De Snelheidscirkel', 'en': 'The Velocity Circle'},
    'v_space': {'nl': r'$v_{ruimte}$ / c', 'en': r'$v_{space}$ / c'},
    'v_time': {'nl': r'$v_{tijd}$ / c', 'en': r'$v_{time}$ / c'},
    'at_rest': {'nl': 'In rust', 'en': 'At rest'},
    'photon': {'nl': 'Foton', 'en': 'Photon'},
    'gamma_title': {'nl': r'Lorentz-factor $\gamma$(v)', 'en': r'Lorentz Factor $\gamma$(v)'},
    'gamma_inv_title': {'nl': r'Tijddilatatie $1/\gamma$(v)', 'en': r'Time Dilation $1/\gamma$(v)'},
    'beta_label': {'nl': r'$\beta = v/c$', 'en': r'$\beta = v/c$'},
    'simultaneity_title': {'nl': 'Relativiteit van Gelijktijdigheid',
                            'en': 'Relativity of Simultaneity'},
    'space': {'nl': 'Ruimte x', 'en': 'Space x'},
    'time_axis': {'nl': 'Tijd t', 'en': 'Time t'},
    'lorentz_title': {'nl': 'Lorentz-transformatie', 'en': 'Lorentz Transformation'},
    'zijn_title': {'nl': 'De Zijn-vector', 'en': 'The Being-vector'},
    'zijn_space': {'nl': r'$S_{ruimte} = p$', 'en': r'$S_{space} = p$'},
    'zijn_time': {'nl': r'$S_{tijd} = m_0 c$', 'en': r'$S_{time} = m_0 c$'},
    'model_comparison_title': {'nl': 'Modelvergelijking', 'en': 'Model Comparison'},
    'newton': {'nl': 'Newton', 'en': 'Newton'},
    'einstein': {'nl': 'Einstein (SRT)', 'en': 'Einstein (SRT)'},
    'spacetime': {'nl': 'ORT', 'en': 'ORT'},
    'c_local_title': {'nl': r'$c_{local}(r)$ profiel', 'en': r'$c_{local}(r)$ Profile'},
    'radius': {'nl': r'Afstand $r / r_s$', 'en': r'Distance $r / r_s$'},
    'c_local_label': {'nl': r'$c_{local} / c$', 'en': r'$c_{local} / c$'},
    'embedding_title': {'nl': 'Flamm-paraboloïde (3D)', 'en': "Flamm's Paraboloid (3D)"},
    'precession_title': {'nl': 'Baanprecissie', 'en': 'Orbital Precession'},
    'deflection_title': {'nl': 'Lichtafbuiging', 'en': 'Light Deflection'},
    'photon_sphere_title': {'nl': 'Fotonsfeer en BH-schaduw',
                             'en': 'Photon Sphere and BH Shadow'},
    'einstein_ring_title': {'nl': 'Einstein-ring', 'en': 'Einstein Ring'},
    'gw_title': {'nl': 'Zwaartekrachtsgolfsignaal', 'en': 'Gravitational Wave Signal'},
    'gw_strain': {'nl': 'Strain $h$', 'en': 'Strain $h$'},
    'time_label': {'nl': 'Tijd [s]', 'en': 'Time [s]'},
    'kerr_title': {'nl': 'Kerr-geometrie', 'en': 'Kerr Geometry'},
    'frame_drag_title': {'nl': 'Frame-dragging veld', 'en': 'Frame-Dragging Field'},
    'isco_title': {'nl': 'ISCO vs spin', 'en': 'ISCO vs Spin'},
    'cosmo_title': {'nl': 'Kosmologische schil', 'en': 'Cosmological Shell'},
    'mass_label': {'nl': 'Massa', 'en': 'Mass'},
    'energy_label': {'nl': 'Energie [J]', 'en': 'Energy [J]'},
    'angle_label': {'nl': r'Hoek $\theta$ [°]', 'en': r'Angle $\theta$ [°]'},
    'prograde': {'nl': 'Prograde', 'en': 'Prograde'},
    'retrograde': {'nl': 'Retrograde', 'en': 'Retrograde'},
    'schwarzschild': {'nl': 'Schwarzschild', 'en': 'Schwarzschild'},
    'event_horizon': {'nl': 'Eventhorizon', 'en': 'Event Horizon'},
    'ergosphere': {'nl': 'Ergosfeer', 'en': 'Ergosphere'},
    'outer_horizon': {'nl': 'Buitenste horizon', 'en': 'Outer Horizon'},
    'inner_horizon': {'nl': 'Binnenste horizon', 'en': 'Inner Horizon'},
    'frequency': {'nl': 'Frequentie [Hz]', 'en': 'Frequency [Hz]'},
    'comparison_table_title': {'nl': 'Eerlijkheidstabel: zes modellen, tien effecten',
                                'en': 'Honesty Table: Six Models, Ten Effects'},
    'zijn_qm_title': {'nl': 'Zijn-vector en QM-conjugaatparen',
                       'en': 'Being-vector and QM Conjugate Pairs'},
    'de_broglie_title': {'nl': 'De Broglie-golflengte vs snelheid',
                          'en': 'De Broglie Wavelength vs Velocity'},
    'wavelength_label': {'nl': r'Golflengte $\lambda$ [m]', 'en': r'Wavelength $\lambda$ [m]'},
    'compton_label': {'nl': r'Comptonlengte $\lambda_C$', 'en': r'Compton wavelength $\lambda_C$'},
}


def _t(key, lang='nl'):
    """Get translated string."""
    return _T.get(key, {}).get(lang, key)


# Style defaults
_COLORS = {
    'newton': '#e74c3c',
    'einstein': '#3498db',
    'spacetime': '#2ecc71',
    'rn': '#9b59b6',
    'kerr': '#e67e22',
    'light': '#f1c40f',
}


def _setup_plot(figsize=(10, 6)):
    """Create a figure with consistent styling."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(True, alpha=0.3)
    return fig, ax


# =============================================================================
# 1. Velocity Circle
# =============================================================================

def velocity_circle(betas=None, lang='nl', figsize=(8, 8)):
    """
    Plot the velocity circle with vectors for given velocities.
    Formula (1): v_r² + v_t² = c²
    """
    if betas is None:
        betas = [0.0, 0.3, 0.5, 0.8, 0.9, 0.99, 1.0]

    fig, ax = plt.subplots(figsize=figsize)

    # Draw quarter circle
    theta_range = np.linspace(0, np.pi/2, 200)
    ax.plot(np.sin(theta_range), np.cos(theta_range), 'k-', linewidth=2)

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(betas)))
    for beta, color in zip(betas, colors):
        stv = ORT.from_beta(beta)
        vr = stv.beta
        vt = stv.v_time / C

        ax.annotate('', xy=(vr, vt), xytext=(0, 0),
                     arrowprops=dict(arrowstyle='->', color=color, lw=2))
        ax.plot(vr, vt, 'o', color=color, markersize=8)

        label = f'$\\beta$={beta:.2f}' if beta < 1 else _t('photon', lang)
        if beta == 0:
            label = _t('at_rest', lang)
        offset = (5, 5) if beta < 0.95 else (-40, 5)
        ax.annotate(label, (vr, vt), textcoords="offset points",
                    xytext=offset, fontsize=9, color=color)

    ax.set_xlim(-0.05, 1.15)
    ax.set_ylim(-0.05, 1.15)
    ax.set_aspect('equal')
    ax.set_xlabel(_t('v_space', lang), fontsize=13)
    ax.set_ylabel(_t('v_time', lang), fontsize=13)
    ax.set_title(_t('velocity_circle_title', lang), fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


def velocity_circle_interactive(lang='nl'):
    """Interactive velocity circle with slider for β."""
    if not HAS_WIDGETS:
        print("ipywidgets not available")
        return

    @widgets.interact(beta=widgets.FloatSlider(min=0, max=0.999, step=0.001,
                                                value=0.5, description='β = v/c',
                                                readout_format='.3f'))
    def _update(beta):
        plt.close('all')
        fig = velocity_circle([0.0, beta, 1.0], lang=lang, figsize=(6, 6))
        stv = ORT.from_beta(beta)
        print(f"β = {beta:.3f}  |  θ = {math.degrees(stv.theta):.2f}°  |  "
              f"γ = {stv.gamma:.4f}  |  v_time/c = {stv.v_time/C:.6f}")
        plt.show()


# =============================================================================
# 2. Gamma Curve
# =============================================================================

def gamma_curve(lang='nl', figsize=(10, 5)):
    """Plot γ(v) and 1/γ(v) — formulas (5)-(6)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    betas = np.linspace(0, 0.999, 1000)
    gammas = 1 / np.sqrt(1 - betas**2)

    ax1.plot(betas, gammas, color=_COLORS['spacetime'], linewidth=2)
    ax1.set_xlabel(_t('beta_label', lang), fontsize=12)
    ax1.set_ylabel(r'$\gamma$', fontsize=14)
    ax1.set_title(_t('gamma_title', lang), fontsize=13, fontweight='bold')
    ax1.set_ylim(0, 15)
    ax1.grid(True, alpha=0.3)

    # Mark key points
    for b, label in [(0.5, r'$\gamma$=1.15'), (0.8, r'$\gamma$=1.67'), (0.9, r'$\gamma$=2.29'), (0.99, r'$\gamma$=7.09')]:
        g = 1/math.sqrt(1-b**2)
        ax1.plot(b, g, 'ko', markersize=5)
        ax1.annotate(label, (b, g), textcoords="offset points", xytext=(10, 5), fontsize=8)

    ax2.plot(betas, 1/gammas, color=_COLORS['einstein'], linewidth=2)
    ax2.set_xlabel(_t('beta_label', lang), fontsize=12)
    ax2.set_ylabel(r'$1/\gamma = v_{time}/c$', fontsize=14)
    ax2.set_title(_t('gamma_inv_title', lang), fontsize=13, fontweight='bold')
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 3. Simultaneity Diagram
# =============================================================================

def simultaneity_diagram(beta=0.6, lang='nl', figsize=(8, 8)):
    """Spacetime diagram showing relativity of simultaneity."""
    fig, ax = plt.subplots(figsize=figsize)

    # Light cone
    ax.plot([-2, 0, 2], [2, 0, 2], 'y-', linewidth=1, alpha=0.5, label='Light')
    ax.fill_between([-2, 0, 2], [2, 0, 2], 3, alpha=0.05, color='yellow')

    # Rest frame axes
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)

    # Moving frame axes (tilted)
    gamma = 1 / math.sqrt(1 - beta**2)
    alpha = math.atan(beta)

    # x' axis (simultaneity line): tilted by alpha from x-axis
    x_range = np.linspace(-2, 2, 100)
    ax.plot(x_range, beta * x_range, 'b-', linewidth=2,
            label=f"x' ({_t('simultaneity_title', lang).split()[0]} $\\beta$={beta})")

    # t' axis (worldline): tilted by alpha from t-axis
    t_range = np.linspace(-0.5, 2, 100)
    ax.plot(beta * t_range, t_range, 'r-', linewidth=2, label=f"t' ($\\beta$={beta})")

    # 45° light line
    ax.plot([0, 2], [0, 2], '--', color=_COLORS['light'], linewidth=1.5, label='c')

    ax.set_xlim(-2, 2)
    ax.set_ylim(-0.5, 2.5)
    ax.set_xlabel(_t('space', lang), fontsize=12)
    ax.set_ylabel(_t('time_axis', lang), fontsize=12)
    ax.set_title(_t('simultaneity_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 4. Lorentz Transform Grid
# =============================================================================

def lorentz_transform_grid(beta=0.5, lang='nl', figsize=(10, 5)):
    """Show how a grid transforms under Lorentz transformation."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    gamma = 1 / math.sqrt(1 - beta**2)

    # Original grid
    for x in np.arange(-2, 2.5, 0.5):
        ax1.axvline(x, color='lightblue', linewidth=0.5)
        ax2.plot([gamma*(x - beta*t) for t in np.linspace(-2, 2, 50)],
                 [gamma*(t - beta*x/1) for t in np.linspace(-2, 2, 50)],
                 'b-', linewidth=0.5, alpha=0.5)
    for t in np.arange(-2, 2.5, 0.5):
        ax1.axhline(t, color='lightyellow', linewidth=0.5)
        ax2.plot([gamma*(x - beta*t) for x in np.linspace(-2, 2, 50)],
                 [gamma*(t - beta*x/1) for x in np.linspace(-2, 2, 50)],
                 'r-', linewidth=0.5, alpha=0.5)

    ax1.set_xlim(-2, 2); ax1.set_ylim(-2, 2)
    ax1.set_xlabel('x/c'); ax1.set_ylabel('t')
    ax1.set_title('S (rust)' if lang == 'nl' else 'S (rest)', fontsize=12)
    ax1.set_aspect('equal'); ax1.grid(True, alpha=0.2)

    ax2.set_xlim(-3, 3); ax2.set_ylim(-3, 3)
    ax2.set_xlabel("x'/c"); ax2.set_ylabel("t'")
    ax2.set_title(f"S' ($\\beta$={beta})", fontsize=12)
    ax2.set_aspect('equal'); ax2.grid(True, alpha=0.2)

    fig.suptitle(_t('lorentz_title', lang) + f" ($\\beta$ = {beta})",
                 fontsize=14, fontweight='bold')
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 5. Zijn-vector Diagram
# =============================================================================

def zijn_vector_diagram(betas=None, m0=1.0, lang='nl', figsize=(8, 8)):
    """Plot Being-vectors in (S_space, S_time) space — formulas (16)-(19)."""
    if betas is None:
        betas = [0.0, 0.3, 0.5, 0.8, 0.9, 0.95]

    fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(betas)))
    max_val = 0

    for beta, color in zip(betas, colors):
        stv = ORT.from_beta(beta)
        s_space, s_time, s_mag = stv.zijn_vector(m0)
        # Normalize to units of m0*c
        sr = s_space / (m0 * C)
        st = s_time / (m0 * C)

        ax.annotate('', xy=(sr, st), xytext=(0, 0),
                     arrowprops=dict(arrowstyle='->', color=color, lw=2.5))
        ax.plot(sr, st, 'o', color=color, markersize=8)
        ax.annotate(f'$\\beta$={beta:.2f}, $\\gamma$={stv.gamma:.2f}',
                     (sr, st), textcoords="offset points",
                     xytext=(8, 5), fontsize=9, color=color)
        max_val = max(max_val, sr, st)

    # S_time = m0*c line (invariant)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5,
               label=r'$S_{time} = m_0 c$ (invariant)')

    ax.set_xlim(-0.2, max_val * 1.15)
    ax.set_ylim(-0.2, max_val * 1.15)
    ax.set_aspect('equal')
    ax.set_xlabel(_t('zijn_space', lang), fontsize=13)
    ax.set_ylabel(_t('zijn_time', lang), fontsize=13)
    ax.set_title(_t('zijn_title', lang), fontsize=15, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 6. Model Comparison Bar Chart
# =============================================================================

def model_comparison_bar(scenario_name, values, lang='nl', figsize=(10, 6)):
    """
    Bar chart comparing Newton, Einstein, and SpaceTime for a given scenario.
    values: dict with keys 'labels', 'newton', 'einstein', 'spacetime'
    """
    fig, ax = plt.subplots(figsize=figsize)

    labels = values['labels']
    x = np.arange(len(labels))
    width = 0.25

    ax.bar(x - width, values['newton'], width, label=_t('newton', lang),
           color=_COLORS['newton'], alpha=0.8)
    ax.bar(x, values['einstein'], width, label=_t('einstein', lang),
           color=_COLORS['einstein'], alpha=0.8)
    ax.bar(x + width, values['spacetime'], width, label=_t('spacetime', lang),
           color=_COLORS['spacetime'], alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_title(f"{_t('model_comparison_title', lang)}: {scenario_name}",
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 7. c_local Profile
# =============================================================================

def c_local_profile(models, labels=None, r_range=None, lang='nl', figsize=(10, 6)):
    """
    Plot c_local(r)/c profile for one or more GravityModels.
    """
    fig, ax = plt.subplots(figsize=figsize)

    if labels is None:
        labels = [f"M = {m.mass:.2e} kg" for m in models]

    colors = plt.cm.Set1(np.linspace(0, 0.8, len(models)))

    for model, label, color in zip(models, labels, colors):
        if r_range is None:
            r_min = model.rs * 1.01
            r_max = model.rs * 20
            r = np.linspace(r_min, r_max, 1000)
        else:
            r = r_range

        c_loc = np.array([model.c_local(ri) / C for ri in r])
        r_norm = r / model.rs

        ax.plot(r_norm, c_loc, label=label, color=color, linewidth=2)

    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='c')
    ax.axhline(0.0, color='red', linestyle='--', alpha=0.3)
    ax.axvline(1.0, color='red', linestyle=':', alpha=0.5,
               label=_t('event_horizon', lang))

    ax.set_xlabel(_t('radius', lang), fontsize=12)
    ax.set_ylabel(_t('c_local_label', lang), fontsize=12)
    ax.set_title(_t('c_local_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_xlim(0.5, 20)
    ax.set_ylim(-0.05, 1.1)
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


def c_local_profile_interactive(lang='nl'):
    """Interactive c_local profile with mass slider."""
    if not HAS_WIDGETS:
        print("ipywidgets not available")
        return

    @widgets.interact(
        log_mass=widgets.FloatSlider(min=24, max=40, step=0.5, value=30,
                                      description='log₁₀(M/kg)'))
    def _update(log_mass):
        plt.close('all')
        mass = 10**log_mass
        model = GravityModel(mass)
        fig = c_local_profile([model], [f'M = 10^{log_mass:.1f} kg'],
                              lang=lang, figsize=(8, 5))
        print(f"r_s = {model.rs:.3e} m")
        plt.show()


# =============================================================================
# 8. Spacetime Embedding (Flamm's Paraboloid)
# =============================================================================

def spacetime_embedding_3d(mass=None, lang='nl'):
    """3D Flamm's paraboloid using plotly."""
    if not HAS_PLOTLY:
        print("plotly not available")
        return

    if mass is None:
        mass = M_SUN
    rs = GravityModel.schwarzschild_radius(mass)

    r = np.linspace(rs * 1.001, rs * 10, 200)
    theta = np.linspace(0, 2 * np.pi, 100)
    R, Theta = np.meshgrid(r, theta)

    # Flamm's paraboloid: z = 2√(r_s(r - r_s))
    Z = 2 * np.sqrt(rs * (R - rs))
    X = R * np.cos(Theta) / rs
    Y = R * np.sin(Theta) / rs

    fig = go.Figure(data=[go.Surface(x=X, y=Y, z=Z/rs,
                                      colorscale='Viridis', opacity=0.8)])
    fig.update_layout(
        title=_t('embedding_title', lang),
        scene=dict(xaxis_title='x/r_s', yaxis_title='y/r_s', zaxis_title='z/r_s',
                   aspectmode='data'),
        width=700, height=600
    )
    return fig


# =============================================================================
# 9. Orbital Precession
# =============================================================================

def orbital_precession_plot(model, a, e, n_orbits=5, lang='nl', figsize=(8, 8)):
    """Plot precessing elliptical orbit."""
    fig, ax = plt.subplots(figsize=figsize)

    delta_phi = model.orbital_precession(a, e)
    # Exaggerate for visibility
    scale = max(1, 0.1 / abs(delta_phi)) if delta_phi != 0 else 1

    colors = plt.cm.coolwarm(np.linspace(0, 1, n_orbits))
    for i in range(n_orbits):
        phi_offset = i * delta_phi * scale
        phi = np.linspace(0, 2*np.pi, 1000) + phi_offset
        r_orbit = a * (1 - e**2) / (1 + e * np.cos(phi - phi_offset))
        x = r_orbit * np.cos(phi)
        y = r_orbit * np.sin(phi)
        ax.plot(x/a, y/a, color=colors[i], linewidth=1.5,
                label=f'Orbit {i+1}' if i < 3 else None)

    ax.plot(0, 0, 'ko', markersize=10)
    ax.set_aspect('equal')
    ax.set_title(f"{_t('precession_title', lang)} ($\\Delta\\varphi$ = {math.degrees(delta_phi)*3600:.4f}\")",
                 fontsize=13, fontweight='bold')
    ax.set_xlabel('x / a'); ax.set_ylabel('y / a')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 10. Light Deflection Diagram
# =============================================================================

def light_deflection_diagram(lang='nl', figsize=(10, 6)):
    """Diagram showing light deflection by a massive object."""
    fig, ax = plt.subplots(figsize=figsize)

    # Mass at center
    circle = plt.Circle((0, 0), 0.3, color='orange', alpha=0.8)
    ax.add_patch(circle)
    ax.annotate('M', (0, 0), ha='center', va='center', fontsize=14, fontweight='bold')

    # Straight path (no deflection)
    ax.plot([-4, 4], [1.5, 1.5], 'b--', alpha=0.5, linewidth=1,
            label=r'Newton ($\frac{1}{2}\alpha$)')

    # Deflected path (GR = 2× Newton)
    x = np.linspace(-4, 4, 200)
    alpha_gr = 0.4  # exaggerated
    y = 1.5 + alpha_gr * np.arctan(-x*2) / np.pi
    ax.plot(x, y, color=_COLORS['spacetime'], linewidth=2.5,
            label=r'GR / SpaceTime ($\alpha$)')

    # Impact parameter
    ax.annotate('', xy=(0, 0.3), xytext=(0, 1.5),
                arrowprops=dict(arrowstyle='<->', color='red'))
    ax.annotate('b', (0.15, 0.9), fontsize=12, color='red')

    # Deflection angle
    ax.annotate(r'$\alpha$', (3.5, 1.7), fontsize=14, color=_COLORS['spacetime'])

    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-1, 3)
    ax.set_aspect('equal')
    ax.set_title(_t('deflection_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.2)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 11. Photon Sphere and Shadow
# =============================================================================

def photon_sphere_shadow(model=None, lang='nl', figsize=(8, 8)):
    """Plot photon sphere and black hole shadow."""
    if model is None:
        model = GravityModel(10 * M_SUN)

    fig, ax = plt.subplots(figsize=figsize)
    rs = model.rs
    r_ph = model.photon_sphere()
    b_crit = model.shadow_radius()

    # Black hole
    bh = plt.Circle((0, 0), rs/rs, color='black')
    ax.add_patch(bh)

    # Photon sphere
    ps = plt.Circle((0, 0), r_ph/rs, fill=False, color='yellow',
                     linewidth=2, linestyle='--',
                     label=f'Photon sphere (r = {r_ph/rs:.2f} r_s)')
    ax.add_patch(ps)

    # Shadow (critical impact parameter)
    shadow = plt.Circle((0, 0), b_crit/rs, fill=False, color='gray',
                         linewidth=2, linestyle=':',
                         label=f'Shadow (b = {b_crit/rs:.2f} r_s)')
    ax.add_patch(shadow)

    # Event horizon label
    ax.annotate(_t('event_horizon', lang), (0.7, -0.3),
                fontsize=10, color='white')

    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_aspect('equal')
    ax.set_facecolor('#1a1a2e')
    ax.set_title(_t('photon_sphere_title', lang), fontsize=14,
                 fontweight='bold', color='white')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.1)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 12. Einstein Ring
# =============================================================================

def einstein_ring_plot(theta_E_arcsec=1.0, lang='nl', figsize=(8, 8)):
    """Plot Einstein ring geometry."""
    fig, ax = plt.subplots(figsize=figsize)

    # Einstein ring
    ring = plt.Circle((0, 0), theta_E_arcsec, fill=False, color='cyan',
                       linewidth=3, label=f'$\\theta_E$ = {theta_E_arcsec:.2f}"')
    ax.add_patch(ring)

    # Lens (central mass)
    ax.plot(0, 0, 'o', color='orange', markersize=15, zorder=5)
    ax.annotate('Lens', (0.1, -0.15), fontsize=11, color='orange')

    # Source (behind lens)
    ax.plot(0, 0, 'x', color='red', markersize=10, markeredgewidth=2, zorder=5)
    ax.annotate('Source', (0.1, 0.1), fontsize=11, color='red')

    lim = theta_E_arcsec * 2
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\theta$ [arcsec]')
    ax.set_ylabel(r'$\theta$ [arcsec]')
    ax.set_title(_t('einstein_ring_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 13. Gravitational Wave Strain
# =============================================================================

def gw_strain_plot(m1=None, m2=None, distance=None, lang='nl', figsize=(12, 5)):
    """Plot gravitational wave inspiral + ringdown signal."""
    if m1 is None: m1 = GW150914_M1
    if m2 is None: m2 = GW150914_M2
    if distance is None: distance = GW150914_DISTANCE

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    m_chirp = GravityModel.chirp_mass(m1, m2)
    m_final = GravityModel.final_mass(m1, m2)
    a_final = GravityModel.final_spin(m1, m2)

    # Inspiral phase
    a0 = 500e3  # initial separation 500 km
    t_merge = GravityModel.time_to_merger(m1, m2, a0)

    t_inspiral = np.linspace(0, t_merge * 0.999, 5000)
    h_inspiral = []
    f_inspiral = []
    for t in t_inspiral:
        a_t = inspiral_separation(t, m1, m2, a0)
        if a_t <= 0:
            break
        f = inspiral_gw_frequency(m1, m2, a_t)
        h = gw_strain(distance, m_chirp, f)
        h_inspiral.append(h * np.sin(2 * np.pi * f * t))
        f_inspiral.append(f)

    t_plot = t_inspiral[:len(h_inspiral)]
    ax1.plot(t_plot, h_inspiral, color=_COLORS['spacetime'], linewidth=0.8)
    ax1.set_xlabel(_t('time_label', lang), fontsize=11)
    ax1.set_ylabel(_t('gw_strain', lang), fontsize=11)
    ax1.set_title('Inspiral', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Ringdown phase
    f_qnm = GravityModel.qnm_frequency(m_final, a_final)
    tau = GravityModel.ringdown_damping_time(m_final, a_final)
    h_peak = max(abs(h) for h in h_inspiral) if h_inspiral else 1e-21

    t_ring = np.linspace(0, 5 * tau, 1000)
    h_ring = [GravityModel.ringdown_amplitude(h_peak, m_final, a_final, t)
              for t in t_ring]

    ax2.plot(t_ring * 1000, h_ring, color=_COLORS['kerr'], linewidth=1.5)
    ax2.set_xlabel(_t('time_label', lang).replace('[s]', '[ms]'), fontsize=11)
    ax2.set_ylabel(_t('gw_strain', lang), fontsize=11)
    ax2.set_title(f'Ringdown (f_QNM = {f_qnm:.0f} Hz)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    fig.suptitle(_t('gw_title', lang), fontsize=14, fontweight='bold')
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


def gw_inspiral_interactive(lang='nl'):
    """Interactive GW inspiral with mass sliders."""
    if not HAS_WIDGETS:
        print("ipywidgets not available")
        return

    @widgets.interact(
        m1_solar=widgets.FloatSlider(min=5, max=100, step=1, value=36,
                                      description='m₁ [M☉]'),
        m2_solar=widgets.FloatSlider(min=5, max=100, step=1, value=29,
                                      description='m₂ [M☉]'))
    def _update(m1_solar, m2_solar):
        plt.close('all')
        m1 = m1_solar * M_SUN
        m2 = m2_solar * M_SUN
        mc = GravityModel.chirp_mass(m1, m2)
        eta = GravityModel.symmetric_mass_ratio(m1, m2)
        af = GravityModel.final_spin(m1, m2)
        print(f"M_chirp = {mc/M_SUN:.1f} M☉  |  η = {eta:.4f}  |  a_f = {af:.3f}")
        fig = gw_strain_plot(m1, m2, lang=lang)
        plt.show()


# =============================================================================
# 14. Kerr Geometry
# =============================================================================

def kerr_geometry_plot(a_star=0.9, mass=10*M_SUN, lang='nl', figsize=(8, 8)):
    """Plot Kerr horizons and ergosphere in (x, z) plane."""
    fig, ax = plt.subplots(figsize=figsize)

    model = GravityModel(mass, spin=a_star)
    rs = model.rs

    # Horizons (circles)
    r_plus = model.kerr_r_plus
    r_minus = model.kerr_r_minus

    theta = np.linspace(0, 2*np.pi, 500)

    # Outer horizon
    ax.plot(r_plus/rs * np.cos(theta), r_plus/rs * np.sin(theta),
            'r-', linewidth=2, label=f"{_t('outer_horizon', lang)} (r₊ = {r_plus/rs:.3f} r_s)")

    # Inner horizon
    if r_minus > 0:
        ax.plot(r_minus/rs * np.cos(theta), r_minus/rs * np.sin(theta),
                'b--', linewidth=2,
                label=f"{_t('inner_horizon', lang)} (r₋ = {r_minus/rs:.3f} r_s)")

    # Ergosphere (θ-dependent)
    th = np.linspace(0, 2*np.pi, 500)
    r_ergo = np.array([model.kerr_ergosphere(t) for t in th])
    x_ergo = r_ergo/rs * np.cos(th)
    z_ergo = r_ergo/rs * np.sin(th)
    ax.plot(x_ergo, z_ergo, color=_COLORS['kerr'], linewidth=2, linestyle='-.',
            label=f"{_t('ergosphere', lang)}")
    ax.fill_between(x_ergo, z_ergo, np.zeros_like(z_ergo), alpha=0.1,
                     color=_COLORS['kerr'])

    # Schwarzschild reference
    ax.plot(0.5 * np.cos(theta), 0.5 * np.sin(theta),
            'k:', linewidth=1, alpha=0.5,
            label=f"{_t('schwarzschild', lang)} (r = r_s/2)")

    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x / r_s$', fontsize=12)
    ax.set_ylabel(r'$z / r_s$', fontsize=12)
    ax.set_title(f"{_t('kerr_title', lang)} (a* = {a_star})",
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 15. Kerr Frame-Dragging Field
# =============================================================================

def kerr_frame_drag_field(model=None, lang='nl'):
    """2D field plot of ω(r,θ) using plotly."""
    if not HAS_PLOTLY:
        print("plotly not available")
        return

    if model is None:
        model = GravityModel(10 * M_SUN, spin=0.9)
    rs = model.rs

    r = np.linspace(rs * 1.5, rs * 10, 100)
    th = np.linspace(0.01, np.pi - 0.01, 100)
    R, TH = np.meshgrid(r, th)

    omega = np.zeros_like(R)
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            omega[i, j] = model.omega_frame_drag(R[i, j], TH[i, j])

    X = R / rs * np.sin(TH)
    Z = R / rs * np.cos(TH)

    fig = go.Figure(data=go.Contour(
        x=X[0], y=Z[:, 0], z=omega,
        colorscale='Inferno',
        colorbar=dict(title='ω [rad/s]')
    ))
    fig.update_layout(
        title=_t('frame_drag_title', lang),
        xaxis_title='r·sin(θ) / r_s',
        yaxis_title='r·cos(θ) / r_s',
        width=700, height=600
    )
    return fig


# =============================================================================
# 16. ISCO Comparison Plot
# =============================================================================

def isco_comparison_plot(lang='nl', figsize=(10, 6)):
    """Plot ISCO radius vs spin for prograde and retrograde orbits."""
    fig, ax = plt.subplots(figsize=figsize)

    spins = np.linspace(0, 0.998, 200)
    isco_pro = []
    isco_ret = []

    for a_star in spins:
        model = GravityModel(10 * M_SUN, spin=a_star)
        rs = model.rs
        isco_pro.append(model.kerr_isco(prograde=True) / rs)
        isco_ret.append(model.kerr_isco(prograde=False) / rs)

    ax.plot(spins, isco_pro, color=_COLORS['spacetime'], linewidth=2.5,
            label=_t('prograde', lang))
    ax.plot(spins, isco_ret, color=_COLORS['newton'], linewidth=2.5,
            label=_t('retrograde', lang))
    ax.axhline(3, color='gray', linestyle='--', alpha=0.5,
               label=f'{_t("schwarzschild", lang)} (3 r_s)')

    ax.set_xlabel(r'$a_*$ (spin)', fontsize=12)
    ax.set_ylabel(r'$r_{ISCO} / r_s$', fontsize=12)
    ax.set_title(_t('isco_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 5)
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 17. Cosmological Shell Diagram
# =============================================================================

def cosmological_shell_diagram(lang='nl', figsize=(10, 6)):
    """Diagram of the cosmological shell (R_H ≈ r_s)."""
    fig, ax = plt.subplots(figsize=figsize)

    cosmo = CosmologicalModel()
    ratio = cosmo.ratio_hubble_schwarzschild()

    # Hubble sphere
    theta = np.linspace(0, 2*np.pi, 200)
    ax.plot(np.cos(theta), np.sin(theta), 'b-', linewidth=2,
            label=f'R_Hubble = c/H₀ = {R_HUBBLE:.3e} m')

    # Schwarzschild radius
    r_s_norm = 1.0 / ratio
    ax.plot(r_s_norm * np.cos(theta), r_s_norm * np.sin(theta),
            'r--', linewidth=2,
            label=f'r_s(M_obs) = {cosmo.r_s_observable:.3e} m')

    ax.annotate(f'R_H / r_s = {ratio:.3f}', (0, 0), ha='center',
                fontsize=14, fontweight='bold')

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.set_title(_t('cosmo_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 18. Comparison Table (HTML)
# =============================================================================

def comparison_table(data=None, lang='nl'):
    """
    Render a comparison table as HTML.
    data: list of dicts with keys 'effect', 'newton', 'michell', 'soldner',
          'art', 'spacetime', 'thomas'
    """
    if data is None:
        data = _default_comparison_data(lang)

    title = _t('comparison_table_title', lang)
    headers = ['Effect', 'Newton', 'Michell/Laplace', 'Soldner/Einstein 1911',
               'ART (1915)', 'SpaceTime', 'Thomas']

    html = f'<h3>{title}</h3>'
    html += '<table style="border-collapse:collapse; width:100%; font-size:14px;">'
    html += '<tr>'
    for h in headers:
        html += f'<th style="border:1px solid #ddd; padding:8px; background:#f0f0f0;">{h}</th>'
    html += '</tr>'

    for row in data:
        html += '<tr>'
        html += f'<td style="border:1px solid #ddd; padding:8px; font-weight:bold;">{row["effect"]}</td>'
        for key in ['newton', 'michell', 'soldner', 'art', 'spacetime', 'thomas']:
            val = row.get(key, '')
            color = '#2ecc71' if val == '✓' else '#e74c3c' if val == '✗' else '#f39c12' if '½' in val or '⚠' in val else '#333'
            html += f'<td style="border:1px solid #ddd; padding:8px; text-align:center; color:{color}; font-size:16px;">{val}</td>'
        html += '</tr>'
    html += '</table>'

    display(HTML(html))


def _default_comparison_data(lang='nl'):
    """Default 10-effect comparison table data."""
    effects_nl = [
        'Gravitationele tijddilatatie', 'Gravitationele roodverschuiving',
        'Lichtafbuiging', 'Baanprecissie', 'Shapiro-vertraging',
        'Geodetische precessie', 'Eventhorizon', 'BH-interieur',
        'Frame-dragging', 'Zwaartekrachtsgolven', 'Kosmologie'
    ]
    effects_en = [
        'Gravitational time dilation', 'Gravitational redshift',
        'Light deflection', 'Orbital precession', 'Shapiro delay',
        'Geodetic precession', 'Event horizon', 'BH interior',
        'Frame-dragging', 'Gravitational waves', 'Cosmology'
    ]
    effects = effects_nl if lang == 'nl' else effects_en

    rows = []
    symbols = [
        ('✗', '✗', '✗', '✓', '✓', '—'),
        ('✗', '✗', '✗', '✓', '✓', '—'),
        ('✗', '✗', '½', '✓', '✓', '—'),
        ('✗', '✗', '✗', '✓', '✓', '—'),
        ('✗', '✗', '½', '✓', '✓', '—'),
        ('✗', '✗', '✗', '✓', '✓', '⅓'),
        ('v_esc', 'r=2GM/c²', '—', '✓', '✓', '—'),
        ('—', '—', '—', '✓', '✓', '—'),
        ('✗', '✗', '✗', '✓', '✓', '—'),
        ('✗', '✗', '✗', '✓', '✓', '—'),
        ('✗', '✗', '✗', '✓', '⚠', '—'),
    ]

    for effect, (n, mi, so, ar, st, th) in zip(effects, symbols):
        rows.append({
            'effect': effect, 'newton': n, 'michell': mi, 'soldner': so,
            'art': ar, 'spacetime': st, 'thomas': th
        })
    return rows


# =============================================================================
# 19. Being-vector QM Diagram
# =============================================================================

def zijn_qm_diagram(m0=None, lang='nl', figsize=(10, 8)):
    """
    Diagram showing the Being-vector with QM conjugate pairs annotated.
    Shows S_space = p (↔ x, de Broglie) and S_time = m₀c (↔ t, Compton).
    """
    if m0 is None:
        m0 = M_ELECTRON

    fig, ax = plt.subplots(figsize=figsize)

    betas = [0.0, 0.3, 0.5, 0.8, 0.95]
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(betas)))

    lambda_C = ORT.compton_wavelength(m0)

    for beta, color in zip(betas, colors):
        stv = ORT.from_beta(beta)
        p = stv.momentum(m0)
        s_time = m0 * C

        # Normalize to m₀c
        sr = p / (m0 * C)
        st = 1.0  # s_time / (m0 * C) = 1

        ax.annotate('', xy=(sr, st), xytext=(0, 0),
                     arrowprops=dict(arrowstyle='->', color=color, lw=2.5))
        ax.plot(sr, st, 'o', color=color, markersize=8)

        # de Broglie wavelength
        lam = stv.de_broglie_wavelength(m0) if beta > 0 else float('inf')
        lam_str = f'{lam:.3e} m' if lam < 1e10 else '∞'
        ax.annotate(f'$\\beta$={beta:.2f}\n$\\lambda$={lam_str}',
                     (sr, st), textcoords="offset points",
                     xytext=(8, 8), fontsize=8, color=color)

    # S_time line
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5,
               label=r'$S_{time} = m_0 c$ (invariant)')

    # QM annotations
    ax.annotate(r'$\Delta x \cdot \Delta p \geq \hbar/2$',
                (0.6, 0.15), fontsize=12, color='navy',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))
    ax.annotate(r'$\Delta t \cdot \Delta E \geq \hbar/2$',
                (0.05, 1.15), fontsize=12, color='darkred',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    # Compton wavelength annotation
    ax.annotate(f'$\\lambda_C$ = h/$S_{{time}}$ = {lambda_C:.3e} m',
                (0.05, 0.85), fontsize=10, color='gray')

    max_sr = max(ORT.from_beta(b).momentum(m0) / (m0 * C)
                 for b in betas)
    ax.set_xlim(-0.2, max_sr * 1.3)
    ax.set_ylim(-0.2, 1.5)
    ax.set_aspect('equal')
    ax.set_xlabel(_t('zijn_space', lang), fontsize=13)
    ax.set_ylabel(_t('zijn_time', lang), fontsize=13)
    ax.set_title(_t('zijn_qm_title', lang), fontsize=15, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig


# =============================================================================
# 20. De Broglie Wavelength Plot
# =============================================================================

def de_broglie_plot(m0=None, lang='nl', figsize=(10, 6)):
    """
    Plot de Broglie wavelength λ = h/p vs velocity (β), with Compton
    wavelength as horizontal asymptote.
    """
    if m0 is None:
        m0 = M_ELECTRON

    fig, ax = plt.subplots(figsize=figsize)

    betas = np.linspace(0.001, 0.999, 1000)
    wavelengths = []
    for beta in betas:
        stv = ORT.from_beta(beta)
        wavelengths.append(stv.de_broglie_wavelength(m0))

    lambda_C = ORT.compton_wavelength(m0)

    ax.semilogy(betas, wavelengths, color=_COLORS['spacetime'], linewidth=2.5,
                label=r'$\lambda = h/p = h/S_{space}$')
    ax.axhline(lambda_C, color='red', linestyle='--', linewidth=1.5,
               label=f'{_t("compton_label", lang)} = {lambda_C:.3e} m')

    # Mark some key points
    for beta_mark in [0.01, 0.1, 0.5, 0.9, 0.99]:
        stv = ORT.from_beta(beta_mark)
        lam = stv.de_broglie_wavelength(m0)
        ax.plot(beta_mark, lam, 'ko', markersize=4)
        ax.annotate(f'{lam:.2e}', (beta_mark, lam),
                    textcoords="offset points", xytext=(5, 5), fontsize=7)

    ax.set_xlabel(_t('beta_label', lang), fontsize=12)
    ax.set_ylabel(_t('wavelength_label', lang), fontsize=12)
    ax.set_title(_t('de_broglie_title', lang), fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3, which='both')
    try:
        plt.tight_layout()
    except ValueError:
        pass
    return fig
