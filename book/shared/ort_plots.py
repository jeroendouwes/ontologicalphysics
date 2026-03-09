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
    'unit_circle_title': {
        'nl': r'Eenheidscirkel: $\cos\theta$ geeft zowel tijddilatatie als lengtecontractie',
        'en': r'Unit circle: $\cos\theta$ gives both time dilation and length contraction'},
    'unit_circle_td': {'nl': r'Tijddilatatie: $\tau = t \cos\theta$',
                        'en': r'Time dilation: $\tau = t \cos\theta$'},
    'unit_circle_lc': {'nl': r'Lengtecontractie: $L = L_0 \cos\theta$',
                        'en': r'Length contraction: $L = L_0 \cos\theta$'},
    'unit_circle_1s': {'nl': '1 seconde', 'en': '1 second'},
    'unit_circle_1ls': {'nl': r'$L_0$', 'en': r'$L_0$'},
    'unit_circle_space': {'nl': 'ruimte', 'en': 'space'},
    'unit_circle_time': {'nl': 'tijd', 'en': 'time'},
    'light_ship_title': {
        'nl': 'Licht beweegt voor iedereen met c',
        'en': 'Light moves at c for everyone'},
    'light_ship_astronaut': {'nl': 'Astronaut (in het schip)',
                              'en': 'Astronaut (in the ship)'},
    'light_ship_platform': {'nl': 'Waarnemer (op het perron)',
                             'en': 'Observer (on the platform)'},
    'light_ship_simultaneous': {'nl': 'Gelijktijdig!', 'en': 'Simultaneous!'},
    'light_ship_not_simul': {'nl': 'Niet gelijktijdig!', 'en': 'Not simultaneous!'},
    'light_ship_back_first': {'nl': 'Achterkant\eerst', 'en': 'Back\nfirst'},
    'light_ship_front_later': {'nl': 'Voorkant\nlater', 'en': 'Front\nlater'},
    'duality_title': {
        'nl': 'Gelijktijdigheid en gelijkplaatsigheid: symmetrisch rond de lichtlijn',
        'en': 'Simultaneity and co-locality: symmetric around the light line'},
    'duality_worldline': {'nl': 'Wereldlijn', 'en': 'Worldline'},
    'duality_simline': {'nl': 'Gelijktijdigheidslijn', 'en': 'Simultaneity line'},
    'duality_lightline': {'nl': 'Lichtlijn', 'en': 'Light line'},
    'duality_colocal': {'nl': 'Gelijkplaatsigheid\n(zelfde plek, andere tijd)',
                         'en': 'Co-locality\n(same place, different time)'},
    'duality_simul': {'nl': 'Gelijktijdigheid\n(zelfde tijd, andere plek)',
                       'en': 'Simultaneity\n(same time, different place)'},
    'lorentz_decomp_title': {
        'nl': r'Lorentz-transformatie: coördinaatmenging via $\theta$',
        'en': r'Lorentz transformation: coordinate mixing via $\theta$'},
    'lorentz_decomp_x': {
        'nl': r"Ruimtecoördinaat $x$: bijdragen van $x'$ en $t'$",
        'en': r"Space coordinate $x$: contributions from $x'$ and $t'$"},
    'lorentz_decomp_t': {
        'nl': r"Tijdcoördinaat $ct$: bijdragen van $x'$ en $ct'$",
        'en': r"Time coordinate $ct$: contributions from $x'$ and $ct'$"},
    'lorentz_decomp_space_part': {'nl': 'ruimtedeel', 'en': 'space part'},
    'lorentz_decomp_time_part': {'nl': 'tijddeel', 'en': 'time part'},
    'lorentz_decomp_event': {'nl': 'Gebeurtenis P', 'en': 'Event P'},
    'obs_axes_title': {
        'nl': r'Eenheidsvectoren van de bewegende waarnemer ($\theta = %d°$)',
        'en': r'Unit vectors of the moving observer ($\theta = %d°$)'},
    'obs_axes_time_dir': {'nl': 'tijdrichting', 'en': 'time direction'},
    'obs_axes_space_dir': {'nl': 'ruimterichting', 'en': 'space direction'},
    'obs_axes_td': {'nl': r'$\cos\theta$ = tijddilatatie', 'en': r'$\cos\theta$ = time dilation'},
    'obs_axes_lc': {'nl': r'$\cos\theta$ = lengtecontractie', 'en': r'$\cos\theta$ = length contraction'},
    'obs_axes_td_cross': {'nl': r'$\sin\theta$ (kruisterm)', 'en': r'$\sin\theta$ (cross-term)'},
    'obs_axes_lc_cross': {'nl': r'$-\sin\theta$ (kruisterm)', 'en': r'$-\sin\theta$ (cross-term)'},
    'obs_axes_space': {'nl': '$x$ (lichtseconden)', 'en': '$x$ (light-seconds)'},
    'obs_axes_time': {'nl': '$ct$ (lichtseconden)', 'en': '$ct$ (light-seconds)'},
    'lc_sym_title': {
        'nl': r'Lengtecontractie uit symmetrie ($\theta = %d°$)',
        'en': r'Length contraction from symmetry ($\theta = %d°$)'},
    'lc_sym_light': {'nl': 'licht ($45°$)', 'en': 'light ($45°$)'},
    'lc_sym_td': {'nl': r'tijddilatatie: $\cos\theta$', 'en': r'time dilation: $\cos\theta$'},
    'lc_sym_lc': {'nl': r'lengtecontractie: $\cos\theta$', 'en': r'length contraction: $\cos\theta$'},
    'lc_sym_time_dir': {'nl': 'tijdrichting', 'en': 'time direction'},
    'lc_sym_space_dir': {'nl': 'ruimterichting', 'en': 'space direction'},
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

        if beta == 0:
            label = _t('at_rest', lang)
        elif beta >= 1:
            label = _t('photon', lang)
        elif lang == 'nl':
            label = f'$v_{{ruimte}}$ = {beta:.2f} c'
        else:
            label = f'$v_{{space}}$ = {beta:.2f} c'
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
        msg = ("Interactieve versie — download het notebook om de slider te gebruiken."
               if lang == 'nl' else
               "Interactive version — download the notebook to use the slider.")
        print(msg)
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
# 3a. Light-in-spaceship diagram (relativity of simultaneity)
# =============================================================================

def light_ship_diagram(beta=0.5, lang='nl', figsize=(9, 9)):
    """Spacetime diagram of light in a moving ship.

    Shows worldlines of front/back walls, light rays at 45 degrees,
    and the simultaneity line connecting the two arrival events.
    Demonstrates that the simultaneity axis (x') and co-locality axis (ct')
    are symmetric around the light line.
    Uses units where c = 1.
    """
    fig, ax = plt.subplots(figsize=figsize)

    L = 1.0  # rest length of ship (in light-seconds)
    alpha = math.atan(beta)  # tilt angle (both axes)

    # Time range
    t_max = L / (1 - beta) * 0.65

    # --- Worldlines of ship walls ---
    t_arr = np.linspace(-0.15, t_max, 200)
    # Back wall: x = -L/2 + beta*t
    ax.plot(-L / 2 + beta * t_arr, t_arr, color='#888888', lw=1.5, ls='-')
    # Front wall: x = L/2 + beta*t
    ax.plot(L / 2 + beta * t_arr, t_arr, color='#888888', lw=1.5, ls='-')
    # Center (worldline = ct' axis): x = beta*t
    wl_end = t_max * 1.05
    ax.plot([0, beta * wl_end], [0, wl_end], color='#2070c0', lw=2.5)

    # --- Light rays from origin at 45 degrees ---
    ax.plot([0, t_max], [0, t_max], color='#f0c020', lw=2, ls='--', alpha=0.8)
    ax.plot([0, -t_max * 0.5], [0, t_max * 0.5], color='#f0c020', lw=2, ls='--', alpha=0.8)

    # --- Intersection points ---
    # Back wall: -t = -L/2 + beta*t  =>  t_back = L / (2(1+beta))
    t_back = L / (2 * (1 + beta))
    x_back = -t_back  # on leftward light ray
    # Front wall: t = L/2 + beta*t  =>  t_front = L / (2(1-beta))
    t_front = L / (2 * (1 - beta))
    x_front = t_front  # on rightward light ray

    ax.plot(x_back, t_back, 'o', color='#c03020', ms=10, zorder=5)
    ax.plot(x_front, t_front, 'o', color='#c03020', ms=10, zorder=5)

    # --- Dashed line connecting the two events (parallel to x'-axis) ---
    ax.plot([x_back, x_front], [t_back, t_front],
            color='#c03020', lw=1.5, ls='--', alpha=0.5, zorder=3)

    # --- Axes ---
    ax_len = t_max * 1.1
    ax.annotate('', xy=(ax_len, 0), xytext=(-ax_len * 0.4, 0),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.annotate('', xy=(0, ax_len), xytext=(0, -0.1),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.text(ax_len + 0.02, -0.06, r'$x$', fontsize=15)
    ax.text(-0.08, ax_len + 0.02, r'$ct$', fontsize=15)

    # --- x'-axis through origin (simultaneity axis) ---
    xp_len = ax_len * 0.85
    ax.plot([0, xp_len * math.cos(alpha)], [0, xp_len * math.sin(alpha)],
            color='#c03020', lw=2, ls='-', alpha=0.5)

    # --- Light line reference (extended) ---
    ax.plot([0, ax_len * 0.9], [0, ax_len * 0.9], color='#f0c020', lw=1.5,
            ls='--', alpha=0.5)

    # --- Angle arcs ---
    arc_r = 0.3
    # Worldline angle alpha from ct-axis
    a_wl = np.linspace(np.pi / 2 - alpha, np.pi / 2, 50)
    ax.plot(arc_r * np.cos(a_wl), arc_r * np.sin(a_wl), color='#2070c0', lw=2)
    mid_wl = np.pi / 2 - alpha / 2
    ax.text(arc_r * 1.35 * math.cos(mid_wl), arc_r * 1.35 * math.sin(mid_wl),
            r'$\alpha$', fontsize=15, color='#2070c0', ha='center', va='center')

    # Simultaneity line angle alpha from x-axis (arc connects to x'-axis at origin)
    a_sl = np.linspace(0, alpha, 50)
    ax.plot(arc_r * np.cos(a_sl), arc_r * np.sin(a_sl), color='#c03020', lw=2)
    mid_sl = alpha / 2
    ax.text(arc_r * 1.35 * math.cos(mid_sl), arc_r * 1.35 * math.sin(mid_sl),
            r'$\alpha$', fontsize=15, color='#c03020', ha='center', va='center')

    # --- Labels ---
    # Worldline label (at 60% along the line, inside diagram)
    wl_lbl_t = wl_end * 0.55
    ax.text(beta * wl_lbl_t + 0.06, wl_lbl_t,
            _t('duality_worldline', lang) + r" ($ct'$)",
            fontsize=11, color='#2070c0', rotation=math.degrees(math.atan2(1, beta)),
            va='bottom')

    # Simultaneity line label (on x'-axis, at 60%)
    xp_lbl = xp_len * 0.55
    ax.text(xp_lbl * math.cos(alpha), xp_lbl * math.sin(alpha) - 0.06,
            _t('duality_simline', lang) + r" ($x'$)",
            fontsize=11, color='#c03020', ha='center', va='top',
            rotation=math.degrees(alpha))

    # Light line label (at 45% along)
    ax.text(ax_len * 0.42, ax_len * 0.47, _t('duality_lightline', lang),
            fontsize=10, color='#c09800', rotation=45, ha='center',
            bbox=dict(boxstyle='round,pad=0.2', fc='#fffde0', ec='none', alpha=0.8))

    # Event labels
    if lang == 'nl':
        ax.text(x_back - 0.04, t_back + 0.03, 'Achterkant',
                fontsize=9, ha='right', color='#c03020')
        ax.text(x_front + 0.04, t_front + 0.03, 'Voorkant',
                fontsize=9, ha='left', color='#c03020')
    else:
        ax.text(x_back - 0.04, t_back + 0.03, 'Back',
                fontsize=9, ha='right', color='#c03020')
        ax.text(x_front + 0.04, t_front + 0.03, 'Front',
                fontsize=9, ha='left', color='#c03020')

    # Lamp event
    ax.plot(0, 0, 'o', color='#f0c020', ms=10, zorder=5)
    if lang == 'nl':
        ax.text(0.06, -0.04, 'Lamp aan', fontsize=9, color='#c09800', va='top')
    else:
        ax.text(0.06, -0.04, 'Lamp on', fontsize=9, color='#c09800', va='top')

    ax.set_xlim(-ax_len * 0.45, ax_len * 1.15)
    ax.set_ylim(-0.15, ax_len * 1.15)
    ax.set_aspect('equal')
    ax.set_title(_t('light_ship_title', lang) +
                 r'  ($v/c$ = %.1f,  $\alpha$ = %.1f°)' % (beta, math.degrees(alpha)),
                 fontsize=13)
    ax.grid(True, alpha=0.15)
    ax.set_xticks([])
    ax.set_yticks([])
    try:
        plt.tight_layout()
    except ValueError:
        pass
    plt.show()
    return fig


# =============================================================================
# 3b. Simultaneity–Co-locality Duality Diagram
# =============================================================================

def duality_diagram(beta=0.5, lang='nl', figsize=(8, 8)):
    """Show that simultaneity and co-locality are symmetric around the light line.

    The worldline (co-locality) tilts from the ct-axis toward the light line.
    The simultaneity line tilts from the x-axis toward the light line by the
    same angle. Both make angle arctan(beta) with their respective axis.
    """
    fig, ax = plt.subplots(figsize=figsize)
    L = 1.6  # axis length

    # Axes
    ax.annotate('', xy=(L, 0), xytext=(-0.1, 0),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.annotate('', xy=(0, L), xytext=(0, -0.1),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.text(L + 0.05, -0.08, r'$x$', fontsize=15)
    ax.text(-0.12, L + 0.03, r'$ct$', fontsize=15)

    # Light line (45 degrees)
    ax.plot([0, L], [0, L], color='#f0c020', lw=2.5, ls='--', alpha=0.8)
    ax.text(L * 0.72, L * 0.78, _t('duality_lightline', lang),
            fontsize=12, color='#c09800', rotation=45, ha='center',
            bbox=dict(boxstyle='round,pad=0.2', fc='#fffde0', ec='none', alpha=0.8))

    # Angle from beta
    alpha = math.atan(beta)  # angle worldline makes with ct-axis
    # Same angle: simultaneity line makes with x-axis

    # --- Worldline (co-locality): tilts from ct-axis ---
    wl_x = L * math.sin(alpha)
    wl_y = L * math.cos(alpha)
    ax.plot([0, wl_x], [0, wl_y], color='#2070c0', lw=2.5)
    ax.text(wl_x + 0.05, wl_y - 0.05,
            _t('duality_worldline', lang) + r" ($ct'$)",
            fontsize=12, color='#2070c0')

    # --- Simultaneity line: tilts from x-axis ---
    sl_x = L * math.cos(alpha)
    sl_y = L * math.sin(alpha)
    ax.plot([0, sl_x], [0, sl_y], color='#c03020', lw=2.5)
    ax.text(sl_x + 0.05, sl_y - 0.05,
            _t('duality_simline', lang) + r" ($x'$)",
            fontsize=12, color='#c03020')

    # --- Rest frame references (faded) ---
    ax.plot([0, 0], [0, L], color='#2070c0', lw=1, ls=':', alpha=0.4)
    ax.plot([0, L], [0, 0], color='#c03020', lw=1, ls=':', alpha=0.4)

    # --- Angle arcs ---
    # Worldline angle from ct-axis
    arc_r = 0.35
    a_wl = np.linspace(np.pi / 2 - alpha, np.pi / 2, 50)
    ax.plot(arc_r * np.cos(a_wl), arc_r * np.sin(a_wl), color='#2070c0', lw=1.5)
    mid_wl = np.pi / 2 - alpha / 2
    ax.text(arc_r * 1.25 * math.cos(mid_wl), arc_r * 1.25 * math.sin(mid_wl),
            r'$\alpha$', fontsize=14, color='#2070c0', ha='center', va='center')

    # Simultaneity line angle from x-axis
    a_sl = np.linspace(0, alpha, 50)
    ax.plot(arc_r * np.cos(a_sl), arc_r * np.sin(a_sl), color='#c03020', lw=1.5)
    mid_sl = alpha / 2
    ax.text(arc_r * 1.25 * math.cos(mid_sl), arc_r * 1.25 * math.sin(mid_sl),
            r'$\alpha$', fontsize=14, color='#c03020', ha='center', va='center')

    # --- Symmetry annotation ---
    # Double-headed arrow between the two lines near the light line
    mid_r = L * 0.45
    ax.annotate('',
                xy=(mid_r * math.cos(alpha), mid_r * math.sin(alpha)),
                xytext=(mid_r * math.sin(alpha), mid_r * math.cos(alpha)),
                arrowprops=dict(arrowstyle='<->', color='purple', lw=1.5,
                                connectionstyle='arc3,rad=0.3'))
    ax.text(mid_r * 0.62, mid_r * 0.62 + 0.12, r'$\alpha = \alpha$',
            fontsize=12, color='purple', ha='center',
            bbox=dict(boxstyle='round,pad=0.2', fc='#f0e8ff', ec='none', alpha=0.8))

    # Labels for the concepts
    ax.text(-0.15, L * 0.6, _t('duality_colocal', lang),
            fontsize=10, color='#2070c0', ha='right', va='center',
            bbox=dict(boxstyle='round,pad=0.3', fc='#e8f0ff', ec='#2070c0', alpha=0.6))
    ax.text(L * 0.6, -0.18, _t('duality_simul', lang),
            fontsize=10, color='#c03020', ha='center', va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='#ffe8e8', ec='#c03020', alpha=0.6))

    ax.set_xlim(-0.4, L + 0.4)
    ax.set_ylim(-0.4, L + 0.4)
    ax.set_aspect('equal')
    ax.set_title(_t('duality_title', lang) +
                 r'  ($v/c$ = %.1f,  $\alpha$ = arctan($v/c$) = %.1f°)' % (
                     beta, math.degrees(alpha)),
                 fontsize=12)
    ax.grid(True, alpha=0.15)
    # Remove ticks for clean look
    ax.set_xticks([])
    ax.set_yticks([])
    try:
        plt.tight_layout()
    except ValueError:
        pass
    plt.show()
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
    Bar chart comparing Newton, Einstein, and ORT for a given scenario.
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
# 6b. Unit-circle diagram (time dilation + length contraction)
# =============================================================================

def unit_circle_diagram(theta_deg=30, lang='nl', figsize=(14, 6)):
    """Unit circle showing time dilation and length contraction as projections.

    The same circle (radius = 1) represents both:
    - 1 second on the time axis -> proper time tau = cos(theta)
    - rest length L_0 perpendicular to velocity -> measured L = L_0 cos(theta)
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    theta = np.radians(theta_deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    def _draw_bg(ax):
        a = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(a), np.sin(a), 'k-', lw=0.8, alpha=0.25)
        ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.3)
        ax.set_aspect('equal')

    # --- Left panel: time dilation ---
    _draw_bg(ax1)
    # Velocity vector at angle theta from time axis (y-axis)
    vx, vy = sin_t, cos_t
    ax1.annotate('', xy=(vx, vy), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='->', color='#2070c0', lw=2.5))
    ax1.text(vx / 2 - 0.14, vy / 2 + 0.08, r'$\vec{c}$', fontsize=16, color='#2070c0')

    # 1 second mark on time axis
    ax1.plot([-0.03, 0.03], [1, 1], 'k-', lw=2)
    ax1.text(-0.12, 1.0, _t('unit_circle_1s', lang), fontsize=11, ha='right', va='center')

    # Projection onto time axis
    ax1.plot([vx, 0], [vy, vy], 'r--', lw=1.5, alpha=0.6)
    ax1.plot([0, 0], [0, vy], 'r-', lw=4, alpha=0.7)
    ax1.plot([-0.03, 0.03], [vy, vy], 'r-', lw=2)
    ax1.text(-0.12, vy / 2, r'$\cos\theta$', fontsize=15, color='red', ha='right', va='center')

    # theta arc
    arc1 = Arc((0, 0), 0.35, 0.35, angle=90, theta1=-theta_deg, theta2=0,
               color='purple', lw=2)
    ax1.add_patch(arc1)
    ax1.text(0.07, 0.22, r'$\theta$', fontsize=14, color='purple')

    ax1.set_xlabel(r'$v_{%s}$ / c' % _t('unit_circle_space', lang), fontsize=13)
    ax1.set_ylabel(r'$v_{%s}$ / c' % _t('unit_circle_time', lang), fontsize=13)
    ax1.set_title(_t('unit_circle_td', lang), fontsize=14)
    ax1.set_xlim(-0.3, 1.2)
    ax1.set_ylim(-0.3, 1.2)

    # --- Right panel: length contraction ---
    _draw_bg(ax2)
    # Velocity vector
    ax2.annotate('', xy=(vx, vy), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='->', color='#2070c0', lw=2.5))
    ax2.text(vx / 2 - 0.14, vy / 2 + 0.08, r'$\vec{c}$', fontsize=16, color='#2070c0')

    # Rod perpendicular to velocity vector (length = 1 = L_0)
    px = cos_t
    py = -sin_t
    ax2.annotate('', xy=(px, py), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='->', color='#e07020', lw=2.5))
    ax2.text(px / 2 + 0.08, py / 2 + 0.08, _t('unit_circle_1ls', lang),
             fontsize=15, color='#e07020')

    # Right-angle mark
    v_hat = np.array([sin_t, cos_t])
    p_hat = np.array([cos_t, -sin_t])
    sq = 0.05
    ax2.plot([sq * v_hat[0], sq * (v_hat[0] + p_hat[0]), sq * p_hat[0]],
             [sq * v_hat[1], sq * (v_hat[1] + p_hat[1]), sq * p_hat[1]], 'k-', lw=1)

    # Projection of rod onto space axis (horizontal)
    ax2.plot([px, px], [py, 0], 'g--', lw=1.5, alpha=0.5)
    ax2.plot([0, px], [0, 0], 'g-', lw=4, alpha=0.7)
    ax2.text(px / 2, -0.12, r'$L_0 \cos\theta$', fontsize=14, color='green', ha='center')

    # theta arc
    arc2 = Arc((0, 0), 0.35, 0.35, angle=90, theta1=-theta_deg, theta2=0,
               color='purple', lw=2)
    ax2.add_patch(arc2)
    ax2.text(0.07, 0.22, r'$\theta$', fontsize=14, color='purple')

    # theta arc at projection
    arc3 = Arc((0, 0), 0.55, 0.55, angle=0,
               theta1=np.degrees(np.arctan2(py, px)), theta2=0,
               color='purple', lw=1.5, linestyle='--')
    ax2.add_patch(arc3)
    ax2.text(0.32, -0.1, r'$\theta$', fontsize=12, color='purple')

    ax2.set_xlabel(r'$x$ (%s)' % _t('unit_circle_space', lang), fontsize=13)
    ax2.set_ylabel(r'$t$ (%s)' % _t('unit_circle_time', lang), fontsize=13)
    ax2.set_title(_t('unit_circle_lc', lang), fontsize=14)
    ax2.set_xlim(-0.3, 1.2)
    ax2.set_ylim(-0.55, 1.2)

    fig.suptitle(_t('unit_circle_title', lang), fontsize=13, y=1.02)
    plt.tight_layout()
    plt.show()
    return fig


# =============================================================================
# 6b1b. Length contraction from light-line symmetry
# =============================================================================

def lc_symmetry_diagram(theta_deg=30, lang='nl', figsize=(7, 7)):
    """Minkowski spacetime diagram showing that length contraction = time
    dilation by the 45° light-line symmetry.

    Both primed axes (ct' and x') tilt by the same angle α toward the 45°
    light line.  This symmetric tilt means: whatever factor applies to time
    (cos θ = time dilation) must equally apply to space (cos θ = length
    contraction).
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    theta = np.radians(theta_deg)
    beta = np.sin(theta)  # v/c
    gamma = 1.0 / np.cos(theta)

    # Minkowski diagram directions (unnormalized):
    #   ct' axis: direction (beta, 1)  — tilts from ct toward light line
    #   x'  axis: direction (1, beta)  — tilts from x  toward light line
    # Both tilt by angle alpha = arctan(beta) toward the 45° diagonal.
    alpha_deg = np.degrees(np.arctan(beta))

    ax_len = 1.4

    # --- Rest-frame axes (gray) ---
    ax.annotate('', xy=(ax_len, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    ax.annotate('', xy=(0, ax_len), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    ax.text(ax_len + 0.03, -0.06, '$x$', fontsize=14, color='gray')
    ax.text(-0.08, ax_len + 0.02, '$ct$', fontsize=14, color='gray')

    # --- 45° light line ---
    ax.plot([0, 1.3], [0, 1.3], color='#e6a817', lw=2, ls='--', alpha=0.7)
    ax.text(1.08, 1.18, _t('lc_sym_light', lang),
            fontsize=11, color='#c49000', rotation=45, ha='center',
            style='italic')

    # --- Moving ct' axis (blue) ---
    ct_end = np.array([beta, 1.0]) * 1.2
    ax.annotate('', xy=ct_end, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='#2070c0', lw=2.5))
    ax.text(ct_end[0] + 0.05, ct_end[1] + 0.02, "$ct'$",
            fontsize=13, color='#2070c0', style='italic')

    # --- Moving x' axis (red) ---
    x_end = np.array([1.0, beta]) * 1.2
    ax.annotate('', xy=x_end, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='#c04020', lw=2.5))
    ax.text(x_end[0] + 0.03, x_end[1] - 0.07, "$x'$",
            fontsize=13, color='#c04020', style='italic')

    # --- Time dilation parallelogram (blue) ---
    # 1 second on ct axis → draw line of simultaneity (parallel to x')
    # to find where it meets ct' axis.
    # Point on ct at height 1: (0, 1)
    # Line of simultaneity through (0,1) parallel to x'-direction (1, beta):
    #   parametric: (t, 1 + t*beta)
    # Intersect with ct' axis direction (beta, 1):
    #   (t, 1 + t*beta) = s*(beta, 1)  →  t = s*beta, 1+t*beta = s
    #   → s = 1 + s*beta^2  → s = 1/(1-beta^2) = gamma^2
    #   → intersection at s*gamma^2*(beta, 1)  ... that's far away.
    # Actually: let me just show the "1 proper second" on ct'.
    # In Minkowski geometry, 1 proper second along ct' is at
    # (x, ct) = (beta*gamma, gamma).  Its ct-coordinate is gamma = 1/cos θ.
    # The REST-FRAME observer reads ct = gamma for that event → time DILATION.
    # Conversely, ct = 1 in rest frame corresponds to tau = cos θ proper time.

    # Show: 1 unit on ct, and the corresponding proper time = cos θ.
    # Mark ct = 1
    ax.plot([-0.02, 0.02], [1, 1], 'k-', lw=1.5)
    ax.text(-0.05, 1.0, '$1$', fontsize=11, color='gray', ha='right', va='center')

    # Line of simultaneity in S' through ct=1 on ct-axis: parallel to x' dir
    # From (0, 1), direction (1, beta): parametric (t, 1+t*beta)
    # This meets ct' at: (t, 1+t*beta) = s*(beta, 1) → t = s*beta, 1+t*beta = s
    # → 1 + s*beta^2 = s → s = 1/(1-beta^2) = gamma^2
    # But we only need the visual line segment. Draw from (0,1) a short way.
    sim_len = 0.4
    sim_end = np.array([sim_len, 1 + sim_len * beta])
    ax.plot([0, sim_end[0]], [1, sim_end[1]], '#2070c0', ls=':', lw=1.5, alpha=0.6)

    # The proper time for this rest-frame second is cos θ.
    # Annotate on the left.
    ax.text(-0.06, 0.5, _t('lc_sym_td', lang),
            fontsize=9, color='#2070c0', ha='right', va='center', rotation=90)

    # --- Length contraction parallelogram (red) ---
    # 1 light-second on x axis → draw world line (parallel to ct')
    # to find where it meets x' axis.
    ax.plot([1, 1], [-0.02, 0.02], 'k-', lw=1.5)
    ax.text(1.0, -0.05, '$1$', fontsize=11, color='gray', ha='center', va='top')

    # World line through (1, 0) parallel to ct' direction (beta, 1):
    wl_len = 0.4
    wl_end = np.array([1 + wl_len * beta, wl_len])
    ax.plot([1, wl_end[0]], [0, wl_end[1]], '#c04020', ls=':', lw=1.5, alpha=0.6)

    ax.text(0.5, -0.06, _t('lc_sym_lc', lang),
            fontsize=9, color='#c04020', ha='center', va='top')

    # --- Angle α labels showing symmetric tilt ---
    # α from ct-axis to ct'
    arc = Arc((0, 0), 0.5, 0.5, angle=90, theta1=-alpha_deg, theta2=0,
              color='purple', lw=1.5)
    ax.add_patch(arc)
    ax.text(0.09, 0.32, r'$\alpha$', fontsize=13, color='purple')

    # α from x-axis to x'
    arc2 = Arc((0, 0), 0.5, 0.5, angle=0, theta1=0, theta2=alpha_deg,
               color='purple', lw=1.5)
    ax.add_patch(arc2)
    ax.text(0.34, 0.08, r'$\alpha$', fontsize=13, color='purple')

    # --- cos θ annotations ---
    # For time: proper time of 1 rest-frame second = cos θ
    cos_t = np.cos(theta)
    # Small bracket or annotation near ct axis
    ax.annotate(r'$\cos\theta$',
                xy=(0, cos_t), xytext=(-0.15, cos_t),
                fontsize=12, color='#2070c0', fontweight='bold',
                ha='right', va='center')

    # For space: proper length of 1 rest-frame light-second = cos θ
    ax.annotate(r'$\cos\theta$',
                xy=(cos_t, 0), xytext=(cos_t, -0.12),
                fontsize=12, color='#c04020', fontweight='bold',
                ha='center', va='top')

    ax.set_title(_t('lc_sym_title', lang) % theta_deg, fontsize=14)
    ax.set_aspect('equal')
    ax.set_xlim(-0.3, 1.5)
    ax.set_ylim(-0.2, 1.5)
    ax.axis('off')

    plt.tight_layout()
    plt.show()
    return fig


# =============================================================================
# 6b2. Observer axes diagram — unit vectors rotated by θ
# =============================================================================

def observer_axes_diagram(theta_deg=30, lang='nl', figsize=(8, 8)):
    """Show the moving observer's unit vectors (time and space directions)
    as two perpendicular arrows rotated by θ inside a unit circle.

    Projections onto the rest-frame axes both give cos(θ), making
    time dilation and length contraction visible at a glance.
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    theta = np.radians(theta_deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Unit circle
    a = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(a), np.sin(a), 'k-', lw=0.8, alpha=0.25)

    # Rest-frame axes
    ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
    ax.axvline(0, color='gray', lw=0.5, alpha=0.3)

    # Tick marks at 1 on the axes
    ax.plot([-0.02, 0.02], [1, 1], 'k-', lw=1.5)
    ax.plot([1, 1], [-0.02, 0.02], 'k-', lw=1.5)

    # --- Time direction of moving observer (blue) ---
    tx, ty = sin_t, cos_t  # on the unit circle, rotated by θ from time axis
    ax.annotate('', xy=(tx, ty), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='#2070c0', lw=2.5))
    ax.text(tx + 0.06, ty + 0.04,
            _t('obs_axes_time_dir', lang), fontsize=12, color='#2070c0',
            style='italic')

    # Projection onto time axis (vertical) — cos θ = time dilation
    ax.plot([tx, 0], [ty, ty], '#2070c0', ls='--', lw=1.5, alpha=0.5)
    ax.plot([0, 0], [0, ty], '#2070c0', lw=4, alpha=0.5)
    ax.plot([-0.03, 0.03], [ty, ty], '#2070c0', lw=2)
    ax.text(-0.08, ty, r'$\cos\theta$', fontsize=14, color='#2070c0',
            ha='right', va='center', fontweight='bold')

    # Projection onto space axis (horizontal) — +sin θ = velocity (v/c)
    ax.plot([tx, tx], [ty, 0], '#2070c0', ls='--', lw=1.5, alpha=0.3)
    ax.plot([0, tx], [0, 0], '#2070c0', lw=4, alpha=0.3)
    ax.plot([tx, tx], [-0.03, 0.03], '#2070c0', lw=2, alpha=0.6)
    ax.text(tx, -0.08, r'$\sin\theta$', fontsize=12, color='#2070c0',
            ha='center', va='top', fontweight='bold', alpha=0.8)

    # --- Space direction of moving observer (red/orange) ---
    sx, sy = cos_t, -sin_t  # perpendicular, on the unit circle
    ax.annotate('', xy=(sx, sy), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='#c04020', lw=2.5))
    ax.text(sx + 0.04, sy - 0.08,
            _t('obs_axes_space_dir', lang), fontsize=12, color='#c04020',
            style='italic')

    # Projection onto space axis (horizontal) — cos θ = length contraction
    ax.plot([sx, sx], [sy, 0], '#c04020', ls='--', lw=1.5, alpha=0.5)
    ax.plot([0, sx], [0, 0], '#c04020', lw=4, alpha=0.5)
    ax.plot([sx, sx], [-0.03, 0.03], '#c04020', lw=2)
    ax.text(sx, -0.08, r'$\cos\theta$', fontsize=14, color='#c04020',
            ha='center', va='top', fontweight='bold')

    # Projection onto time axis (vertical) — -sin θ = simultaneity shift
    ax.plot([sx, 0], [sy, sy], '#c04020', ls='--', lw=1.5, alpha=0.3)
    ax.plot([0, 0], [0, sy], '#c04020', lw=4, alpha=0.3)
    ax.plot([-0.03, 0.03], [sy, sy], '#c04020', lw=2, alpha=0.6)
    ax.text(-0.08, sy, r'$-\sin\theta$', fontsize=12, color='#c04020',
            ha='right', va='center', fontweight='bold', alpha=0.8)

    # Right-angle mark between the two arrows
    sq = 0.06
    t_hat = np.array([tx, ty])
    s_hat = np.array([sx, sy])
    ax.plot([sq * t_hat[0], sq * (t_hat[0] + s_hat[0]), sq * s_hat[0]],
            [sq * t_hat[1], sq * (t_hat[1] + s_hat[1]), sq * s_hat[1]],
            'k-', lw=1)

    # θ arc from vertical axis to time direction
    arc = Arc((0, 0), 0.4, 0.4, angle=90, theta1=-theta_deg, theta2=0,
              color='purple', lw=2)
    ax.add_patch(arc)
    ax.text(0.08, 0.25, r'$\theta$', fontsize=15, color='purple')

    # Annotation labels for the diagonal projections
    ax.text(-0.22, ty / 2, _t('obs_axes_td', lang),
            fontsize=10, color='#2070c0', ha='right', va='center',
            rotation=90)
    ax.text(sx / 2, -0.28, _t('obs_axes_lc', lang),
            fontsize=10, color='#c04020', ha='center', va='top')

    # Annotation labels for the cross-term projections
    ax.text(tx / 2, -0.42, _t('obs_axes_td_cross', lang),
            fontsize=9, color='#2070c0', ha='center', va='top', alpha=0.8)
    ax.text(-0.36, sy / 2, _t('obs_axes_lc_cross', lang),
            fontsize=9, color='#c04020', ha='right', va='center',
            rotation=90, alpha=0.8)

    # Axis labels
    ax.set_xlabel(_t('obs_axes_space', lang), fontsize=14)
    ax.set_ylabel(_t('obs_axes_time', lang), fontsize=14)
    ax.set_title(_t('obs_axes_title', lang) % theta_deg, fontsize=14)

    ax.set_aspect('equal')
    ax.set_xlim(-0.5, 1.3)
    ax.set_ylim(-0.7, 1.3)

    plt.tight_layout()
    plt.show()
    return fig


# =============================================================================
# 6c. Lorentz decomposition diagram (coordinate mixing)
# =============================================================================

def lorentz_decomposition_diagram(theta_deg=30, lang='nl', figsize=(14, 6)):
    """Show how x and ct each mix space and time parts via cos(theta), sin(theta).

    Inverse Lorentz transformation in ORT notation:
        x  = (1/cos theta) x' + c tan theta . t'
        ct = (tan theta) x'   + (1/cos theta) ct'

    Left panel:  x  decomposed into x'-part (blue) and t'-part (red)
    Right panel: ct decomposed into x'-part (blue) and ct'-part (red)
    """
    theta = np.radians(theta_deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    tan_t = np.tan(theta)

    # Pick an event P in the primed frame
    xp = 0.6   # x' value
    ctp = 0.8  # ct' value

    # Inverse Lorentz: compute rest-frame coordinates
    x_total = xp / cos_t + ctp * tan_t
    ct_total = xp * tan_t + ctp / cos_t

    # Decompose
    x_from_xp = xp / cos_t       # space part of x
    x_from_tp = ctp * tan_t      # time part of x
    ct_from_xp = xp * tan_t      # space part of ct
    ct_from_tp = ctp / cos_t     # time part of ct

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    space_color = '#2070c0'  # blue for space (x') contributions
    time_color = '#c03020'   # red for time (t') contributions
    bar_h = 0.25

    # --- Left panel: x coordinate ---
    ax1.barh(0, x_from_xp, height=bar_h, color=space_color, alpha=0.8,
             label=r"$\frac{1}{\cos\theta}\,x'$ = %.3f  (%s)" % (
                 x_from_xp, _t('lorentz_decomp_space_part', lang)))
    ax1.barh(0, x_from_tp, height=bar_h, left=x_from_xp, color=time_color, alpha=0.8,
             label=r"$\tan\theta \cdot ct'$ = %.3f  (%s)" % (
                 x_from_tp, _t('lorentz_decomp_time_part', lang)))

    # Annotation on the bars
    if x_from_xp > 0.15:
        ax1.text(x_from_xp / 2, 0, r'$\frac{x\prime}{\cos\theta}$',
                 ha='center', va='center', fontsize=14, color='white', fontweight='bold')
    if x_from_tp > 0.15:
        ax1.text(x_from_xp + x_from_tp / 2, 0, r'$\tan\theta \cdot ct\prime$',
                 ha='center', va='center', fontsize=14, color='white', fontweight='bold')

    # Total marker
    ax1.axvline(x_total, color='black', lw=2, ls='--', alpha=0.6)
    ax1.text(x_total + 0.02, bar_h * 0.8,
             r'$x = %.3f$' % x_total, fontsize=12, va='bottom')

    ax1.set_xlim(0, x_total * 1.25)
    ax1.set_yticks([])
    ax1.set_xlabel(r'$x$', fontsize=14)
    ax1.set_title(_t('lorentz_decomp_x', lang), fontsize=13)
    ax1.legend(fontsize=11, loc='upper left')
    ax1.grid(True, alpha=0.2, axis='x')

    # Formula below
    ax1.text(0.5, -0.15, r'$x = \frac{1}{\cos\theta}\,x\prime + \tan\theta \cdot ct\prime$',
             transform=ax1.transAxes, ha='center', fontsize=15, style='italic')

    # --- Right panel: ct coordinate ---
    ax2.barh(0, ct_from_xp, height=bar_h, color=space_color, alpha=0.8,
             label=r"$\tan\theta \cdot x'$ = %.3f  (%s)" % (
                 ct_from_xp, _t('lorentz_decomp_space_part', lang)))
    ax2.barh(0, ct_from_tp, height=bar_h, left=ct_from_xp, color=time_color, alpha=0.8,
             label=r"$\frac{1}{\cos\theta}\,ct'$ = %.3f  (%s)" % (
                 ct_from_tp, _t('lorentz_decomp_time_part', lang)))

    if ct_from_xp > 0.15:
        ax2.text(ct_from_xp / 2, 0, r'$\tan\theta \cdot x\prime$',
                 ha='center', va='center', fontsize=14, color='white', fontweight='bold')
    if ct_from_tp > 0.15:
        ax2.text(ct_from_xp + ct_from_tp / 2, 0, r'$\frac{ct\prime}{\cos\theta}$',
                 ha='center', va='center', fontsize=14, color='white', fontweight='bold')

    ax2.axvline(ct_total, color='black', lw=2, ls='--', alpha=0.6)
    ax2.text(ct_total + 0.02, bar_h * 0.8,
             r'$ct = %.3f$' % ct_total, fontsize=12, va='bottom')

    ax2.set_xlim(0, ct_total * 1.25)
    ax2.set_yticks([])
    ax2.set_xlabel(r'$ct$', fontsize=14)
    ax2.set_title(_t('lorentz_decomp_t', lang), fontsize=13)
    ax2.legend(fontsize=11, loc='upper left')
    ax2.grid(True, alpha=0.2, axis='x')

    ax2.text(0.5, -0.15, r'$ct = \tan\theta \cdot x\prime + \frac{1}{\cos\theta}\,ct\prime$',
             transform=ax2.transAxes, ha='center', fontsize=15, style='italic')

    # Suptitle with event info
    fig.suptitle(
        _t('lorentz_decomp_title', lang) +
        r"  ($\theta$ = %g°,  $x'$ = %.1f,  $ct'$ = %.1f)" % (theta_deg, xp, ctp),
        fontsize=13, y=1.02)
    plt.tight_layout()
    plt.show()
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
        msg = ("Interactieve versie — download het notebook om de slider te gebruiken."
               if lang == 'nl' else
               "Interactive version — download the notebook to use the slider.")
        print(msg)
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
# 8a. Isotropic vs Anisotropic Spatial Stretching
# =============================================================================

def spatial_stretching_comparison(lang='nl', figsize=(12, 6)):
    """Two-panel plot comparing Newton, ART and ORT spatial stretching.

    Left panel: proper radial distance vs coordinate r
    Right panel: physical circumference vs coordinate r
    """
    r_over_rs = np.linspace(1.01, 50, 500)

    # --- Proper radial distance (integrated from r_s to r) ---
    # Newton: dl = dr  →  L = r - r_s
    radial_newton = r_over_rs - 1.0

    # ART (Schwarzschild): dl = dr/sqrt(1 - r_s/r)
    # Integral: L = sqrt(r(r-r_s)) + r_s * ln(sqrt(r/r_s - 1) + sqrt(r/r_s))
    x = r_over_rs  # r/r_s
    radial_art = np.sqrt(x * (x - 1)) + np.log(np.sqrt(x - 1) + np.sqrt(x))

    # ORT (isotropic): dl = dr/sqrt(1 - r_s/r)  (same radial integral as ART
    # for the radial direction — the difference shows in the circumference)
    radial_ort = radial_art  # radial stretching is the same

    # --- Physical circumference ---
    # Newton: C = 2πr
    circ_newton = 2 * np.pi * r_over_rs

    # ART (Schwarzschild coords): C = 2πr  (unchanged!)
    circ_art = 2 * np.pi * r_over_rs

    # ORT (isotropic): C = 2πr / sqrt(1 - r_s/r)
    circ_ort = 2 * np.pi * r_over_rs / np.sqrt(1 - 1.0 / r_over_rs)

    # --- Labels ---
    if lang == 'nl':
        title_radial = 'Eigenafstand (radiaal)'
        title_circ = 'Fysieke omtrek'
        xlabel = r'Coördinaat $r / r_s$'
        ylabel_radial = r'Eigenafstand vanaf $r_s$ [$r_s$]'
        ylabel_circ = r'Omtrek [$r_s$]'
    else:
        title_radial = 'Proper radial distance'
        title_circ = 'Physical circumference'
        xlabel = r'Coordinate $r / r_s$'
        ylabel_radial = r'Proper distance from $r_s$ [$r_s$]'
        ylabel_circ = r'Circumference [$r_s$]'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left panel: radial distance (ART and ORT are identical here)
    ax1.plot(r_over_rs, radial_newton, 'b--', linewidth=2, label='Newton')
    ax1.plot(r_over_rs, radial_art, 'r-', linewidth=2, label='ART = ORT')
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel_radial, fontsize=12)
    ax1.set_title(title_radial, fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(1, 50)

    # Right panel: circumference
    ax2.plot(r_over_rs, circ_newton, 'b--', linewidth=2, label='Newton')
    ax2.plot(r_over_rs, circ_art, 'r-', linewidth=2, label='ART')
    ax2.plot(r_over_rs, circ_ort, 'g-', linewidth=2, label='ORT')
    ax2.set_xlabel(xlabel, fontsize=12)
    ax2.set_ylabel(ylabel_circ, fontsize=12)
    ax2.set_title(title_circ, fontsize=13, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
    return fig


# =============================================================================
# 8b. Spacetime Embedding (Flamm's Paraboloid)
# =============================================================================

def spacetime_embedding_3d(mass=None, lang='nl', figsize=(9, 7)):
    """3D Flamm's paraboloid using matplotlib (static, works in jupyter-book)."""
    from mpl_toolkits.mplot3d import Axes3D

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

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z / rs, cmap='viridis', alpha=0.8,
                    edgecolor='none', rcount=100, ccount=100)
    # Symmetric axis limits so the funnel is visually centered
    lim = 10
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel('x / $r_s$', fontsize=11)
    ax.set_ylabel('y / $r_s$', fontsize=11)
    ax.set_zlabel('z / $r_s$', fontsize=11)
    ax.set_title(_t('embedding_title', lang), fontsize=14, fontweight='bold')
    ax.view_init(elev=25, azim=-60)
    plt.tight_layout()
    plt.show()
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
            label=r'GR / ORT ($\alpha$)')

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
        msg = ("Interactieve versie — download het notebook om de slider te gebruiken."
               if lang == 'nl' else
               "Interactive version — download the notebook to use the slider.")
        print(msg)
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
               'ART (1915)', 'ORT', 'Thomas']

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
