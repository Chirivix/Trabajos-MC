"""
Modelos progresivos de la lata parcialmente llena que rueda (Reto IV-2).

Cada nivel agrega UNA idea física al anterior y se implementa con la misma receta:
    ligaduras -> energías T y V -> lagrangiano L -> ecuaciones de Euler-Lagrange
    -> forma M(q)·q'' = rhs(q, q') -> solve_ivp.

    Nivel 1  Lata homogénea (sin efecto del líquido)      1 GDL: theta
    Nivel 2  Líquido rígido pegado (medio cilindro)       1 GDL: theta   (CM excéntrico)
    Nivel 3  Líquido libre: modo de oleaje pendular       2 GDL: theta, phi
    Nivel 4  Nivel 3 + disipación (rodadura y viscosa)    2 GDL, función de Rayleigh
    Hamilton Nivel 3 en variables (q, p), para contrastar con Euler-Lagrange

El Nivel 4 es el modelo del cuaderno final (notebooks/Mecanica_clasica2.ipynb).
Unidades SI. theta se mide desde el inicio del video (theta(0) = 0) en los niveles 1, 3 y 4.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

G = 9.80665

# Parámetros del modo de oleaje en función de la fracción de llenado, tomados tal cual del
# cuaderno final. ORIGEN NO VERIFICADO: citar la fuente (oleaje lineal en cilindro horizontal)
# antes de usarlos en un informe.
_F_TAB = np.array([0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95])
_KAPPA_TAB = np.array([1.02, 1.05, 1.11, 1.18, 1.26, 1.3557, 1.48, 1.65, 1.92, 2.45, 3.15])
_MU_TAB = np.array([0.92, 0.85, 0.73, 0.62, 0.53, 0.458, 0.38, 0.30, 0.21, 0.11, 0.05])


@dataclass
class Lata:
    R: float            # radio exterior (m)
    m_c: float          # masa de la lata vacía (kg)
    m_l: float          # masa del líquido (kg)
    fill: float = 0.5   # fracción de llenado (0.05 a 0.95)
    nu: float = 1.004e-6  # viscosidad cinemática del líquido (m^2/s); agua a ~20 °C
    mu_override: float | None = None  # solo para pruebas de límite

    # ---- magnitudes comunes a todos los niveles ----
    @property
    def M(self): return self.m_c + self.m_l
    @property
    def I_c(self): return self.m_c * self.R**2            # envoltura (kappa = 1)
    @property
    def J_bar(self): return self.I_c + self.M * self.R**2  # inercia respecto al contacto, sin excentricidad

    # ---- Nivel 2: líquido rígido (medio cilindro, solo 50 %) ----
    @property
    def e(self): return 4.0 * self.R / (3.0 * np.pi)       # centroide del semidisco
    @property
    def a(self): return self.m_l * self.e / self.M          # excentricidad del CM total
    @property
    def I_O(self): return self.I_c + 0.5 * self.m_l * self.R**2

    # ---- Nivel 3: modo de oleaje ----
    @property
    def kappa_s(self): return float(np.interp(self.fill, _F_TAB, _KAPPA_TAB))
    @property
    def mu_s(self):
        return self.mu_override if self.mu_override is not None else float(np.interp(self.fill, _F_TAB, _MU_TAB))
    @property
    def omega_s(self): return np.sqrt(self.kappa_s * G / self.R)   # oleaje con la lata fija
    @property
    def l_s(self): return self.R / self.kappa_s                    # longitud del péndulo equivalente
    @property
    def m1(self): return self.mu_s * self.m_l                      # masa efectiva del oleaje
    @property
    def J_s(self): return self.m1 * self.l_s**2
    @property
    def omega_acoplada(self):
        """Frecuencia del modo de oleaje acoplado a la lata, linealizado en reposo."""
        return self.omega_s / np.sqrt(1.0 - self.m1 * self.R**2 / self.J_bar)
    @property
    def zeta_teo(self): return min(0.85, 2.0 * np.sqrt(self.nu / (self.omega_s * self.R**2)))
    @property
    def c_visc_teo(self): return max(1e-7, 2.0 * self.zeta_teo * self.m1 * self.l_s**2 * self.omega_s)


# ----------------------------------------------------------------------------------
# Nivel 1: lata homogénea.  L = 1/2 J_bar theta'^2   ->   J_bar theta'' = 0
# ----------------------------------------------------------------------------------
def simular_nivel1(L, omega0, t, c_rod=0.0):
    f = lambda _t, s: [s[1], -c_rod * s[1] / L.J_bar]
    sol = solve_ivp(f, (t[0], t[-1]), [0.0, omega0], t_eval=t, method="DOP853", rtol=1e-10, atol=1e-12)
    return sol.y[0], sol.y[1]


# ----------------------------------------------------------------------------------
# Nivel 2: líquido rígido.  L = 1/2 J(theta) theta'^2 + M g a cos(theta)
#   J(theta) theta'' + M R a sin(theta) theta'^2 + M g a sin(theta) = 0
# theta = 0 cuando el CM del sistema está justo debajo del eje.
# ----------------------------------------------------------------------------------
def simular_nivel2(L, theta0, omega0, t):
    M, R, a, g = L.M, L.R, L.a, G
    Jbar = L.I_O + M * R**2
    J = lambda th: Jbar - 2.0 * M * R * a * np.cos(th)
    f = lambda _t, s: [s[1], -(M * R * a * np.sin(s[0]) * s[1]**2 + M * g * a * np.sin(s[0])) / J(s[0])]
    sol = solve_ivp(f, (t[0], t[-1]), [theta0, omega0], t_eval=t, method="DOP853", rtol=1e-10, atol=1e-12)
    return sol.y[0], sol.y[1]


# ----------------------------------------------------------------------------------
# Niveles 3 y 4: lata + péndulo de oleaje.  q = (theta, phi)
#   T = 1/2 J_bar theta'^2 + m1 R l_s cos(phi) theta' phi' + 1/2 J_s phi'^2
#   V = - m1 g l_s cos(phi)
#   Rayleigh: F = 1/2 c_rod theta'^2 + 1/2 c_visc phi'^2
#   M(q) q'' = rhs :   [J_bar, M12; M12, J_s] [theta''; phi''] = [m1 R l_s sin(phi) phi'^2 - c_rod theta';
#                                                              - m1 g l_s sin(phi) - c_visc phi']
#   con M12 = m1 R l_s cos(phi).   Nivel 3: c_rod = c_visc = 0.   Nivel 4: c's > 0.
# ----------------------------------------------------------------------------------
def simular_oleaje(L, omega0, phi0, dphi0, t, c_visc=0.0, c_rod=0.0, rtol=1e-10, atol=1e-12):
    R, m1, ls, Js, Jb = L.R, L.m1, L.l_s, L.J_s, L.J_bar

    def f(_t, s):
        th, om, ph, dph = s
        if m1 <= 0:
            return [om, -(c_rod * om) / Jb, 0.0, 0.0]
        cp, sp = np.cos(ph), np.sin(ph)
        M12 = m1 * R * ls * cp
        rhs1 = m1 * R * ls * sp * dph**2 - c_rod * om
        rhs2 = -m1 * G * ls * sp - c_visc * dph
        det = Jb * Js - M12**2
        return [om, (Js * rhs1 - M12 * rhs2) / det, dph, (Jb * rhs2 - M12 * rhs1) / det]

    sol = solve_ivp(f, (t[0], t[-1]), [0.0, omega0, phi0, dphi0], t_eval=t, method="DOP853", rtol=rtol, atol=atol)
    return sol.y[0], sol.y[1], sol.y[2], sol.y[3]


def energia_oleaje(L, om, ph, dph):
    """E = T + V del sistema lata + oleaje."""
    M12 = L.m1 * L.R * L.l_s * np.cos(ph)
    return 0.5 * L.J_bar * om**2 + M12 * om * dph + 0.5 * L.J_s * dph**2 - L.m1 * G * L.l_s * np.cos(ph)


def momento_theta(L, om, ph, dph):
    """p_theta = dL/dtheta' = J_bar theta' + M12 phi'. Se conserva si c_rod = 0 (theta es cíclica)."""
    return L.J_bar * om + L.m1 * L.R * L.l_s * np.cos(ph) * dph


# ----------------------------------------------------------------------------------
# Hamilton (sin disipación): H = 1/2 p^T M^{-1} p + V
#   theta' =  (J_s p_th - M12 p_ph)/det        phi' = (J_bar p_ph - M12 p_th)/det
#   p_th'  =  0                                 p_ph' = - dH/dphi
# ----------------------------------------------------------------------------------
def simular_hamilton(L, omega0, phi0, dphi0, t, rtol=1e-10, atol=1e-12):
    R, m1, ls, Js, Jb = L.R, L.m1, L.l_s, L.J_s, L.J_bar
    k = m1 * R * ls
    M12_0 = k * np.cos(phi0)
    p_th0 = Jb * omega0 + M12_0 * dphi0
    p_ph0 = M12_0 * omega0 + Js * dphi0

    def f(_t, s):
        th, ph, pth, pph = s
        M12 = k * np.cos(ph)
        det = Jb * Js - M12**2
        N = Js * pth**2 - 2.0 * M12 * pth * pph + Jb * pph**2
        dth = (Js * pth - M12 * pph) / det
        dph = (Jb * pph - M12 * pth) / det
        dH_dM12 = -pth * pph / det + N * M12 / det**2
        dM12_dph = -k * np.sin(ph)
        dpph = -(dH_dM12 * dM12_dph + m1 * G * ls * np.sin(ph))
        return [dth, dph, 0.0, dpph]

    sol = solve_ivp(f, (t[0], t[-1]), [0.0, phi0, p_th0, p_ph0], t_eval=t, method="DOP853", rtol=rtol, atol=atol)
    th, ph, pth, pph = sol.y
    M12 = k * np.cos(ph)
    det = Jb * Js - M12**2
    om = (Js * pth - M12 * pph) / det
    dph = (Jb * pph - M12 * pth) / det
    return th, om, ph, dph


# ----------------------------------------------------------------------------------
# Observable: punto marcado en el borde.  alpha0 se mide desde la vertical hacia arriba.
#   alpha0 = pi/2  -> la marca empieza hacia delante, en la horizontal.
# ----------------------------------------------------------------------------------
def posicion_marca(R, dtheta, alpha0):
    dx = R * (dtheta + np.sin(alpha0 + dtheta) - np.sin(alpha0))
    dy = R * (np.cos(alpha0 + dtheta) - np.cos(alpha0))
    return dx, dy
