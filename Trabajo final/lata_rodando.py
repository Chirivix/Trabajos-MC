#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lata_rodando.py
===============
Reto IV-2 · Simulación de una lata con medio cilindro de líquido que rueda
(modelo de las Etapas I-III: líquido rígido y solidario a la lata, 1 g.d.l.).

MODELO (resumen)
----------------
  e   = 4 R_int / (3 pi)                    centroide del semidisco de líquido
  M   = m_c + m_l
  a   = m_l e / M                           excentricidad del CM del sistema
  I_O = kappa m_c R^2 + (1/2) m_l R_int^2   inercia respecto al centro geométrico
  J(th) = I_O + M R^2 - 2 M R a cos(th)     inercia respecto al punto de contacto

  Ecuación de movimiento (Lagrange, rodadura sin deslizar, x_O = R th):
      J(th) th'' + M R a sin(th) th'^2 + M g a sin(th) = 0

  th = 0  <->  CM del sistema justo debajo del centro (equilibrio estable).
  th crece en el sentido de la rodadura (x_O = R th).

USO (desde terminal)
--------------------
  1) Verificar el código (sin datos):
        python lata_rodando.py --selftest

  2) Simular con parámetros medidos:
        python lata_rodando.py --R <m> --mc <kg> --ml <kg> --kappa <adim> \
               --theta0 <rad> --omega0 <rad/s> --tmax <s>

  3) Contrastar con Tracker (CSV/TSV exportado, columnas t y x del centro):
        python lata_rodando.py --R <m> --mc <kg> --ml <kg> --kappa <adim> \
               --csv datos.csv --tcol 0 --xcol 1 --theta0 <rad>
     Añadir --fit theta0,omega0   para refinar SOLO las condiciones iniciales,
     o   --fit theta0,omega0,kappa   si además se quiere estimar kappa.

No hay ningún valor experimental dentro de este archivo: todo entra por
parámetros. Las únicas constantes internas son g y las del autotest, que
usa unidades adimensionales (R = 1, g = 1).

Dependencias: numpy, scipy, matplotlib.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.signal import savgol_filter

G_STD = 9.80665  # m/s^2


# =============================================================================
# 1. GEOMETRÍA Y PARÁMETROS DERIVADOS (Etapa I)
# =============================================================================
@dataclass
class Lata:
    """Parámetros de la lata. Unidades SI.

    R      : radio exterior (el que rueda: x_O = R th)
    m_c    : masa de la envoltura (lata vacía)
    m_l    : masa del líquido (medio cilindro, cara plana por el eje)
    kappa  : I_c = kappa m_c R^2   (1 = sólo pared lateral; 0.5 = disco macizo)
    R_int  : radio interior (por defecto = R)
    g      : gravedad
    """
    R: float
    m_c: float
    m_l: float
    kappa: float = 1.0
    R_int: float | None = None
    g: float = G_STD

    def __post_init__(self):
        if self.R_int is None:
            self.R_int = self.R

    # --- magnitudes geométricas -------------------------------------------
    @property
    def e(self) -> float:          # centroide del semidisco respecto a O
        return 4.0 * self.R_int / (3.0 * np.pi)

    @property
    def M(self) -> float:
        return self.m_c + self.m_l

    @property
    def a(self) -> float:          # excentricidad del CM total
        return self.m_l * self.e / self.M

    @property
    def I_O(self) -> float:        # inercia respecto al centro geométrico
        return self.kappa * self.m_c * self.R**2 + 0.5 * self.m_l * self.R_int**2

    @property
    def I_cm(self) -> float:       # inercia respecto al CM (Steiner)
        return self.I_O - self.M * self.a**2

    @property
    def J_bar(self) -> float:      # parte constante de J(th)
        return self.I_O + self.M * self.R**2

    # --- pequeñas oscilaciones (Etapa III) --------------------------------
    @property
    def eps(self) -> float:        # parámetro de modulación
        return self.M * self.R * self.a / self.J_bar

    @property
    def omega0(self) -> float:     # frecuencia de balanceo alrededor de th=0
        return np.sqrt(self.M * self.g * self.a
                       / (self.J_bar - 2.0 * self.M * self.R * self.a))

    def J(self, th):
        return self.J_bar - 2.0 * self.M * self.R * self.a * np.cos(th)

    def resumen(self) -> str:
        return (f"  R={self.R:.5g} m  R_int={self.R_int:.5g} m  m_c={self.m_c:.5g} kg  "
                f"m_l={self.m_l:.5g} kg  kappa={self.kappa:.4g}\n"
                f"  e={self.e:.5g} m  a={self.a:.5g} m  I_O={self.I_O:.5g} kg m^2  "
                f"I_cm={self.I_cm:.5g} kg m^2\n"
                f"  eps={self.eps:.4f}  omega0={self.omega0:.5g} rad/s  "
                f"f0={self.omega0 / (2 * np.pi):.5g} Hz  T0={2 * np.pi / self.omega0:.5g} s")


# =============================================================================
# 2. DINÁMICA (Etapa II)
# =============================================================================
def aceleracion(lata: Lata, th, om):
    """th'' a partir de la ecuación de Lagrange."""
    s = np.sin(th)
    return -(lata.M * lata.R * lata.a * s * om**2 + lata.M * lata.g * lata.a * s) / lata.J(th)


def _rhs(t, y, lata: Lata):
    th, om = y
    return [om, aceleracion(lata, th, om)]


def energia(lata: Lata, th, om):
    """E = 1/2 J th'^2 - M g a cos(th)   (constante en el modelo)."""
    return 0.5 * lata.J(th) * om**2 - lata.M * lata.g * lata.a * np.cos(th)


def regimen(lata: Lata, th0: float, om0: float) -> str:
    """Clasifica el movimiento por energía: ¿llega el CM al punto más alto (th = pi)?"""
    E = energia(lata, th0, om0)
    E_top = lata.M * lata.g * lata.a          # E en th = pi con om = 0
    if E < E_top:
        return "vaivén (el CM no llega arriba): oscilación alrededor de th = 0"
    return "rodadura continua (el CM pasa por arriba): la velocidad oscila a la frecuencia de giro"


def simular(lata: Lata, th0: float, om0: float, t_eval: np.ndarray,
            rtol: float = 1e-10, atol: float = 1e-12) -> dict:
    """Integra la ecuación de movimiento con DOP853 (Runge-Kutta explícito de orden 8).

    Justificación: sistema conservativo, suave, de 1 g.d.l. y sin rigidez (stiffness);
    un RK de alto orden con paso adaptativo y tolerancias estrictas mantiene la
    deriva de energía ~1e-8 relativa. (Si en la Etapa IV aparece viscosidad con
    escalas de tiempo muy distintas, habrá que revisar y usar un método para
    sistemas rígidos, p. ej. Radau o LSODA.)
    """
    t_eval = np.asarray(t_eval, float)
    sol = solve_ivp(_rhs, (t_eval[0], t_eval[-1]), [th0, om0], args=(lata,),
                    method="DOP853", t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    th, om = sol.y
    al = aceleracion(lata, th, om)
    R, a, M, g = lata.R, lata.a, lata.M, lata.g
    return dict(
        t=sol.t, theta=th, omega=om, alpha=al,
        x_O=R * (th - th0),                       # desplazamiento del centro geométrico
        v_O=R * om,
        X_cm=R * th - a * np.sin(th) - R * th0 + a * np.sin(th0),
        Y_cm=R - a * np.cos(th),
        N=M * (g + a * np.cos(th) * om**2 + a * np.sin(th) * al),          # normal
        F=M * ((R - a * np.cos(th)) * al + a * np.sin(th) * om**2),        # fricción
        E=energia(lata, th, om),
    )


def solucion_lineal(lata: Lata, th0: float, om0: float, t: np.ndarray) -> np.ndarray:
    """th(t) linealizado alrededor de th=0 (válido para amplitudes pequeñas)."""
    w = lata.omega0
    return th0 * np.cos(w * t) + (om0 / w) * np.sin(w * t)


# =============================================================================
# 3. AUTOTEST (unidades adimensionales: R = 1, g = 1)
# =============================================================================
def autotest() -> bool:
    ok = True

    def check(nombre, cond, detalle=""):
        nonlocal ok
        print(f"  [{'OK' if cond else 'FALLA'}] {nombre} {detalle}")
        ok &= bool(cond)

    print("Autotest (unidades adimensionales R=1, g=1)")

    # 1. Control analítico: m_c -> 0 => omega0^2 = 0.652 g/R
    l0 = Lata(R=1, m_c=0, m_l=1, kappa=1, g=1)
    check("omega0^2 (m_c=0) = 0.652 g/R", abs(l0.omega0**2 - 0.6517) < 2e-3,
          f"-> {l0.omega0**2:.4f}")

    # 2. Si m_l -> 0 no hay excentricidad: th'' = 0 (rodadura uniforme)
    l1 = Lata(R=1, m_c=1, m_l=1e-12, kappa=1, g=1)
    r = simular(l1, 0.3, 2.0, np.linspace(0, 20, 400))
    check("a -> 0 : velocidad uniforme", np.ptp(r["omega"]) < 1e-9)

    # 3. Conservación de la energía en gran amplitud (vaivén y rodadura continua)
    l2 = Lata(R=1, m_c=0.5, m_l=1, kappa=0.8, g=1)
    for om0, etiqueta in [(1.5, "vaivén grande"), (4.0, "rodadura continua")]:
        r = simular(l2, 0.0, om0, np.linspace(0, 100, 4000))
        deriva = np.ptp(r["E"]) / abs(r["E"][0])
        check(f"energía constante ({etiqueta})", deriva < 1e-6, f"-> deriva {deriva:.1e}")

    # 4. Ecuación lineal vs numérica con amplitud pequeña
    t = np.linspace(0, 30, 3000)
    th0, om0 = 0.0, 0.01 * l2.omega0
    r = simular(l2, th0, om0, t)
    err = np.max(np.abs(r["theta"] - solucion_lineal(l2, th0, om0, t)))
    check("lineal vs numérico (amplitud 0.01 rad)", err < 1e-4,
          f"-> error máx {err:.1e} rad")

    # 5. Periodo numérico vs 2*pi/omega0 (amplitud pequeña)
    th = r["theta"]
    cruces = t[:-1][(th[:-1] < 0) & (th[1:] >= 0)]
    if len(cruces) >= 2:
        T_num = np.mean(np.diff(cruces))
        T_teo = 2 * np.pi / l2.omega0
        check("periodo numérico = 2pi/omega0", abs(T_num / T_teo - 1) < 1e-3,
              f"-> {T_num:.4f} vs {T_teo:.4f}")

    # 6. Coherencia Lagrange <-> Newton-Euler: momento respecto al CM
    #    I_cm th'' + a sin(th) N + (R - a cos(th)) F = 0   (th horario)
    r = simular(l2, 0.7, 1.0, np.linspace(0, 20, 500))
    res = (l2.I_cm * r["alpha"] + l2.a * np.sin(r["theta"]) * r["N"]
           + (l2.R - l2.a * np.cos(r["theta"])) * r["F"])
    check("torque respecto al CM (Newton-Euler)", np.max(np.abs(res)) < 1e-9,
          f"-> residuo {np.max(np.abs(res)):.1e}")

    print("Autotest:", "TODO CORRECTO" if ok else "HAY FALLOS")
    return ok


# =============================================================================
# 4. LECTURA DE DATOS DE TRACKER
# =============================================================================
def leer_tracker(ruta: str) -> np.ndarray:
    """Lee un CSV/TSV exportado desde Tracker y devuelve una matriz numérica.

    Se saltan las líneas que no son numéricas (cabeceras, unidades). Detecta
    separador ';' , tabulador o ',' y admite coma decimal (configuración regional
    en español).
    """
    filas = []
    with open(ruta, encoding="utf-8-sig", errors="replace") as f:
        texto = f.read()
    if ";" in texto:
        sep, coma_decimal = ";", True
    elif "\t" in texto:
        sep, coma_decimal = "\t", True
    else:
        sep, coma_decimal = ",", False
    for linea in texto.splitlines():
        linea = linea.strip()
        if not linea:
            continue
        if coma_decimal:
            linea = linea.replace(",", ".")
        trozos = [p.strip() for p in linea.split(sep)]
        try:
            filas.append([float(p) for p in trozos])
        except ValueError:
            continue          # cabecera o texto
    if not filas:
        raise ValueError("No se encontraron filas numéricas en el archivo.")
    n = min(len(f) for f in filas)
    return np.array([f[:n] for f in filas])


def velocidad_suavizada(t: np.ndarray, u: np.ndarray, ventana: int = 11, orden: int = 3):
    """Derivada de u(t) con Savitzky-Golay (requiere muestreo uniforme, como en Tracker).

    Derivar datos ruidosos amplifica el error ~ sigma/dt; la ventana debe abarcar
    bastante menos que un periodo de oscilación. Se imprime una advertencia si no.
    """
    dt = np.median(np.diff(t))
    ventana = min(ventana, len(t) - (1 - len(t) % 2))
    if ventana % 2 == 0:
        ventana += 1
    return savgol_filter(u, ventana, orden, deriv=1, delta=dt)


def frecuencia_dominante(t: np.ndarray, y: np.ndarray) -> float:
    """Frecuencia (Hz) del pico principal de la FFT de y(t) sin su tendencia lineal."""
    dt = np.median(np.diff(t))
    y = y - np.polyval(np.polyfit(t, y, 1), t)
    n = 8 * len(y)
    espectro = np.abs(np.fft.rfft(y * np.hanning(len(y)), n))
    frec = np.fft.rfftfreq(n, dt)
    return float(frec[1:][np.argmax(espectro[1:])])


# =============================================================================
# 5. COMPARACIÓN MODELO - EXPERIMENTO
# =============================================================================
def metricas(y_exp: np.ndarray, y_mod: np.ndarray) -> dict:
    res = y_exp - y_mod
    ss_res = np.sum(res**2)
    ss_tot = np.sum((y_exp - np.mean(y_exp))**2)
    return dict(RMSE=float(np.sqrt(np.mean(res**2))),
                max_abs=float(np.max(np.abs(res))),
                NRMSE=float(np.sqrt(np.mean(res**2)) / np.ptp(y_exp)),
                R2=float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan)


def desplazamiento_modelo(lata: Lata, th0, om0, t):
    return simular(lata, th0, om0, t)["x_O"]


def ajustar(lata: Lata, t, u, th0, om0, a_ajustar: list[str]):
    """Ajusta por mínimos cuadrados SOLO los parámetros indicados.

    Posibles: 'theta0', 'omega0' (condiciones iniciales) y 'kappa' (único parámetro
    incierto de la envoltura). La masa, la geometría y el líquido NO se ajustan.
    Como la fase th0 hace el problema no convexo, se prueban varios arranques.
    """
    permitidos = {"theta0", "omega0", "kappa"}
    if not set(a_ajustar) <= permitidos:
        raise ValueError(f"Sólo se puede ajustar {permitidos}")

    def desempaquetar(p):
        d = dict(theta0=th0, omega0=om0, kappa=lata.kappa)
        d.update(dict(zip(a_ajustar, p)))
        return d

    def residuo(p):
        d = desempaquetar(p)
        l = Lata(lata.R, lata.m_c, lata.m_l, d["kappa"], lata.R_int, lata.g)
        try:
            return desplazamiento_modelo(l, d["theta0"], d["omega0"], t) - u
        except RuntimeError:
            return np.full_like(u, 1e3)

    lim_inf = {"theta0": -np.pi, "omega0": -np.inf, "kappa": 0.5}
    lim_sup = {"theta0": np.pi, "omega0": np.inf, "kappa": 1.0}
    inicio = {"theta0": th0, "omega0": om0, "kappa": lata.kappa}

    arranques_th = np.linspace(-np.pi, np.pi, 9)[:-1] if "theta0" in a_ajustar else [th0]
    mejor = None
    for th_ini in arranques_th:
        p0 = [th_ini if k == "theta0" else inicio[k] for k in a_ajustar]
        sol = least_squares(residuo, p0,
                            bounds=([lim_inf[k] for k in a_ajustar],
                                    [lim_sup[k] for k in a_ajustar]),
                            x_scale="jac")
        if mejor is None or sol.cost < mejor.cost:
            mejor = sol
    return desempaquetar(mejor.x), mejor


def comparar(lata: Lata, ruta_csv: str, tcol: int, xcol: int, thcol: int | None,
             th_en_grados: bool, th0: float, om0: float | None,
             a_ajustar: list[str], ventana: int, n_ic: int,
             salida_png: str, mostrar: bool):
    import matplotlib
    if not mostrar:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    datos = leer_tracker(ruta_csv)
    t = datos[:, tcol] - datos[0, tcol]
    x = datos[:, xcol]
    sentido = 1.0 if x[-1] >= x[0] else -1.0           # sentido de la rodadura
    u = sentido * (x - x[0])                            # desplazamiento del centro (>0 al rodar)

    # --- condiciones iniciales a partir de los datos -----------------------
    if om0 is None:
        n = max(3, n_ic)
        v_ini = np.polyfit(t[:n], u[:n], 1)[0]          # pendiente de los primeros n puntos
        om0 = v_ini / lata.R
        print(f"omega0 estimada de los primeros {n} puntos: {om0:.4g} rad/s "
              f"(v0 = {v_ini:.4g} m/s)")

    print("\nPARÁMETROS DEL MODELO\n" + lata.resumen())
    print("Régimen predicho:", regimen(lata, th0, om0))
    dt = np.median(np.diff(t))
    T_pred = 2 * np.pi / lata.omega0
    print(f"Muestreo: dt = {dt:.4g} s ({1/dt:.4g} fps equivalentes); "
          f"T0 (pequeñas oscilaciones) = {T_pred:.4g} s -> {T_pred/dt:.1f} muestras/periodo")

    # --- ajuste opcional ----------------------------------------------------
    lata_uso = lata
    if a_ajustar:
        params, sol = ajustar(lata, t, u, th0, om0, a_ajustar)
        print("\nAJUSTE (solo " + ", ".join(a_ajustar) + "):", {k: round(float(v), 5) for k, v in params.items()})
        th0, om0 = params["theta0"], params["omega0"]
        lata_uso = Lata(lata.R, lata.m_c, lata.m_l, params["kappa"], lata.R_int, lata.g)
        print("  AVISO: los valores ajustados son estimaciones; no son una predicción del modelo.")

    # --- simulación y métricas ---------------------------------------------
    r = simular(lata_uso, th0, om0, t)
    v_exp = velocidad_suavizada(t, u, ventana)
    m_x = metricas(u, r["x_O"])
    m_v = metricas(v_exp, r["v_O"])
    print("\nMÉTRICAS  (x_O: desplazamiento del centro, v: velocidad suavizada)")
    print(f"  x : RMSE = {m_x['RMSE']:.4g} m   máx|res| = {m_x['max_abs']:.4g} m   "
          f"NRMSE = {m_x['NRMSE']:.2%}   R2 = {m_x['R2']:.4f}")
    print(f"  v : RMSE = {m_v['RMSE']:.4g} m/s máx|res| = {m_v['max_abs']:.4g} m/s "
          f"NRMSE = {m_v['NRMSE']:.2%}   R2 = {m_v['R2']:.4f}")
    f_exp = frecuencia_dominante(t, v_exp)
    f_mod = frecuencia_dominante(t, r["v_O"])
    print(f"  Frecuencia dominante de v(t): experimento {f_exp:.4g} Hz | modelo {f_mod:.4g} Hz")
    print(f"  Condición de no deslizamiento: mu_s mínimo = max|F/N| = "
          f"{np.max(np.abs(r['F'] / r['N'])):.3g};  N mínimo = {np.min(r['N']):.4g} N "
          f"({'sin despegue' if np.min(r['N']) > 0 else 'HAY DESPEGUE: modelo no válido'})")

    # --- ángulo opcional ----------------------------------------------------
    th_exp = None
    if thcol is not None:
        th_exp = datos[:, thcol]
        th_exp = np.deg2rad(th_exp) if th_en_grados else th_exp
        th_exp = np.unwrap(th_exp)
        th_exp = sentido * (th_exp - th_exp[0])         # giro en el sentido de la rodadura
        m_th = metricas(th_exp, r["theta"] - th0)
        print(f"  theta : RMSE = {m_th['RMSE']:.4g} rad   R2 = {m_th['R2']:.4f}")

    # --- gráficas -----------------------------------------------------------
    fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    ax[0].plot(t, u, "k.", ms=4, label="Tracker")
    ax[0].plot(t, r["x_O"], "r-", label="modelo")
    ax[0].set_ylabel("desplazamiento del centro [m]")
    ax[0].legend()
    ax[0].set_title("Lata con medio cilindro de líquido: modelo vs. experimento")
    ax[1].plot(t, v_exp, "k.", ms=4, label="Tracker (Savitzky-Golay)")
    ax[1].plot(t, r["v_O"], "r-", label="modelo")
    ax[1].set_ylabel("velocidad del centro [m/s]")
    ax[1].legend()
    ax[2].plot(t, u - r["x_O"], "b-", label="residuo en x")
    ax[2].axhline(0, color="gray", lw=0.8)
    ax[2].set_ylabel("residuo x [m]")
    ax[2].set_xlabel("t [s]")
    ax[2].legend()
    for a_ in ax:
        a_.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(salida_png, dpi=150)
    print(f"\nFigura guardada en: {salida_png}")
    if mostrar:
        plt.show()
    return r


def solo_simular(lata: Lata, th0: float, om0: float, tmax: float, salida_png: str, mostrar: bool):
    import matplotlib
    if not mostrar:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("PARÁMETROS DEL MODELO\n" + lata.resumen())
    print("Régimen:", regimen(lata, th0, om0))
    t = np.linspace(0, tmax, 4000)
    r = simular(lata, th0, om0, t)
    print(f"Deriva relativa de energía: {np.ptp(r['E']) / abs(r['E'][0]):.1e}")
    print(f"mu_s mínimo = {np.max(np.abs(r['F'] / r['N'])):.3g};  N mínimo = {np.min(r['N']):.4g} N")

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    ax[0, 0].plot(t, r["x_O"]);                       ax[0, 0].set_ylabel("x_O [m]")
    ax[0, 1].plot(t, r["v_O"])
    ax[0, 1].set_ylabel("v_O [m/s]")
    ax[1, 0].plot(t, r["theta"]);                     ax[1, 0].set_ylabel("theta [rad]")
    ax[1, 1].plot(t, r["N"], label="N")
    ax[1, 1].plot(t, r["F"], label="F (fricción)")
    ax[1, 1].set_ylabel("fuerzas [N]")
    ax[1, 1].legend()
    for a_ in ax.ravel():
        a_.set_xlabel("t [s]")
        a_.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(salida_png, dpi=150)
    print(f"Figura guardada en: {salida_png}")
    if mostrar:
        plt.show()


# =============================================================================
# 6. LÍNEA DE COMANDOS
# =============================================================================
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--selftest", action="store_true", help="ejecuta las verificaciones y termina")
    p.add_argument("--R", type=float, help="radio exterior [m]")
    p.add_argument("--Rint", type=float, default=None, help="radio interior [m] (def. = R)")
    p.add_argument("--mc", type=float, help="masa de la lata vacía [kg]")
    p.add_argument("--ml", type=float, help="masa del líquido [kg]")
    p.add_argument("--kappa", type=float, default=1.0, help="I_c = kappa m_c R^2 (def. 1)")
    p.add_argument("--g", type=float, default=G_STD)
    p.add_argument("--theta0", type=float, default=0.0, help="ángulo inicial del CM respecto a la vertical [rad]")
    p.add_argument("--omega0", type=float, default=None, help="velocidad angular inicial [rad/s]")
    p.add_argument("--tmax", type=float, default=5.0, help="duración de la simulación (sin CSV) [s]")
    p.add_argument("--csv", type=str, help="archivo exportado de Tracker")
    p.add_argument("--tcol", type=int, default=0, help="columna del tiempo (0 = primera)")
    p.add_argument("--xcol", type=int, default=1, help="columna de x del centro de la lata")
    p.add_argument("--thcol", type=int, default=None, help="columna del ángulo (opcional)")
    p.add_argument("--deg", action="store_true", help="el ángulo del CSV está en grados")
    p.add_argument("--fit", type=str, default="", help="parámetros a ajustar: theta0,omega0,kappa")
    p.add_argument("--ventana", type=int, default=11, help="ventana Savitzky-Golay (impar)")
    p.add_argument("--nic", type=int, default=5, help="puntos iniciales para estimar v0")
    p.add_argument("--out", type=str, default="resultado.png")
    p.add_argument("--show", action="store_true", help="mostrar la figura en pantalla")
    args = p.parse_args()

    if args.selftest:
        sys.exit(0 if autotest() else 1)

    if None in (args.R, args.mc, args.ml):
        p.error("Se necesitan --R, --mc y --ml (o use --selftest).")

    lata = Lata(args.R, args.mc, args.ml, args.kappa, args.Rint, args.g)
    ajustar_lista = [s.strip() for s in args.fit.split(",") if s.strip()]

    if args.csv:
        comparar(lata, args.csv, args.tcol, args.xcol, args.thcol, args.deg,
                 args.theta0, args.omega0, ajustar_lista, args.ventana, args.nic,
                 args.out, args.show)
    else:
        if args.omega0 is None:
            p.error("Sin CSV hay que dar --omega0.")
        solo_simular(lata, args.theta0, args.omega0, args.tmax, args.out, args.show)


if __name__ == "__main__":
    main()
