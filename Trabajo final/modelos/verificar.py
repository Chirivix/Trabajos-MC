"""Verificaciones de los niveles de modelos_progresivos.py. Ejecutar: python verificar.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modelos_progresivos import *

L = Lata(R=0.054, m_c=0.0122, m_l=0.193, fill=0.5)   # datos del ensayo del cuaderno final
ok = True
_trapz = getattr(np, 'trapezoid', None) or np.trapz
def chequeo(nombre, valor, tol, texto=""):
    global ok
    bien = valor <= tol
    ok &= bool(bien)
    print(f"[{'OK' if bien else 'FALLA'}] {nombre}: {valor:.3e} (tolerancia {tol:.0e}) {texto}")

t = np.linspace(0, 1.5, 3001)

# 1) Nivel 3 -> Nivel 1 cuando la masa de oleaje tiende a cero: theta = omega0 * t
Lc = Lata(R=L.R, m_c=L.m_c, m_l=L.m_l, mu_override=1e-9)
th, om, ph, dph = simular_oleaje(Lc, 6.0, 0.1, 0.0, t)
chequeo("Nivel 3 -> Nivel 1 (m1 -> 0): max|theta - omega0 t|", np.max(np.abs(th - 6.0 * t)), 1e-6)

# 2) Conservación de la energía y de p_theta (sin disipación)
th, om, ph, dph = simular_oleaje(L, 6.0, 0.2, 1.0, t)
E = energia_oleaje(L, om, ph, dph); p = momento_theta(L, om, ph, dph)
chequeo("Energía: deriva relativa", np.max(np.abs(E - E[0])) / abs(E[0]), 1e-8)
chequeo("p_theta: deriva relativa", np.max(np.abs(p - p[0])) / abs(p[0]), 1e-8)

# 3) Hamilton frente a Euler-Lagrange
thH, omH, phH, dphH = simular_hamilton(L, 6.0, 0.2, 1.0, t)
chequeo("Hamilton vs Euler-Lagrange: max|dtheta|", np.max(np.abs(th - thH)), 1e-7)
chequeo("Hamilton vs Euler-Lagrange: max|dphi|", np.max(np.abs(ph - phH)), 1e-7)

# 4) Frecuencia del modo acoplado (pequeñas oscilaciones, lata en reposo): w^2 = w_s^2 / (1 - m1 R^2 / J_bar)
tl = np.linspace(0, 20, 40001)
th, om, ph, dph = simular_oleaje(L, 0.0, 1e-3, 0.0, tl)
idx = np.where(np.sign(ph[:-1]) * np.sign(ph[1:]) < 0)[0]          # cruces por cero de phi
tc = tl[idx] - ph[idx] * (tl[idx + 1] - tl[idx]) / (ph[idx + 1] - ph[idx])
w_num = np.pi / np.mean(np.diff(tc))                                  # medio periodo entre cruces
chequeo("Frecuencia acoplada: |w_num - w_teo|/w_teo", abs(w_num - L.omega_acoplada) / L.omega_acoplada, 1e-3,
        f"(w_teo = {L.omega_acoplada:.3f} rad/s, w_s = {L.omega_s:.3f} rad/s)")

# 5) Disipación: dE/dt = - c_rod theta'^2 - c_visc phi'^2 (función de Rayleigh)
cv, cr = 3e-5, 2e-5
th, om, ph, dph = simular_oleaje(L, 6.0, 0.2, 1.0, t, c_visc=cv, c_rod=cr)
E = energia_oleaje(L, om, ph, dph)
perdida = _trapz(cr * om**2 + cv * dph**2, t)
chequeo("Balance de energía con disipación: |dE + trabajo disipado| / |dE|", abs((E[0] - E[-1]) - perdida) / abs(E[0] - E[-1]), 2e-5)
chequeo("Energía no creciente (máximo aumento)", max(0.0, np.max(np.diff(E))), 1e-9)

# 6) Cicloide ideal (Nivel 1): la marca llega a 2*R*Omega arriba y a 0 en el contacto
Om = 8.0
th1, om1 = simular_nivel1(L, Om, t)
dx, dy = posicion_marca(L.R, th1, np.pi / 2)
vx = np.gradient(dx, t)
chequeo("Cicloide: |v_max - 2 R Omega| / (2 R Omega)", abs(vx.max() - 2 * L.R * Om) / (2 * L.R * Om), 1e-3)
chequeo("Cicloide: v_min (contacto)", abs(vx.min()) / (2 * L.R * Om), 1e-3)

# 7) Nivel 2 frente al script lata_rodando.py (si está disponible)
try:
    import lata_rodando as lr
    Lr = lr.Lata(R=L.R, m_c=L.m_c, m_l=L.m_l, kappa=1.0)
    tt = np.linspace(0, 1.0, 400)
    thA, _ = simular_nivel2(L, 0.4, 14.0, tt)
    thB = lr.simular(Lr, 0.4, 14.0, tt)["theta"]
    chequeo("Nivel 2 vs lata_rodando.py: max|dtheta|", np.max(np.abs(thA - thB)), 1e-6)
except ImportError:
    print("[--] lata_rodando.py no encontrado en la carpeta superior: se omite la prueba 7")

print("\nTODO CORRECTO" if ok else "\nHAY VERIFICACIONES QUE FALLAN")
sys.exit(0 if ok else 1)
