# Del modelo al código

`01_problema_parcial.md` termina en ecuaciones. Este documento muestra cómo cada ecuación se convierte en código, nivel por nivel, hasta el cuaderno final. El código de cada nivel está en `modelos/modelos_progresivos.py` y sus comprobaciones en `modelos/verificar.py`.

## 1. La receta: de Euler-Lagrange a `solve_ivp`

Un integrador numérico necesita ecuaciones de **primer orden**: $\dot s = f(t, s)$. El camino es siempre el mismo:

1. Las ecuaciones de Euler-Lagrange se ordenan como $M(q)\,\ddot q = \text{rhs}(q,\dot q)$ (matriz de masa por aceleraciones = todo lo demás).
2. Se define el **vector de estado** $s = (q, \dot q)$. Con dos grados de libertad, $s = (\theta, \dot\theta, \varphi, \dot\varphi)$.
3. Se resuelve el sistema lineal para $\ddot q$. Con una matriz de $2\times2$ se hace con la regla de Cramer, sin invertir nada.
4. La función `f(t, s)` devuelve $(\dot\theta,\ \ddot\theta,\ \dot\varphi,\ \ddot\varphi)$ y se integra con `solve_ivp`.

Se usa el método `DOP853` (Runge-Kutta de orden 8, paso adaptativo). `modelos_progresivos.py` usa tolerancias estrictas (`rtol = 1e-10`) para que las pruebas de conservación de energía den una deriva de ~10⁻¹¹; el cuaderno final usa `rtol = 1e-8`, más rápido para el ajuste. Si más adelante aparecen escalas de tiempo muy distintas (viscosidad muy alta), habrá que revisar el método (por ejemplo `Radau` o `LSODA`).

## 2. Nivel 1: una línea de física

$\bar J\,\ddot\theta = 0$ con la ligadura ya incorporada. La disipación opcional entra como $\ddot\theta = -c_{rod}\dot\theta/\bar J$.

```python
def simular_nivel1(L, omega0, t, c_rod=0.0):
    f = lambda _t, s: [s[1], -c_rod * s[1] / L.J_bar]      # s = (theta, omega)
    sol = solve_ivp(f, (t[0], t[-1]), [0.0, omega0], t_eval=t, method="DOP853", rtol=1e-10, atol=1e-12)
    return sol.y[0], sol.y[1]
```

El resultado se convierte en lo que mide la cámara con `posicion_marca(R, dtheta, alpha0)`:

```python
def posicion_marca(R, dtheta, alpha0):
    dx = R * (dtheta + np.sin(alpha0 + dtheta) - np.sin(alpha0))
    dy = R * (np.cos(alpha0 + dtheta) - np.cos(alpha0))
    return dx, dy
```

## 3. Nivel 2: ecuación no lineal de un grado de libertad

$J(\theta)\,\ddot\theta = -M R a\sin\theta\,\dot\theta^2 - M g a\sin\theta$. Se despeja $\ddot\theta$ dividiendo por $J(\theta)$:

```python
J = lambda th: Jbar - 2.0 * M * R * a * np.cos(th)
f = lambda _t, s: [s[1], -(M*R*a*np.sin(s[0])*s[1]**2 + M*g*a*np.sin(s[0])) / J(s[0])]
```

Aquí $\theta = 0$ significa «CM abajo». En los niveles 1, 3 y 4 el código toma $\theta(0) = 0$ al inicio del video.

## 4. Nivel 3: dos grados de libertad, matriz de masa

Las dos ecuaciones de Euler-Lagrange del Nivel 3 son

$$\bar J\ddot\theta + M_{12}\ddot\varphi = \text{rhs}_1,\qquad M_{12}\ddot\theta + J_s\ddot\varphi = \text{rhs}_2 .$$

Por Cramer, con $\det = \bar J J_s - M_{12}^2$:

$$\ddot\theta = \frac{J_s\,\text{rhs}_1 - M_{12}\,\text{rhs}_2}{\det},\qquad \ddot\varphi = \frac{\bar J\,\text{rhs}_2 - M_{12}\,\text{rhs}_1}{\det}.$$

```python
def f(_t, s):
    th, om, ph, dph = s
    cp, sp = np.cos(ph), np.sin(ph)
    M12  = m1 * R * ls * cp
    rhs1 = m1 * R * ls * sp * dph**2 - c_rod * om        # Nivel 3: c_rod = 0
    rhs2 = -m1 * G * ls * sp - c_visc * dph               # Nivel 3: c_visc = 0
    det  = Jb * Js - M12**2
    return [om, (Js*rhs1 - M12*rhs2)/det, dph, (Jb*rhs2 - M12*rhs1)/det]
```

Los parámetros salen de la clase `Lata`: $m_1 = \mu_s m_l$, $l_s = R/\kappa_s$, $J_s = m_1 l_s^2$, $\bar J = I_c + M R^2$, con $\kappa_s$ y $\mu_s$ interpolados de la tabla según la fracción de llenado.

## 5. Nivel 4: dos términos más

La disipación de Rayleigh agrega **un término a cada lado derecho**: `- c_rod * om` en `rhs1` y `- c_visc * dph` en `rhs2`. Es la misma función `simular_oleaje`; el Nivel 3 es el caso `c_rod = c_visc = 0`. Para la comprobación, la energía se calcula con

```python
def energia_oleaje(L, om, ph, dph):
    M12 = L.m1 * L.R * L.l_s * np.cos(ph)
    return 0.5*L.J_bar*om**2 + M12*om*dph + 0.5*L.J_s*dph**2 - L.m1*G*L.l_s*np.cos(ph)
```

## 6. Hamilton en código

En variables canónicas el estado es $(\theta, \varphi, p_\theta, p_\varphi)$. Se integra $\dot q = \partial H/\partial p$ y $\dot p = -\partial H/\partial q$, y las velocidades se recuperan al final con las fórmulas de $\dot\theta$ y $\dot\varphi$ del documento 01. La única derivada que hace falta es $\partial H/\partial\varphi$, que sale de que $M_{12}$ depende de $\varphi$:

```python
dH_dM12  = -pth*pph/det + N*M12/det**2          # N = Js*pth**2 - 2*M12*pth*pph + Jb*pph**2
dM12_dph = -k*np.sin(ph)                        # k = m1*R*ls
dpph     = -(dH_dM12*dM12_dph + m1*G*ls*np.sin(ph))
```

`simular_hamilton` reproduce la trayectoria de `simular_oleaje` con una diferencia de ~10⁻¹⁰, lo que valida a la vez las dos derivaciones.

## 7. Del modelo a la predicción: ajuste con Tracker

Lo que mide la cámara no es $\theta$, sino la posición del punto marcado. El ajuste sigue estos pasos:

1. **Leer** el archivo de Tracker (tiempo, $x$ y opcionalmente $y$), detectar separador y unidad, y pasar a metros.
2. **Corregir el tiempo.** El archivo está en tiempo de reproducción; se multiplica por (FPS de reproducción)/(FPS reales de grabación) y se recorta lo que ya no es cámara lenta (`tiempo_fin_camara_lenta_s`).
3. **Restar la posición inicial** de $x$ e $y$: el modelo predice desplazamientos.
4. **Estimar $\alpha_0$** (posición inicial de la marca) ajustando la forma geométrica de la trayectoria $(x, y)$. Si no hay $y$, se usa $\pi/2$.
5. **Separar** el primer 45 % de los datos para calibrar y el 55 % restante para evaluar.
6. **Ajustar** por mínimos cuadrados (`least_squares`) los cinco parámetros $(\dot\theta_0, \varphi_0, \dot\varphi_0, \log_{10}c_{visc}, \log_{10}c_{rod})$ con cotas, partiendo de varios puntos iniciales. Los dos coeficientes de amortiguamiento se ajustan en escala logarítmica porque pueden variar en varios órdenes de magnitud.
7. **Predecir** con los parámetros ajustados todo el video y reportar RMSE, error máximo y $R^2$ solo en el 55 % no usado.

El vector de residuos concatena los errores en $x$ y en $y$ de la marca, `np.concatenate([xm - x_cal, ym - y_cal])`.

**Medido, calculado o ajustado.**

| Magnitud | Tipo |
|---|---|
| $R$, $m_c$, $m_l$, $f$, FPS reales | Medidas |
| $\nu$ del líquido | Tabulada (o del usuario) |
| $\kappa_s$, $\mu_s$, $m_1$, $l_s$, $J_s$, $\omega_s$ | Calculadas (tabla en $f$; origen por verificar) |
| $\alpha_0$ | Estimada de la geometría de la trayectoria |
| $\dot\theta_0$, $\varphi_0$, $\dot\varphi_0$, $c_{visc}$, $c_{rod}$ | **Ajustadas** |

## 8. Mapa del cuaderno final, celda por celda

| Celda | Qué hace | Equivale a |
|---|---|---|
| 1. Datos | Masas, $R$, llenado (5 a 95 %), FPS, fin de cámara lenta, unidad, líquido y $\nu$. Convierte a SI y valida. | Parámetros de `Lata` |
| 2. Modelo | Calcula $\kappa_s$, $\mu_s$, $m_1$, $l_s$, $J_s$, $\zeta$ y $c_{visc}$ inicial. | Propiedades de `Lata` |
| 2. `simular` | Resuelve las ecuaciones de Euler-Lagrange del Nivel 4. | `simular_oleaje` |
| 2. `posicion_marca` | Posición de la marca del borde. | `posicion_marca` |
| 2. `leer_tracker` | Lee CSV/TSV, detecta columnas, separador y unidad. | Sección 7, paso 1 |
| 2. Ajuste | Corrección de tiempo, $\alpha_0$, calibración 45 %, `least_squares`. | Sección 7, pasos 2 a 6 |
| 2. Métricas y gráficas | RMSE, error máximo, $R^2$, tres paneles y cinco instantáneas del líquido. | Sección 7, paso 7 |

**Equivalencia con `modelos_progresivos.py`.** Se generaron datos con `simular_oleaje` y el cuaderno, con los parámetros verdaderos, los reproduce con un RMSE de 0,19 mm (el nivel del ruido añadido). Es decir, las ecuaciones y la posición de la marca coinciden.

## 9. Qué comprueba cada prueba de `verificar.py`

| Prueba | Valida |
|---|---|
| Nivel 3 con $m_1\to0$ da $\theta = \omega_0 t$ | Que el Nivel 3 contiene al Nivel 1 |
| Deriva de $E$ y de $p_\theta$ ~10⁻¹¹ | La derivación de $T$, $V$ y de la simetría cíclica de $\theta$ |
| Hamilton frente a Euler-Lagrange | Las dos formulaciones, y la expresión de $\partial H/\partial\varphi$ |
| Frecuencia numérica frente a $\omega_{ac}$ | La linealización y la matriz de masa |
| Balance de energía con disipación | La función de Rayleigh y los signos de $c_{rod}$ y $c_{visc}$ |
| Velocidad de la marca: $2R\Omega$ y $0$ | La fórmula de la cicloide y la convención de $\alpha_0$ |
| Nivel 2 frente a `lata_rodando.py` | Que ambas implementaciones del líquido rígido coinciden |

Ejecución: `python modelos/verificar.py` (desde la carpeta «Trabajo final»; la prueba 7 usa `lata_rodando.py` de esa carpeta).

## 10. Cómo agregar un grado de libertad

1. Escribir la ligadura nueva y las coordenadas de cada masa.
2. Calcular $T$ y $V$; escribir $L$.
3. Derivar Euler-Lagrange a mano y ordenar como $M(q)\ddot q = \text{rhs}$. Si hay disipación, sumar $-\partial\mathcal F/\partial\dot q$ a `rhs`.
4. Pasar a código con el vector de estado ampliado.
5. **Antes de ajustar datos**, comprobar los límites: que el modelo nuevo contenga al anterior cuando el parámetro nuevo tiende a cero, que la energía se conserve sin disipación, y que el balance de energía cierre con disipación. Si hay Hamilton, compararlo con Euler-Lagrange.
