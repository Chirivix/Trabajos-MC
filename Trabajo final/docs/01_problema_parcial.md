# Problema de parcial: la lata parcialmente llena que rueda

Este documento plantea el reto como un problema de Mecánica Clásica y lo resuelve en cuatro niveles. Cada nivel agrega **una** idea física al anterior y se resuelve con la misma receta. El Nivel 4 es el modelo del cuaderno final (`notebooks/Mecanica_clasica2.ipynb`). El paso a código está en `02_del_modelo_al_codigo.md`, y cada resultado numérico de aquí se comprueba con `modelos/verificar.py`.

## Enunciado

Una lata cilíndrica de radio $R$ y masa vacía $m_c$ contiene un líquido de masa $m_l$ que ocupa una fracción $f$ de su volumen. La lata rueda sin deslizar, en movimiento plano, sobre un plano horizontal. Un punto del borde está marcado y se graba de lado: la cámara mide su posición $(x, y)$ en función del tiempo.

1. Plantee las ligaduras e identifique los grados de libertad.
2. Escriba la energía cinética, la potencial y el lagrangiano.
3. Obtenga las ecuaciones de Euler-Lagrange (y las de Hamilton cuando haya más de un grado de libertad). Identifique coordenadas cíclicas y cantidades conservadas.
4. Estudie los límites: equilibrio, pequeñas oscilaciones y régimen de rodadura continua.
5. Prediga la trayectoria del punto marcado y compárela con el video.
6. Explique qué produce las oscilaciones de velocidad y cómo dependen del líquido.

## Receta general

1. **Ligaduras** → coordenadas generalizadas. Se cuentan los grados de libertad como coordenadas menos ligaduras.
2. **Energías** $T$ y $V$ en esas coordenadas.
3. **Lagrangiano** $L = T - V$.
4. **Euler-Lagrange**: $\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot q_i} - \dfrac{\partial L}{\partial q_i} = Q_i$, donde $Q_i$ son las fuerzas no conservativas.
5. **Simetrías**: coordenada cíclica ⇒ momento conjugado conservado. Si $L$ no depende de $t$, se conserva la energía.
6. **Linealización** alrededor del equilibrio ⇒ frecuencias.
7. **Código**: llevar las ecuaciones a la forma $M(q)\,\ddot q = \text{rhs}(q,\dot q)$ e integrar numéricamente.

---

## Nivel 1. Lata homogénea (primera aproximación)

**Hipótesis.** El líquido no afecta al movimiento: la lata se comporta como un cilindro homogéneo, con el centro de masa en el eje. Rodadura sin deslizar. Envoltura con $I_c = m_c R^2$ ($\kappa = 1$).

**Ligaduras.** El cuerpo rígido plano tiene tres coordenadas $(x_O, y_O, \theta)$, con $O$ el eje. Hay dos ligaduras holónomas:

$$y_O = R, \qquad x_O = R\,\theta .$$

Queda **1 grado de libertad**, $\theta$.

**Energías.** Con $M = m_c + m_l$ y $\bar J = I_c + M R^2$:

$$T = \tfrac12 M \dot x_O^2 + \tfrac12 I_c \dot\theta^2 = \tfrac12 \bar J\,\dot\theta^2, \qquad V = \text{cte}.$$

**Lagrangiano y ecuación.**

$$L = \tfrac12 \bar J \dot\theta^2 \quad\Rightarrow\quad \bar J\,\ddot\theta = 0 .$$

$\theta$ es cíclica, así que $p_\theta = \bar J \dot\theta$ se conserva: rodadura uniforme, $\dot\theta = \Omega$.

**Observable.** Se mide el desplazamiento del punto marcado. Sea $\alpha_0$ su ángulo inicial medido desde la vertical hacia arriba, en el sentido de la rodadura ($\alpha_0 = \pi/2$ es hacia delante, en la horizontal). Con $\Delta\theta$ el giro desde el primer fotograma:

$$\Delta x = R\left[\Delta\theta + \sin(\alpha_0 + \Delta\theta) - \sin\alpha_0\right], \qquad \Delta y = R\left[\cos(\alpha_0 + \Delta\theta) - \cos\alpha_0\right].$$

Con $\dot\theta = \Omega$ constante es una **cicloide ideal**. La velocidad horizontal de la marca es $R\Omega\,[1 + \cos(\alpha_0 + \Delta\theta)]$: vale $2R\Omega$ arriba y $0$ en el contacto (el $2V_{CM}$ de la diapositiva de Physics 207).

**Verificado:** la marca alcanza $2R\Omega$ y $0$ con error relativo de ~10⁻⁶ (`verificar.py`, prueba 6).

**Qué falla.** No hay modulación de velocidad: no explica el fenómeno. Sirve como **referencia**: la diferencia entre la cicloide medida y la ideal es la huella del líquido. Además, en datos simulados el $R^2$ de la cicloide ideal puede ser alto (0,98 a 0,997 en rodadura rápida) aunque su error sea de 15 a 26 mm, así que se comparan **RMSE**, no solo $R^2$.

---

## Nivel 2. Líquido rígido pegado a la lata (un grado de libertad)

**Hipótesis.** El líquido es un sólido pegado a la lata: medio cilindro del mismo radio y longitud, con la cara plana pasando por el eje (llenado del 50 %). Es el límite de viscosidad muy alta.

**Ligaduras.** Las mismas del Nivel 1. Un grado de libertad. Ahora $\theta$ se mide desde la posición en que el CM del sistema está justo debajo del eje.

**Geometría.** El centroide de un semidisco está a $e = \dfrac{4R}{3\pi}$ del eje. El CM total queda a una distancia

$$a = \frac{m_l\,e}{M}, \qquad I_O = I_c + \tfrac12 m_l R^2$$

del eje, con $I_O$ la inercia respecto al eje geométrico.

**Energías.** El cuerpo rígido gira instantáneamente alrededor del punto de contacto. La distancia de este al CM es $d^2 = R^2 + a^2 - 2Ra\cos\theta$, y por Steiner:

$$J(\theta) = (I_O - M a^2) + M d^2 = I_O + M R^2 - 2 M R a\cos\theta,$$

$$T = \tfrac12 J(\theta)\,\dot\theta^2, \qquad V = -M g a\cos\theta .$$

**Euler-Lagrange.** Con $L = \tfrac12 J\dot\theta^2 + M g a\cos\theta$ y $J' = 2 M R a \sin\theta$:

$$\frac{d}{dt}(J\dot\theta) - \tfrac12 J'\dot\theta^2 + M g a\sin\theta = 0
\;\Rightarrow\;
J(\theta)\,\ddot\theta + M R a \sin\theta\,\dot\theta^2 + M g a \sin\theta = 0 .$$

$L$ no depende explícitamente de $t$ y la ligadura es esclerónoma, así que se conserva

$$E = \tfrac12 J(\theta)\dot\theta^2 - M g a\cos\theta .$$

**Fuerzas** (Newton, a posteriori): fricción estática $F = M[(R - a\cos\theta)\ddot\theta + a\sin\theta\,\dot\theta^2]$ y normal $N = M[g + a\cos\theta\,\dot\theta^2 + a\sin\theta\,\ddot\theta]$. Validez: $|F| \le \mu_s N$ (no desliza) y $N > 0$ (no despega).

**Régimen de vaivén (pequeñas oscilaciones).** Linealizando en $\theta = 0$:

$$\omega_0^2 = \frac{M g a}{J(0)} = \frac{M g a}{I_O + M R^2 - 2 M R a} .$$

Control analítico: con $m_c \to 0$ resulta $\omega_0^2 \approx 0{,}652\,g/R$, la fórmula conocida del semicilindro macizo que se mece (`lata_rodando.py --selftest`).

**Régimen de rodadura continua** (el de los videos). El CM pasa por arriba ($\theta = \pi$, donde $V = +Mga$) si $E > M g a$:

$$\omega_0^2 > \frac{2 M g a\,(1 + \cos\theta_0)}{J(\theta_0)} .$$

Con $\varepsilon = \dfrac{M R a}{I_O + M R^2}$ y $\Omega^2 \equiv \dfrac{2E}{I_O + M R^2}$, a primer orden en $\varepsilon$:

$$\dot\theta \approx \Omega\left[1 + \varepsilon\left(1 + \frac{g}{R\Omega^2}\right)\cos\theta\right] .$$

La velocidad es máxima cuando el CM pasa por abajo y mínima cuando pasa por arriba. Esa es la oscilación de velocidad del reto.

**Verificado** (`lata_rodando.py`, `notebooks/lata_tracker_colab_sencillo.ipynb`): el modelo coincide entre ambos códigos y con `modelos_progresivos.py` (prueba 7, diferencia 0); el umbral de rodadura se comprobó por simulación (por debajo se mece, por encima rueda sin parar); la fórmula de primer orden tiene error de ~3 % de la amplitud con $\varepsilon \approx 0{,}04$ y de ~16 % con $\varepsilon \approx 0{,}17$. Con $\varepsilon$ mayor, usar la simulación numérica.

**Qué falla.** Un líquido real mantiene la superficie libre casi horizontal y no gira como sólido: este modelo es una **cota superior** de la excentricidad y puede sobrestimar las oscilaciones. Para la lata del ensayo (12,2 g, 193 g de líquido, $R = 5{,}4$ cm) daría $\varepsilon \approx 0{,}26$, bien fuera del rango donde la fórmula de primer orden es confiable. Tampoco hay disipación y solo cubre el 50 %.

---

## Nivel 3. Líquido libre: modo de oleaje pendular (dos grados de libertad)

**Idea.** El líquido no gira con la lata. Su respuesta a las aceleraciones de la lata se describe con el **modo fundamental de oleaje**, cuyo equivalente mecánico es un péndulo:

- una masa efectiva $m_1 = \mu_s\, m_l$, colgada del eje con una varilla de longitud $l_s = R/\kappa_s$, cuyo ángulo $\varphi$ se mide **respecto a la vertical**;
- el resto del líquido, $m_l - m_1$, se traslada con el eje pero **no rota** (sin inercia rotacional propia).

$\kappa_s$ y $\mu_s$ dependen de la fracción de llenado $f$; el código las interpola de una tabla (para $f = 0{,}5$: $\kappa_s = 1{,}3557$ y $\mu_s = 0{,}458$). **El origen de esa tabla no está verificado**: hay que citar la fuente (oleaje lineal en un cilindro horizontal) antes de usarla en un informe.

**Ligaduras.** Rodadura sin deslizar, $x_O = R\theta$, $y_O = R$, y varilla rígida entre el eje y la masa del péndulo, $|\vec r_1 - \vec r_O| = l_s$. Las coordenadas de la masa del péndulo son

$$x_1 = R\theta + l_s\sin\varphi, \qquad y_1 = R - l_s\cos\varphi .$$

Quedan **2 grados de libertad**, $q = (\theta, \varphi)$.

**Energía cinética.** Suma de tres partes: la lata ($\tfrac12 m_c R^2\dot\theta^2 + \tfrac12 I_c\dot\theta^2$), el líquido que solo se traslada ($\tfrac12 (m_l - m_1) R^2\dot\theta^2$) y la masa del péndulo, con velocidad $(R\dot\theta + l_s\cos\varphi\,\dot\varphi,\; l_s\sin\varphi\,\dot\varphi)$:

$$T = \tfrac12 \bar J\,\dot\theta^2 + m_1 R\, l_s\cos\varphi\;\dot\theta\dot\varphi + \tfrac12 J_s\,\dot\varphi^2, \qquad \bar J = I_c + M R^2,\quad J_s = m_1 l_s^2 .$$

**Energía potencial.** Solo la masa del péndulo cambia de altura:

$$V = -m_1 g\, l_s\cos\varphi .$$

**Lagrangiano.** Con $M_{12}(\varphi) = m_1 R\, l_s\cos\varphi$:

$$L = \tfrac12\bar J\dot\theta^2 + M_{12}(\varphi)\,\dot\theta\dot\varphi + \tfrac12 J_s\dot\varphi^2 + m_1 g\, l_s\cos\varphi .$$

**Euler-Lagrange.**

*Para $\theta$* ($\partial L/\partial\theta = 0$):

$$\frac{d}{dt}\big(\bar J\dot\theta + M_{12}\dot\varphi\big) = 0
\;\Rightarrow\;
\bar J\,\ddot\theta + M_{12}\,\ddot\varphi = m_1 R\, l_s\sin\varphi\;\dot\varphi^2 .$$

*Para $\varphi$* ($\partial L/\partial\varphi = -m_1 R\, l_s \sin\varphi\,\dot\theta\dot\varphi - m_1 g\, l_s\sin\varphi$):

$$M_{12}\,\ddot\theta + J_s\,\ddot\varphi = -m_1 g\, l_s\sin\varphi .$$

(Los términos $\mp m_1 R\, l_s\sin\varphi\,\dot\theta\dot\varphi$ de ambos lados se cancelan.) En forma matricial:

$$\begin{pmatrix}\bar J & M_{12}\\ M_{12} & J_s\end{pmatrix}\begin{pmatrix}\ddot\theta\\ \ddot\varphi\end{pmatrix} = \begin{pmatrix} m_1 R\, l_s\sin\varphi\,\dot\varphi^2\\ -m_1 g\, l_s\sin\varphi\end{pmatrix}.$$

**Cantidades conservadas.** $\theta$ es **cíclica**, así que se conserva

$$p_\theta = \frac{\partial L}{\partial\dot\theta} = \bar J\dot\theta + M_{12}\dot\varphi .$$

Como $L$ no depende de $t$, también se conserva $E = T + V$. La coordenada $\varphi$ no es cíclica: el oleaje intercambia energía con la rotación a través del acoplamiento $M_{12}$.

**Hamilton.** Con $p_\theta = \bar J\dot\theta + M_{12}\dot\varphi$ y $p_\varphi = M_{12}\dot\theta + J_s\dot\varphi$, y $D = \bar J J_s - M_{12}^2$:

$$H = \frac{J_s\,p_\theta^2 - 2 M_{12}\,p_\theta p_\varphi + \bar J\,p_\varphi^2}{2D} - m_1 g\, l_s\cos\varphi ,$$

$$\dot\theta = \frac{J_s p_\theta - M_{12}p_\varphi}{D},\quad
\dot\varphi = \frac{\bar J p_\varphi - M_{12}p_\theta}{D},\quad
\dot p_\theta = 0,\quad
\dot p_\varphi = -\frac{\partial H}{\partial\varphi}.$$

**Pequeñas oscilaciones.** Alrededor de $\varphi = 0$ (equilibrio, con la lata en reposo o rodando a velocidad constante) los términos cuadráticos en $\dot\varphi$ se desprecian y $M_{12} \approx m_1 R\, l_s$. La ecuación característica $\det(K - \omega^2 M) = 0$, con $K = \mathrm{diag}(0,\, m_1 g\, l_s)$, da dos modos:

- $\omega = 0$: rodadura uniforme (traslación libre);
- el modo de oleaje acoplado

$$\omega_{ac}^2 = \frac{\omega_s^2}{1 - m_1 R^2/\bar J}, \qquad \omega_s^2 = \frac{g}{l_s} = \frac{\kappa_s\, g}{R}.$$

$\omega_s$ es la frecuencia con la lata fija; la lata que gira al compás del líquido la sube. Para el ensayo del cuaderno ($m_c = 12{,}2$ g, $m_l = 193$ g, $R = 5{,}4$ cm, $f = 50\,\%$): $m_1 = 88{,}4$ g, $l_s = 3{,}98$ cm, $\omega_s = 15{,}69$ rad/s (2,50 Hz) y $\omega_{ac} = 20{,}37$ rad/s (3,24 Hz).

**Límite.** Con $m_1 \to 0$ el modelo vuelve al Nivel 1.

**Verificado** (`verificar.py`): límite $m_1\to 0$ (error 10⁻¹⁰); energía y $p_\theta$ constantes (deriva ~10⁻¹¹); Hamilton y Euler-Lagrange dan la misma trayectoria (diferencia ~10⁻¹⁰); la frecuencia numérica coincide con $\omega_{ac}$ (error relativo 2·10⁻⁷).

**Qué falla.** No hay disipación: el oleaje nunca se apaga, y la viscosidad del líquido no interviene.

---

## Nivel 4. Disipación: rodadura y viscosidad (modelo del cuaderno final)

**Idea.** Se agregan dos pérdidas, ambas lineales en la velocidad, mediante la **función de disipación de Rayleigh**:

$$\mathcal F = \tfrac12 c_{rod}\,\dot\theta^2 + \tfrac12 c_{visc}\,\dot\varphi^2 .$$

- $c_{rod}$: resistencia de rodadura de la lata (par proporcional a $\dot\theta$).
- $c_{visc}$: amortiguamiento viscoso del oleaje.

**Ecuaciones.** Las fuerzas generalizadas son $Q_i = -\partial\mathcal F/\partial\dot q_i$, es decir $Q_\theta = -c_{rod}\dot\theta$ y $Q_\varphi = -c_{visc}\dot\varphi$:

$$\begin{pmatrix}\bar J & M_{12}\\ M_{12} & J_s\end{pmatrix}\begin{pmatrix}\ddot\theta\\ \ddot\varphi\end{pmatrix} = \begin{pmatrix} m_1 R\, l_s\sin\varphi\,\dot\varphi^2 - c_{rod}\dot\theta\\ -m_1 g\, l_s\sin\varphi - c_{visc}\dot\varphi\end{pmatrix}.$$

La energía ya no se conserva: $\dot E = -2\mathcal F = -c_{rod}\dot\theta^2 - c_{visc}\dot\varphi^2 \le 0$. $p_\theta$ tampoco, porque $c_{rod}$ rompe la simetría.

**Relación con la viscosidad del líquido.** Un modo lineal con la lata fija cumple $J_s\ddot\varphi + c_{visc}\dot\varphi + m_1 g\, l_s\varphi = 0$, es decir $\ddot\varphi + 2\zeta\omega_s\dot\varphi + \omega_s^2\varphi = 0$ con

$$c_{visc} = 2\,\zeta\, J_s\,\omega_s .$$

El código toma como valor inicial $\zeta = \min\!\big(0{,}85,\; 2\sqrt{\nu/(\omega_s R^2)}\big)$. El argumento $\nu/(\omega_s R^2)$ es adimensional (viscosidad cinemática sobre frecuencia por radio al cuadrado), y la dependencia $\zeta\propto\sqrt{\nu/(\omega_s R^2)}$ es la escala típica del amortiguamiento por capa límite laminar. **El coeficiente 2 no está respaldado por una derivación ni una referencia verificadas**: es un valor inicial, y el ajuste con los datos decide el $c_{visc}$ final.

Valores iniciales para la lata del ensayo ($\zeta$ con tope en 0,85, $c_{visc}$ en N·m·s):

| Líquido | $\nu$ (m²/s) | $\zeta$ | $c_{visc}$ inicial |
|---|---|---|---|
| Agua | 1,004·10⁻⁶ | 0,0094 | 4,1·10⁻⁵ |
| Aceite vegetal | 5,5·10⁻⁵ | 0,069 | 3,1·10⁻⁴ |
| Aceite de motor | 2,4·10⁻⁴ | 0,145 | 6,4·10⁻⁴ |
| Glicerina | 1,18·10⁻³ | 0,321 | 1,4·10⁻³ |

**Verificado** (`verificar.py`): el balance de energía $E(0) - E(T) = \int (c_{rod}\dot\theta^2 + c_{visc}\dot\varphi^2)\,dt$ se cumple con error relativo de 4·10⁻⁸, y la energía no crece nunca.

---

## Del modelo a la predicción (resumen)

| Nivel | Grados de libertad | Qué explica | Qué no explica |
|---|---|---|---|
| 1 | 1 ($\theta$) | Cicloide ideal; referencia | Ninguna modulación de velocidad |
| 2 | 1 ($\theta$) | Modulación de velocidad por el CM excéntrico; vaivén y rodadura continua | Líquido que no gira; disipación; solo 50 % |
| 3 | 2 ($\theta, \varphi$) | Líquido libre con modo de oleaje; frecuencia acoplada; llenado variable | Disipación |
| 4 | 2 ($\theta, \varphi$) | Amortiguamiento por rodadura y por viscosidad | Oleaje no lineal, superficie libre real, otros modos |

**Interpretación física de las oscilaciones de velocidad.** En todos los niveles con líquido, la rotación de la lata intercambia energía con una masa que sube y baja respecto al eje. En el Nivel 2 esa masa va pegada a la lata y el intercambio es potencial. En los niveles 3 y 4 el líquido oscila libremente y el intercambio pasa por el acoplamiento $M_{12}(\varphi)$ entre la rotación y el oleaje. La amplitud crece con la fracción de masa que oscila ($m_1/M$) y baja al aumentar la velocidad de giro.

## Aviso sobre el ajuste con datos

El Nivel 4 tiene cinco parámetros ajustables ($\dot\theta_0$, $\varphi_0$, $\dot\varphi_0$, $c_{visc}$, $c_{rod}$), más el ángulo inicial de la marca $\alpha_0$. En una prueba con **datos simulados por el propio modelo** (ruido de 0,2 mm, 0,65 s de video, calibración con el primer 45 %), el cuaderno final **no recuperó los parámetros** (RMSE de predicción de 55 mm). La causa y las pruebas están en `04_hoja_de_ruta.md`. Hasta resolverlo, conviene leer los parámetros ajustados con cautela y comparar siempre con los niveles 1 y 2.
