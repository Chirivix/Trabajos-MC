# Reto IV-2: lata parcialmente llena que rueda

Trabajo final · de una primera aproximación al modelo con oleaje y viscosidad, contrastado con videos analizados en Tracker.

## El problema

Un recipiente cilíndrico parcialmente lleno de líquido que rueda presenta oscilaciones en su velocidad. El objetivo es **predecir el movimiento con un modelo mecánico** y contrastarlo con videos de Tracker, no solo explicarlo cualitativamente. En los videos la lata **rueda de forma continua** (no son pequeñas oscilaciones): un punto marcado en el borde describe una cicloide, y el líquido la deforma.

Preguntas del reto: ¿qué produce las oscilaciones de velocidad?, ¿de qué depende su amplitud?, ¿cómo varían con el nivel de llenado?, ¿cómo interviene la viscosidad?

## El recorrido: cuatro niveles

Cada nivel agrega **una** idea física, se resuelve con la misma receta (ligaduras → $T$ y $V$ → lagrangiano → Euler-Lagrange o Hamilton → código) y se comprueba antes de pasar al siguiente. El desarrollo completo, como problema de parcial, está en [`docs/01_problema_parcial.md`](docs/01_problema_parcial.md).

| Nivel | Idea física que se agrega | Grados de libertad | Qué aprendimos o por qué se pasó al siguiente |
|---|---|---|---|
| 1 | Lata homogénea: rodadura sin deslizar, $L = \tfrac12\bar J\dot\theta^2$ | 1 ($\theta$) | Da una cicloide ideal, sin modulación de velocidad. Sirve de **referencia**: la diferencia con el video es la huella del líquido. |
| 2 | Líquido **rígido** pegado a la lata: CM excéntrico, $L = \tfrac12 J(\theta)\dot\theta^2 + Mga\cos\theta$ | 1 ($\theta$) | Explica la modulación de velocidad, el vaivén y la rodadura continua (con umbral de energía). Pero un líquido real no gira como sólido: es una cota superior. |
| 3 | Líquido **libre**: modo de oleaje como péndulo de masa $m_1$ y longitud $l_s$ | 2 ($\theta$, $\varphi$) | $\theta$ cíclica ⇒ $p_\theta$ conservado; acoplamiento $M_{12}(\varphi)$; frecuencia del oleaje acoplado. No hay disipación. |
| 4 | **Disipación** de Rayleigh: rodadura ($c_{rod}$) y viscosidad ($c_{visc}$, ligada a $\nu$) | 2 ($\theta$, $\varphi$) | Es el modelo del **cuaderno final**. Ecuaciones verificadas; el ajuste con datos aún tiene limitaciones (ver abajo). |

Hallazgos del camino:

- **El $R^2$ engaña.** Una cicloide ideal ya tiene $R^2$ alto cuando la lata rueda rápido; se comparan RMSE, no solo $R^2$.
- **Orden de los regímenes.** Si el empujón es débil, la lata solo se mece (pequeñas oscilaciones); para dar vueltas hace falta superar un umbral de energía.
- **El sentido de $x$** no se puede tomar de $x_{final} - x_{inicial}$ cuando la lata se mece; el cuaderno sencillo usa la pendiente inicial.
- **El ángulo inicial de la marca** no es un detalle: un error de ~10° produce errores de varios milímetros (ver `docs/04_hoja_de_ruta.md`).

## Cómo se pasa a código

[`docs/02_del_modelo_al_codigo.md`](docs/02_del_modelo_al_codigo.md) muestra el paso de las ecuaciones al código nivel por nivel: de Euler-Lagrange a $M(q)\ddot q = \text{rhs}$, del vector de estado a `solve_ivp`, de Hamilton a código, y del cuaderno final celda por celda.

## Contenido

```
Trabajo final/
├── README.md
├── requirements.txt
├── lata_rodando.py                       # Nivel 2 de referencia: autotest, fuerzas N y F, ajuste de kappa
├── modelos/
│   ├── modelos_progresivos.py            # Niveles 1 a 4 y Hamilton, con la misma interfaz
│   └── verificar.py                      # pruebas de cada nivel
├── notebooks/
│   ├── lata_tracker_colab_sencillo.ipynb # Nivel 2 con datos de Tracker (marca del borde)
│   └── Mecanica_clasica2.ipynb           # Nivel 4: modelo final con oleaje y viscosidad
├── docs/
│   ├── 01_problema_parcial.md            # el problema resuelto, nivel por nivel
│   ├── 02_del_modelo_al_codigo.md        # de las ecuaciones al código
│   ├── 03_experimento_tracker.md         # protocolo de medición
│   └── 04_hoja_de_ruta.md                # limitaciones, próximos pasos, preguntas abiertas
├── data/raw/ , data/processed/           # exportaciones de Tracker sin tocar / procesadas
└── results/                              # figuras y tablas de cada experimento
```

## Uso

Requiere Python 3.10 o superior. `pip install -r requirements.txt`.

```bash
python modelos/verificar.py        # comprueba los niveles 1 a 4 (y el Nivel 2 contra lata_rodando.py)
python lata_rodando.py --selftest  # autotest del Nivel 2
```

**Modelo final con datos de Tracker.** Abrir en Colab [`Mecanica_clasica2.ipynb`](https://colab.research.google.com/github/Chirivix/Trabajos-MC/blob/main/Trabajo%20final/notebooks/Mecanica_clasica2.ipynb), escribir los datos del ensayo en la celda 1 (masas, radio, llenado, FPS, líquido) y subir en la celda 2 el archivo de Tracker con el tiempo y la posición $x$, $y$ del punto marcado. Protocolo completo en [`docs/03_experimento_tracker.md`](docs/03_experimento_tracker.md).

**Comparación con el Nivel 2.** [`lata_tracker_colab_sencillo.ipynb`](https://colab.research.google.com/github/Chirivix/Trabajos-MC/blob/main/Trabajo%20final/notebooks/lata_tracker_colab_sencillo.ipynb) ajusta el modelo de líquido rígido y, en su celda 3, lo compara con una lata sin líquido.

**Modelos por niveles en Python.**

```python
import numpy as np
from modelos.modelos_progresivos import Lata, simular_nivel1, simular_nivel2, simular_oleaje, posicion_marca
t = np.linspace(0, 0.65, 157)                                      # tiempos en s
L = Lata(R=0.054, m_c=0.0122, m_l=0.193, fill=0.5)             # lata del ensayo del cuaderno final
th, om, ph, dph = simular_oleaje(L, omega0=9.0, phi0=0.2, dphi0=0.0, t=t,
                                 c_visc=L.c_visc_teo, c_rod=1e-6)  # Nivel 4
dx, dy = posicion_marca(L.R, th, alpha0=1.2)                       # lo que mide la cámara
```

## Qué está verificado

- Nivel 3 contiene al Nivel 1: con $m_1\to0$, $\theta = \omega_0 t$ (error de 10⁻¹⁰).
- La energía y $p_\theta$ se conservan sin disipación (deriva ~10⁻¹¹).
- Hamilton y Euler-Lagrange dan la misma trayectoria (diferencia ~10⁻¹⁰).
- La frecuencia del oleaje acoplado coincide con $\omega_{ac}$ (error relativo 2·10⁻⁷).
- Con disipación, el balance de energía cierra (error 4·10⁻⁸) y la energía no crece.
- La marca alcanza $2R\Omega$ arriba y $0$ en el contacto (cicloide ideal).
- El Nivel 2 coincide con `lata_rodando.py`.
- El cuaderno final, con los parámetros verdaderos, reproduce datos simulados con 0,19 mm de error (el nivel del ruido): las ecuaciones y la marca están bien implementadas.

**Limitaciones conocidas del ajuste del cuaderno final** (detalle en `docs/04_hoja_de_ruta.md`). Con datos simulados por el propio modelo, el cuaderno **no recuperó los parámetros** (RMSE de predicción de 55 mm, $R^2$ = 0,31): el ángulo inicial de la marca sale sesgado (~10°) y los amortiguamientos no quedan determinados con 0,3 s de calibración. Es una sola configuración simulada, pero obliga a leer los parámetros ajustados con cautela.

## Reglas del trabajo

1. Trabajar por niveles: ligaduras, $T$ y $V$, lagrangiano, ecuaciones, código. No pasar al siguiente hasta que el actual cumpla sus verificaciones.
2. Cada nivel debe contener al anterior como caso límite.
3. Distinguir siempre parámetros **medidos**, **calculados** y **ajustados**. Masas, radio y llenado nunca se ajustan.
4. No inventar datos experimentales ni bibliografía.
5. Guardar con cada experimento un archivo de parámetros (medidas, incertidumbres, FPS, fecha) junto al CSV de Tracker, y no modificar nunca los datos crudos.
