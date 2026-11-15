# Reto IV-2: lata parcialmente llena que rueda

Trabajo final · modelo matemático de una lata con líquido que rueda, contrastado con videos analizados en Tracker.

## Resumen

**Problema.** Un recipiente cilíndrico parcialmente lleno de líquido que rueda presenta oscilaciones en su velocidad.

**Objetivo.** Predecir el movimiento con un modelo matemático y contrastarlo con videos analizados en Tracker (ya no basta con la explicación cualitativa).

**Sistema.** Una lata con un medio cilindro de líquido del mismo radio `R` y la misma longitud `L` que la lata, con la cara plana pasando por el eje: llenado del 50 %.

**Modelo vigente.** Un grado de libertad (el ángulo de giro `θ`). El líquido se trata como un sólido pegado a la lata. El centro de masa queda desplazado del eje, y eso produce las oscilaciones: al rodar, el líquido sube y baja; cuando sube la lata frena y cuando baja acelera.

**Preguntas del reto (respuesta provisional).**

| Pregunta | Respuesta del modelo actual |
|---|---|
| ¿Qué produce las oscilaciones? | La altura del CM cambia con el giro: modula la energía potencial y la inercia efectiva. |
| ¿De qué depende la amplitud? | Crece con `ε` (excentricidad y masa de líquido) y disminuye al aumentar la velocidad de giro `Ω`. |
| ¿Cómo varía con el nivel de llenado? | El modelo está fijo al 50 %. Otros niveles requieren el segmento circular (ver `docs/modelo.md`, sec. 6). |
| ¿Cómo interviene la viscosidad? | Todavía no interviene; entra en la Etapa IV, solo si el contraste con datos lo exige. |

## Estado por etapas

| Etapa | Contenido | Estado |
|---|---|---|
| I | Geometría: masa, centro de masa, inercia | Formulada (medio cilindro, 50 %) |
| II | Mecánica: Lagrange, rodadura sin deslizar, fuerzas | Formulada y verificada numéricamente |
| III | Pequeñas oscilaciones y frecuencia natural | Formulada |
| IV | Viscosidad y movimiento interno del líquido | Pendiente (solo si el contraste lo exige) |
| V | Conexión con mediciones de Tracker | Protocolo escrito, sin datos reales |
| VI | Simulación y comparación con experimento | Script listo; sin contraste real aún |

## Qué está verificado y qué no

**Verificado**
- `python lata_rodando.py --selftest` pasa las 7 comprobaciones: `ω₀² = 0,652·g/R` con `m_c → 0` (semicilindro macizo), rodadura uniforme con `a → 0`, conservación de la energía (vaivén grande y rodadura continua), solución lineal frente a numérica, periodo numérico frente a `2π/ω₀`, y coherencia Lagrange ↔ Newton-Euler.
- La aproximación de primer orden de la rodadura continua (`docs/modelo.md`, sec. 5) se contrastó con la integración numérica: error máximo de ~3 % de la amplitud pico a pico con `ε ≈ 0,04` y de ~16 % con `ε ≈ 0,17`.

**No verificado**
- El script **nunca se ha corrido con un CSV real de Tracker**; solo con datos sintéticos generados por el propio modelo.
- Las fórmulas del segmento circular (otros niveles de llenado) son derivaciones propias; faltan contrastarlas con una referencia estándar de centroides y momentos de áreas.
- No hay datos experimentales ni referencias bibliográficas verificadas en este repositorio.

## Uso

```bash
pip install -r requirements.txt

# 1. Verificar el entorno (sin datos)
python lata_rodando.py --selftest

# 2. Simular con parámetros medidos
python lata_rodando.py --R <m> --mc <kg> --ml <kg> --kappa <k> \
       --theta0 <rad> --omega0 <rad/s> --tmax <s>

# 3. Contrastar con Tracker (CSV/TSV con columnas t y x del centro)
python lata_rodando.py --R <m> --mc <kg> --ml <kg> --kappa <k> \
       --csv data/raw/datos.csv --tcol 0 --xcol 1 --theta0 0 \
       --out results/experimento1.png
# Añadir --fit theta0,omega0 para refinar solo las condiciones iniciales,
# o --fit theta0,omega0,kappa si además se quiere estimar kappa.
```

Salida: curvas `x(t)`, `v(t)`, `θ(t)`, periodo, régimen (vaivén o rodadura continua), `μ_s` mínimo y `N` mínimo (validez de no deslizar y no despegar), métricas (RMSE, error máximo, NRMSE, R²), frecuencia dominante por FFT y figura de residuos.

## Estructura

```
trabajo final/
├── README.md
├── requirements.txt
├── lata_rodando.py          # modelo, simulación, lectura de Tracker, ajuste, autotest, CLI
├── docs/
│   ├── modelo.md            # hipótesis, ecuaciones, oscilaciones, otros llenados
│   ├── experimento_tracker.md  # variables, protocolo, cuidados con los datos
│   └── hoja_de_ruta.md      # lo que falta, criterio para la Etapa IV, preguntas abiertas
├── data/raw/                # exportaciones de Tracker sin tocar
├── data/processed/
└── results/                 # figuras y tablas de cada experimento
```

## Reglas del trabajo

1. Trabajar por etapas: primero hipótesis, luego desarrollo matemático, al final interpretación física.
2. Complejidad progresiva: viscosidad o grados de libertad internos solo cuando el contraste con datos lo justifique, y con justificación física explícita.
3. Distinguir siempre parámetros **medidos**, **calculados** y **ajustados**. Masas y geometría nunca se ajustan (el modelo debe seguir siendo predictivo); solo condiciones iniciales y `κ`.
4. No inventar datos experimentales ni bibliografía.
5. Verificar cada fórmula nueva con el autotest o con un caso límite antes de usarla.
6. Guardar con cada experimento un archivo de parámetros (medidas, incertidumbres, fps, fecha) junto al CSV de Tracker, y no modificar nunca los datos crudos.
