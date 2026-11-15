# Modelo matemático (Etapas I–III)

## 1. Geometría y decisiones

- **Lectura de la geometría.** «Cilindro + semicilindro» se interpreta como una composición: un medio cilindro de líquido dentro de la lata, con las mismas medidas que ella (no es suma ni resta de masas). Un segundo medio cilindro encima la dejaría llena. Esto fija el llenado en el 50 % y no deja parámetros geométricos libres. *(Pendiente de confirmar con el profesor.)*
- **Formalismo.** Lagrange: un solo grado de libertad y la fricción estática no hace trabajo. Las fuerzas se recuperan después con Newton.

## 2. Hipótesis

- **H1.** El líquido es rígido y solidario a la lata (límite de viscosidad muy alta).
- **H2.** Movimiento plano.
- **H3.** Rodadura sin deslizamiento sobre un plano horizontal.
- **H4.** Sin resistencia de rodadura ni disipación.
- **H5.** La lata no rebota.

## 3. Definiciones y magnitudes

- `R`: radio exterior (el que rueda). `R_int`: radio interior.
- `m_c`: masa de la lata vacía. `m_l`: masa del líquido.
- `κ`: inercia de la envoltura, `I_c = κ·m_c·R²` (entre 0,5 y 1; 1 si fuera solo pared lateral, baja un poco por las tapas). Es el único parámetro de la lata que no se mide directamente.
- `θ`: ángulo girado en el sentido de la rodadura. `θ = 0` cuando el CM del sistema está justo debajo del eje (líquido abajo, cara plana arriba, posición de reposo natural). Se mide desde el líquido, no desde una marca cualquiera.
- Rodadura sin deslizar: `x_O = R·θ`.

```
e   = 4·R_int / (3π)                    # centroide del semidisco de líquido respecto al eje
M   = m_c + m_l
a   = m_l·e / M                         # excentricidad del CM total
I_O = κ·m_c·R² + ½·m_l·R_int²           # inercia respecto al eje geométrico
I_cm = I_O − M·a²                       # inercia respecto al CM (Steiner)
J(θ) = I_O + M·R² − 2·M·R·a·cos θ      # inercia respecto al punto de contacto
```

Si `R_int ≠ R`: `R` exterior en la condición de rodadura y en la inercia de la envoltura; `R_int` para el líquido.

## 4. Ecuación de movimiento

`T = ½·J(θ)·θ'²`, `V = −M·g·a·cos θ`. Lagrange:

```
J(θ)·θ'' + M·R·a·sin θ·θ'² + M·g·a·sin θ = 0
```

Se conserva `E = ½·J·θ'² − M·g·a·cos θ`. Si `a = 0` (lata homogénea), `θ'' = 0`: rodadura uniforme.

Fuerzas (Newton, a posteriori):

```
F = M·[(R − a·cos θ)·θ'' + a·sin θ·θ'²]     # fricción estática
N = M·[g + a·cos θ·θ'² + a·sin θ·θ'']       # normal
```

Validez: `|F| ≤ μ_s·N` (no desliza) y `N > 0` (no despega).

## 5. Dos regímenes

**Pequeñas oscilaciones (vaivén).** Alrededor de `θ = 0`:

```
θ'' + ω₀²·θ = 0,     ω₀² = M·g·a / (I_O + M·R² − 2·M·R·a)
θ(t) = θ₀·cos(ω₀t) + (θ'₀/ω₀)·sin(ω₀t),     v_O = R·θ'
```

Control analítico: con `m_c → 0`, `ω₀² ≈ 0,652·g/R`, que coincide con la fórmula conocida del semicilindro macizo que se mece. Validez: `θ_max = θ'₀/ω₀ ≪ 1`. Quedan fuera amplitud grande, deslizamiento, despegue y disipación.

**Rodadura continua (el CM pasa por arriba).** Ocurre si la energía inicial supera `M·g·a`. Con `ε = M·R·a / (I_O + M·R²)` y `Ω² ≡ 2E/(I_O + M·R²)`:

```
θ' ≈ Ω·[1 + ε·(1 + g/(R·Ω²))·cos θ]
```

La velocidad es máxima cuando el CM pasa por abajo y mínima cuando pasa por arriba.

*Verificación numérica (R = g = 1):* el error máximo de esta fórmula frente a la integración es ~3 % de la amplitud pico a pico con `ε ≈ 0,04` y ~16 % con `ε ≈ 0,17`. Es una aproximación de primer orden en `ε`: para latas con `ε` grande conviene usar la simulación numérica y no la fórmula.

## 6. Otros niveles de llenado (no implementado)

Las ecuaciones de las secciones 3 a 5 valen igual; cambian `m_l`, la distancia `d` del centroide y la inercia del líquido. Con semiángulo central `θ_s` y fracción de llenado `φ_f`:

```
φ_f = (θ_s − sin θ_s·cos θ_s) / π
d   = 2·R·sin³θ_s / [3·(θ_s − sin θ_s·cos θ_s)]     # centroide respecto al centro
```

Casos límite que se cumplen: `θ_s = π/2` da `φ_f = ½` y `d = 4R/(3π)` (el modelo actual); `θ_s → π` da lata llena y `d → 0`. La inercia del líquido para un segmento circular también hay que recalcularla (ya no es `½·m_l·R²`).

**Pendiente:** son derivaciones propias; contrastarlas con tablas de centroides y momentos de áreas de un texto de estática (p. ej. Beer & Johnston o Meriam & Kraige, comprobando la edición usada). No se citan páginas ni ediciones porque no se han verificado.

## 7. Limitaciones del modelo

- **Líquido rígido.** Un líquido real (agua) mantiene la superficie libre horizontal y no gira como sólido. El modelo es una cota superior de la excentricidad (límite de viscosidad muy alta); en el límite opuesto, líquido completamente libre, casi no hay oscilación. Puede sobrestimar las oscilaciones.
- **Sin disipación.** La amplitud medida decaerá y los residuos crecerán con el tiempo; es esperable y motiva la Etapa IV.
- **Llenado fijo** al 50 %.
- **Deslizamiento y despegue.** Si el script reporta `N ≤ 0` o un `μ_s` mínimo mayor que el real, el modelo no es válido en esas condiciones.
- **Plano horizontal.** Un plano inclinado añadiría el término de la pendiente.

## 8. Sobre el código (`lata_rodando.py`)

- Clase `Lata`: parámetros y magnitudes derivadas (`e`, `a`, `I_O`, `I_cm`, `ε`, `ω₀`).
- `simular`: integra con DOP853 (Runge-Kutta de orden 8, paso adaptativo, `rtol = 1e-10`); devuelve `θ`, `ω`, `x_O`, `v_O`, posición del CM, `N`, `F` y energía. Si en la Etapa IV aparece viscosidad con escalas de tiempo muy distintas, habrá que revisar el método (Radau o LSODA).
- `solucion_lineal`, `regimen` (vaivén o rodadura continua por energía), `leer_tracker` (CSV/TSV, admite coma decimal), `velocidad_suavizada` (Savitzky-Golay), `comparar` (RMSE, error máximo, NRMSE, R², FFT, `μ_s` mínimo, residuos), `ajustar` (mínimos cuadrados solo sobre `θ₀`, `ω₀` y `κ`), `autotest`.
