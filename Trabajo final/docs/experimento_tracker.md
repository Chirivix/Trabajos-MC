# Protocolo experimental y mediciones en Tracker

## 1. Variables

| Variable | Cómo se obtiene | Tipo |
|---|---|---|
| `R`, `R_int`, `L` | Calibre | Medida directa |
| `m_c`, `m_l` | Balanza (lata vacía y lata llena; `m_l` = llena − vacía) | Medida directa |
| Escala del video | Una longitud conocida visible en el cuadro (p. ej. el diámetro) | Calibración |
| `t` | Tiempo del video (fps conocidos) | Medida directa |
| `x(t)` del centro de la lata | Tracker, un punto en el centro geométrico | Medida directa |
| `v(t)` | Derivada numérica de `x(t)`, suavizada | Calculada |
| `θ(t)` | Opcional: rotación de una marca en la lata | Medida o calculada |
| `e`, `a`, `I_O`, `I_cm`, `ω₀`, `ε` | Fórmulas de `modelo.md` | Calculadas |
| `κ` | Geometría de la lata, o ajuste | Estimada o ajustada |
| `θ₀`, `ω₀` | Condiciones iniciales | Medidas o ajustadas (se declara cuál) |

Se ajustan únicamente las condiciones iniciales y `κ`. Masas y geometría nunca se ajustan.

## 2. Primer experimento

1. Colocar la lata sobre una superficie horizontal y esperar a que quede en reposo. Se asienta con el líquido abajo, así que `θ₀ = 0` sin necesidad de verlo.
2. Dar un empujón suave para que se meza o ruede despacio.
3. Grabar de lado, con la cámara fija y perpendicular al movimiento.
4. En Tracker, calibrar la escala y marcar el centro de la lata en cada cuadro.
5. Exportar tiempo y `x` a `data/raw/`.
6. Correr el script con las medidas (ver README).

## 3. Primera comparación

El resultado más limpio es el **periodo de balanceo**: se mide en el video y el modelo lo predice solo con `R`, `m_c`, `m_l` y `κ` (`T₀ = 2π/ω₀`). Se compara primero este valor, antes de mirar nada más fino.

## 4. Cómo leer el resultado

- Coinciden periodo y amplitud: el modelo simple explica el fenómeno.
- La frecuencia coincide pero la amplitud decae en el video y no en el modelo: falta disipación.
- La frecuencia no coincide: probablemente el líquido no se mueve pegado a la lata (sloshing).

## 5. Cuidados con los datos

- Derivar datos ruidosos amplifica el error: la ventana de suavizado debe cubrir bastante menos que un periodo.
- Conviene tener decenas de muestras por periodo (el script informa cuántas hay).
- La velocidad inicial estimada con pocos puntos es imprecisa: en una prueba sintética con ruido dio 2,29 rad/s cuando el valor real era 2,5.
- Registrar la incertidumbre de cada medida (calibre, balanza, escala del video).
- Guardar con cada experimento un archivo de parámetros (medidas, incertidumbres, fps, fecha). No modificar nunca los datos crudos.
