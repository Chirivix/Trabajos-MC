# Hoja de ruta y limitaciones conocidas

## Estado

| Nivel | Modelo | Código | Estado |
|---|---|---|---|
| 1 | Lata homogénea (cicloide ideal) | `modelos_progresivos.py` | Verificado con datos simulados |
| 2 | Líquido rígido, 1 GDL | `lata_rodando.py`, cuaderno sencillo, `modelos_progresivos.py` | Verificado con datos simulados |
| 3 | Oleaje pendular, 2 GDL | `modelos_progresivos.py` | Verificado (energía, $p_\theta$, Hamilton, frecuencia) |
| 4 | Oleaje + disipación | `modelos_progresivos.py`, cuaderno final | Ecuaciones verificadas; **ajuste con datos no validado** |

Falta, para todos: correr con videos reales de Tracker.

## Limitaciones del cuaderno final, probadas con datos simulados

Se generaron datos con el propio modelo del Nivel 4 (lata de 12,2 g con 193 g de agua, $R = 5{,}4$ cm, 240 fps, 0,65 s, ruido de 0,2 mm en $x$ e $y$, marca inicial a 70°) y se ejecutó el cuaderno final tal como está. Es **una sola configuración simulada**: sirve como alerta, no como medida general.

| Prueba | Resultado |
|---|---|
| Parámetros verdaderos y $\alpha_0$ verdadero | RMSE de predicción 0,19 mm (el ruido): el modelo y la marca están bien implementados |
| Ajuste completo del cuaderno | RMSE 54,8 mm, $R^2$ = 0,31; no recupera los parámetros |
| Parámetros verdaderos pero $\alpha_0$ del cuaderno | RMSE 5,6 mm |

### Problema A: el ángulo inicial de la marca sale sesgado

El cuaderno estima $\alpha_0$ ajustando una cicloide **ideal** a la trayectoria $(x, y)$, pero la trayectoria real está deformada por el líquido. En la prueba dio 79,8° en lugar de 70°, y ese error de ~10° propaga un error de 5,6 mm aun con los demás parámetros exactos. Ajustando $\alpha_0$ **junto con los demás** (seis parámetros), la prueba recuperó 70,0°.

### Problema B: los amortiguamientos y las condiciones iniciales del oleaje no están bien determinados

Con $\alpha_0$ correcto, el ajuste terminó con RMSE de predicción de 3,5 mm (frente a 0,19 mm con los parámetros verdaderos). Los 36 puntos iniciales probados convergieron a un costo (error de ajuste en la ventana de calibración) menor que el doble del de los parámetros verdaderos, y el mejor fue incluso ligeramente menor que ese, pero con valores muy distintos:

| Parámetro | Verdadero | Ajustado |
|---|---|---|
| $\varphi_0$ | 14,3° | 11,4° |
| $\dot\varphi_0$ | −2,0 rad/s | −1,7 rad/s |
| $c_{visc}$ | 3·10⁻⁴ N·m·s | 3,2·10⁻⁷ (en la cota inferior) |
| $c_{rod}$ | 2·10⁻⁶ N·m·s | 3,2·10⁻⁹ (en la cota inferior) |

Interpretación: en los 0,3 s de calibración el amortiguamiento casi no deja huella (sus efectos se acumulan después), así que muchas combinaciones ajustan igual de bien y se separan recién en el tramo de predicción. Es una limitación de **identificabilidad**, no un error de código.

### Qué hacer al respecto

1. Ajustar $\alpha_0$ junto con los demás parámetros, o fijarlo con una medida directa (el ángulo de la marca en el primer fotograma).
2. Usar una ventana de calibración más larga, o ajustar con varios videos a la vez.
3. Si los datos no determinan los amortiguamientos, **fijar uno** (por ejemplo $c_{rod} = 0$, o $c_{visc}$ de la teoría) y ajustar solo las condiciones iniciales; después liberarlos de uno en uno.
4. Imprimir una alerta cuando un parámetro quede en una cota.
5. Comparar siempre con los niveles 1 y 2 sobre el mismo video.

Estas propuestas **no se han probado todas**: solo la primera se probó (recupera $\alpha_0$), y la identificabilidad del resto sigue pendiente.

## Otras limitaciones

- **Tabla de $\kappa_s$ y $\mu_s$.** Se tomó tal cual del código final; su fuente no está verificada. Citar la referencia (oleaje lineal en un cilindro horizontal) o validarla contra experimentos.
- **Coeficiente 2 en $\zeta = 2\sqrt{\nu/(\omega_s R^2)}$.** Sin derivación ni referencia verificadas; es un valor inicial, no una predicción.
- **Líquido que no gira.** El modelo supone que el resto del líquido ($m_l - m_1$) se traslada sin rotar. Es un límite; un líquido viscoso real gira parcialmente con la lata.
- **Oleaje lineal.** Un solo modo y péndulo simple: no hay oleaje no lineal ni otros modos.
- **Sentido de $x$.** El cuaderno final lo toma de $x_{final} - x_{inicial}$; en vaivén puede fallar (ver el cuaderno sencillo, que usa la pendiente inicial).
- **Sin $y$ en el archivo**, $\alpha_0$ se fija en 90°.
- **Nivel 2 solo cubre el 50 %** (medio cilindro). Los niveles 3 y 4 cubren 5 a 95 % mediante la tabla.
- **Plano horizontal.** Un plano inclinado añadiría el término de la pendiente.
- **Deslizamiento y despegue** no se comprueban en los niveles 3 y 4 (las fuerzas $N$ y $F$ del Nivel 2 están en `lata_rodando.py`).

## Próximos pasos

1. Correr el cuaderno final con un video real y guardar parámetros, RMSE y gráficas en `results/`.
2. Resolver los problemas A y B (arriba) y volver a verificar con datos simulados antes de volver a datos reales.
3. Ajustar con un video y predecir otro distinto (otras condiciones iniciales).
4. Repetir con otros líquidos (el cuaderno trae agua, aceite vegetal, aceite de motor y glicerina) y comprobar si $c_{visc}$ ajustado crece con $\nu$.
5. Repetir con otros niveles de llenado.
6. Confirmar con el profesor la lectura de la geometría del Nivel 2 (medio cilindro de líquido, mismas medidas que la lata).

## Preguntas abiertas

1. ¿Superficie horizontal o plano inclinado? (los modelos suponen horizontal)
2. ¿Se confirmó con el profesor la geometría del Nivel 2?
3. ¿Cuál es el llenado real de las latas usadas?
4. ¿De dónde sale la tabla de $\kappa_s$ y $\mu_s$?
5. FPS y resolución de la cámara en cada video: *(anotar los valores reales)*.
6. ¿La lata es transparente o hay forma de observar el movimiento del líquido? Permitiría contrastar el ángulo $\varphi$ del modelo.

## Si se entrega a otra IA

Retomar así: leer el README y `01_problema_parcial.md`; correr `python modelos/verificar.py`; revisar las limitaciones de este documento; con el primer CSV real, correr el cuaderno final y comparar con los niveles 1 y 2. Reglas: las del README; no avanzar de nivel hasta que el actual cumpla sus verificaciones.
