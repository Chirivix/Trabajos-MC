# Hoja de ruta

## Antes de medir
- Responder las preguntas abiertas (sección final).
- Opcional: dividir el script en módulos (geometría, dinámica, lectura de Tracker, comparación) y pasar el autotest a pruebas formales (`tests/`).

## Etapa V: medir y conectar con Tracker
1. Medir `R`, `L`, `m_c` y `m_l`; comprobar que el llenado real esté cerca del 50 %.
2. Grabar y analizar el movimiento (`experimento_tracker.md`).
3. Comparar primero el periodo medido con `2π/ω₀`.

## Etapa IV: viscosidad y movimiento interno (solo si el contraste lo exige)
Criterio: si el modelo reproduce la frecuencia pero no el decaimiento de la amplitud, falta disipación; si ni siquiera reproduce la frecuencia, el líquido probablemente no se mueve pegado a la lata.

Opciones, de menor a mayor complejidad:
1. Añadir resistencia de rodadura (disipación de la lata, no del líquido) y ver si basta.
2. Añadir un segundo grado de libertad `φ`, el ángulo del líquido respecto a la lata, como péndulo interno acoplado a la pared con amortiguamiento viscoso (candidato natural al modo de sloshing).
3. Estimar los parámetros del amortiguamiento (coeficiente viscoso, brazo efectivo) a partir del decaimiento observado.

La viscosidad entra solo con una justificación física explícita para el tipo de amortiguamiento elegido.

## Otros niveles de llenado
Implementar el segmento circular real (`modelo.md`, sec. 6) y verificarlo contra una referencia de estática.

## Etapa VI final
- Contrastar con varios videos y distintas condiciones iniciales.
- Reportar residuos y métricas por experimento.
- Probar la capacidad predictiva: ajustar con un experimento y predecir otro distinto.

## Preguntas abiertas
1. ¿Superficie horizontal o plano inclinado? (el modelo supone horizontal)
2. ¿Se confirmó con el profesor la lectura de la geometría (medio cilindro de líquido, mismas medidas que la lata)?
3. ¿Cuál es el llenado real de las latas usadas? ¿Realmente 50 %?
4. ¿Qué líquido se usa (agua u otro) y se planea variar la viscosidad?
5. ¿La lata es transparente o hay forma de observar el movimiento del líquido? Permitiría decidir si es solidario, libre o intermedio.
6. ¿Qué fps y qué resolución tiene la cámara?
7. ¿Se medirá `θ(t)` con una marca en la lata o solo `x(t)`?
8. ¿Se quiere variar el nivel de llenado?

## Si se entrega a otra IA
Retomar así: revisar las preguntas abiertas con el estudiante; correr `--selftest`; ayudar a organizar el primer experimento; con el primer CSV real, correr la comparación y analizar residuos; decidir con ellos si se pasa a la Etapa IV. Reglas: las del README, y mantener la geometría de `modelo.md` salvo que el profesor la corrija.
