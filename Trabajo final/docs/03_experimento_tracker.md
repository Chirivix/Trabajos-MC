# Protocolo experimental y mediciones en Tracker

Se mide con el cuaderno final (`notebooks/Mecanica_clasica2.ipynb`, Nivel 4). El cuaderno sencillo (`notebooks/lata_tracker_colab_sencillo.ipynb`, Nivel 2, líquido rígido) usa el mismo tipo de archivo y sirve como comparación. Régimen de los videos: la lata **rueda de forma continua** y un punto marcado en el borde describe una cicloide deformada por el líquido.

## 1. Variables

| Variable | Cómo se obtiene | Tipo |
|---|---|---|
| $R$ | Calibre | Medida directa |
| $m_c$, $m_l$ | Balanza (lata vacía y llena; $m_l$ = llena − vacía) | Medida directa |
| Fracción de llenado $f$ | Volumen de líquido sobre volumen interior de la lata (el cuaderno acepta 5 a 95 %) | Medida |
| Líquido y $\nu$ | Menú del cuaderno (agua, aceite vegetal, aceite de motor, glicerina) o valor propio | Tabulada o del usuario |
| FPS reales | Ajustes de la cámara (cámara lenta) | Dato del equipo |
| Fin de la cámara lenta | Instante del video donde termina el tramo en cámara lenta (`tiempo_fin_camara_lenta_s`) | Medida del video |
| Escala del video | Una longitud conocida visible en el cuadro (p. ej. el diámetro) | Calibración |
| $x(t)$ y $y(t)$ del **punto marcado en el borde** | Tracker | Medida directa |
| $\dot\theta_0$, $\varphi_0$, $\dot\varphi_0$, $c_{visc}$, $c_{rod}$ | Primer 45 % de los datos | **Ajustadas** |
| $\alpha_0$ | Forma geométrica de la trayectoria $(x, y)$ | Estimada |

Las masas, el radio y el llenado nunca se ajustan.

## 2. Primer experimento

1. Pegar una **marca visible en el borde de la tapa**, de color que contraste con la lata. Pesar la lata vacía y llena y medir $R$.
2. Dar un empujón **firme**: el movimiento debe ser rodadura continua (la lata da varias vueltas), no un vaivén.
3. Grabar de lado, en cámara lenta, con la cámara fija y perpendicular al movimiento. Anotar los FPS reales y el instante en que termina la cámara lenta.
4. En Tracker, calibrar la escala y marcar el punto del borde fotograma por fotograma. **Exportar tiempo, $x$ e $y$.** Con $y$ el cuaderno estima el ángulo inicial de la marca; sin $y$ supone que empieza hacia delante, en la horizontal ($\alpha_0 = \pi/2$).
5. Guardar el CSV en `data/raw/` sin modificarlo y escribir al lado un archivo de parámetros (masas, $R$, llenado, líquido, FPS, incertidumbres, fecha).
6. Abrir el cuaderno, escribir los datos en la celda 1, ejecutar y subir el archivo en la celda 2.

## 3. Cómo leer el resultado

El cuaderno calibra con el primer 45 % y mide el error en el 55 % restante.

- **Mira el RMSE de predicción, no solo el $R^2$.** Una cicloide ideal ya tiene $R^2$ alto cuando la lata rueda rápido.
- **Compara con los niveles anteriores** sobre el mismo video: el cuaderno sencillo (Nivel 2, líquido rígido) y la cicloide ideal (Nivel 1). El Nivel 4 solo se justifica si predice mejor que ellos.
- **Revisa si algún parámetro quedó en el límite de su cota.** El cuaderno acota $c_{visc}$ entre 3·10⁻⁷ y 10⁻² N·m·s y $c_{rod}$ entre 3·10⁻⁹ y 3·10⁻⁵ N·m·s. Un parámetro pegado a una cota indica que los datos no lo determinan (ver `04_hoja_de_ruta.md`).
- **Compara $c_{visc}$ ajustado con el inicial** que imprime el cuaderno (a partir de $\nu$); ese valor inicial no es una predicción confiable.
- Si el error es grande desde el inicio de la predicción: revisar FPS, escala de Tracker, unidad (cm, mm, m) y el punto seguido.

## 4. Cuidados con los datos

- **FPS y cámara lenta.** El cuaderno estima los FPS de reproducción del archivo y corrige el tiempo con los FPS reales. Si el seguimiento saltó fotogramas, esa corrección no vale. Los puntos posteriores al fin de la cámara lenta se descartan.
- **Unidad.** El cuaderno lee la unidad del encabezado y, si no la encuentra, usa la elegida en la celda 1.
- **Sentido de $x$.** Se toma del signo de $x_{final} - x_{inicial}$. Con rodadura continua es fiable; en un vaivén puede equivocarse.
- **Cantidad de datos.** Hace falta al menos 12 puntos válidos; el tramo de calibración usa como máximo 250.
- **Duración de la calibración.** Con un video corto (0,3 s de calibración) los coeficientes de amortiguamiento casi no afectan a los datos de calibración, y el ajuste no los determina. Conviene más tiempo de video o varios videos.
- **Incertidumbres.** Registrar la de cada medida (calibre, balanza, escala del video).
- **No modificar nunca los datos crudos.**
