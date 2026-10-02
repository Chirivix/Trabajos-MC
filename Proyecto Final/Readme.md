# Reto IV-2: Dinámica de una botella parcialmente llena que rueda

Bienvenido a la carpeta del **Reto IV-2** del curso de Mecánica Clásica. Este proyecto aborda la resolución de un problema físico complejo en el que no existe una única respuesta "correcta" predefinida. Nuestro objetivo es formular un modelo consistente, establecer hipótesis claras y contrastar resultados computacionales con datos experimentales, apoyados por asistentes de Inteligencia Artificial bajo una estricta auditoría.

## Equipo de Trabajo y Roles
Siguiendo los lineamientos del curso, el proyecto se desarrolla distribuyendo e intercambiando funciones fundamentales. 
* **David Felipe Chirivi Carreño** - [Líder de Equipo]
* **[Andrés Santiago Santander Fonseca]** - [Coordinador general]

*Roles a rotar:* 
1. **Modelado y análisis:** Formulación física, coordenadas generalizadas, Lagrangiano, aproximaciones y cantidades conservadas.
2. **Cálculo, simulación y datos:** Integración numérica, análisis de datos y control de errores.
3. **Validación y auditoría de IA:** Comprobación independiente, análisis dimensional, límites y revisión crítica de respuestas generadas por IA.

## Descripción del Problema
Cuando una botella cilíndrica parcialmente llena de líquido rueda sobre una superficie horizontal, su velocidad de traslación presenta oscilaciones perceptibles. 

Buscamos determinar qué produce estas oscilaciones, de qué depende su amplitud, cómo interviene el nivel de llenado y qué papel juega la viscosidad. Especialmente, nos enfocamos en encontrar una **magnitud observable discriminante** para separar los efectos del desplazamiento del centro de masa (CoM) respecto al *sloshing* interno del líquido, e investigar escenarios de resonancia.

## Metodología de Resolución
Nuestra solución no depende ciegamente del ajuste de parámetros. El trabajo está respaldado por:
1. **Modelado Teórico:** Planteamiento de un Lagrangiano para el sistema (cilindro rodante + péndulo equivalente para el *sloshing*), definiendo claramente qué efectos incluir y cuáles despreciar.
2. **Diseño Experimental Controlado:** Pruebas con un cilindro rígido aislando variables: 3 líquidos (distintas viscosidades), 3 niveles de llenado y 3 velocidades iniciales.
3. **Computación:** Extracción de datos mediante videometría con **Tracker** y resolución de ecuaciones de movimiento vía integración numérica (RK4) en **Python/Colab**.

## Uso de IA y Protocolo de Validación
El uso de asistentes de IA en este repositorio es activo (exploración de derivaciones, álgebra intermedia, generación de código y sugerencia de aproximaciones). Sin embargo, **ninguna respuesta de la IA constituye una verificación**. 

Toda conclusión importante de este proyecto ha sido validada independientemente mediante los siguientes mecanismos:
* Comparación con los datos de nuestro propio experimento.
* Análisis dimensional de las ecuaciones formuladas.
* Comportamiento en límites conocidos (ej. comportamiento de cuerpo rígido con botella al 0% y 100% de llenado).
* Conservación de la energía mecánica total en las simulaciones (para el caso sin viscosidad).

Todo uso de IA sigue la **metodología 4D** discutida en el curso.

## structura del Repositorio

* `📁 bitacora_IA/` : Registro de interacciones con la Inteligencia Artificial. Incluye los prompts, respuestas originales y la auditoría/análisis crítico donde se demuestra qué partes de la solución de la IA son confiables y cuáles fueron descartadas.
* `📁 data/` : Archivos `.csv` exportados desde Tracker con observaciones crudas y procesadas.
* `📁 src/` : Código fuente, simulaciones computacionales y cuadernos (notebooks) usados para análisis de datos y resolución de ecuaciones.
* `📁 docs/` : Documentos teóricos, deducción analítica de modelos, derivaciones matemáticas en LaTeX y reportes.
* `📁 media/` : Fotografías del montaje experimental, videos de las pruebas e imágenes/gráficos generados.
