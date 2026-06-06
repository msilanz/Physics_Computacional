En términos de código, en esta segunda entrega se hacen las siguientes modificiones respecto al código de la primera:
- Se generan dos códigos, en vez de uno.
- Cada código genera múltiples simulaciones al compilarlas y luego correrlas, cada una de estas asociadas a un número de reynolds diferente. Esto se hace
usando una clásica lógica de "maestro y esclavo" para tener un ciclo continuo que vaya realizando las animaciones. Además, cada simulación que se obtiene
como resultado se nombra en función de los FPS que posee y del número de reynolds que modela.
- Se hacen modificaciones técnicas a los vídeos que se obtienen como outputs de los códigos, esencialmente que ahora tienen más definición y que el dominio
alteró un poco para que fuese más fácil de poder observar y analizar (ya no es tan delgado como antes).
- El código para el modelo laminar tiene un flujo laminar incidiendo sobre el objeto fijo que se utiliza (un cilindro), el cuál es el mismo que el del
código de la entrega anterior.
- Por otro lado, el código para el modelo turbulento es igual al laminar solo que ahora se añade como condición de borde unos torbellinos que provocan que
el flujo incidente sea turbulento.

En términos del paper, en esta segunda entrega se hacen las siguientes modificiones respecto al pdf de la primera:
- Ahora sí se establece una estructura de paper para todo el documento, con cada parte claramente identificada y todos los componentes bien claros.
- Se establecen los contenidos que debe tener cada parte del paper, colocándolos explícitamente en cursiva
- Se fijan los objetivos a conseguir con el proyecto, los cuáles se mencionan en la introducción.
- Se determinan las subdivisiones de las secciones de la metodología y los resultados.

La mayor parte del avance realizado en esta fase consistió en desarrollar los códigos entregados, para poder generar todas las simulaciones que se encuentran
adjuntas en la carpeta de esta entrega. Las cuáles serán analizadas en detalle en la entrega final, junto con añadir varias más. Se cree que ahora está mucho
más claro "el norte" del proyecto.
