# Clasificación de formas parecidas
## El armario
Para organizar tu armario necesitas clasificar todos los tiliches, y quién sabe qué cosas pueda haber ahí. Los zapatos deben ir en la caja de zapatos, las camisas en la caja de camisas y las GTX en la caja de tarjetas gráficas GTX.

La pregunta entonces es, ¿podemos automatizar este proceso? Replanteando, ¿habrá un algoritmo que nos ayude a enseñarle a la computadora a clasificar objetos misceláneos? Hay que tomar en cuenta que ese armario puede tener cualquier cosa, por lo que usar una base de datos puede no ser una buena opción.

Así, podemos basar nuestro criterio en una sola característica: la forma de las cosas. De todos modos, ¿cuál es la probabilidad de que un zapato tome la forma de una tarjeta gráfica?

## A considerar
Ya de una vez, y para aclarar cualquier confusión, es necesario decir que el algoritmo principal que construiremos solamente nos ayudará a cuantificar la similitud entre 2 formas. Si le haces a las mates, vamos a construir una norma, una distancia, en el espacio de contornos, la cual cuantificará la similitud entre 2 figuras (contornos).

Pero bueno, una posible norma sería una que cuente la cantidad de pixeles negros y ya. Aunque esa norma sólo nos va a medir qué tan negra es la imágen, lo cual puede no ser muy útil. Algo que sí queremos que cumpla nuestra distancia es que si estamos comparando un coche con sí mismo, pero movido 10 cm a la izquiera, la distancia sea 0. Del mismo modo si el coche está rotado o, idealmente, si se encuentra en un plano desfasado.
