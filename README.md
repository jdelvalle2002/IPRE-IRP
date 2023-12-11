# IPRE-IRP

## Supuestos generales
Tabajaremos con un solo vehículo, con capacidad de carga ilimitada y que cada vez que se visita un cliente se llena completamente su bodega. Se determina que debe ser visitado un cliente cuando su inventario cae bajo cierto valor umbral en función de la distribución de la demanda. 

## Comentarios sobre Resultados

En primer lugar concretamos la implementación de la política reactiva. En este caso, para cada periodo se obtenía, en función de una predicción de demanda basada en los datos históricos, los clientes que necesitaban ser visitados. A continuación, se calculaba una ruta inicial con Nearest Neighbour y se aplicaba la optimización de la ruta con 2-opt. Finalmente esta ruta se almacenaba y se pasaba al siguiente periodo.

Del análisis de los resultados de la implementación se observó que el comportamiento para cada local, representado en la frecuencia de visitas, estaba fuertemente correlacionado con la proporción de demanda inicial respecto a la capacidad de la bodega. Esto es, si la demanda inicial era muy alta, el local era visitado con mayor frecuencia.

Se implementó una visualización que permite apreciar la periodicidad (o no) de las visitas a cada local. Se puede apreciar que cuando los valores de frecuencia de visitas sobrepasan 0.5 pero aún están lejos de 1 (es decir rondan 0.65) tienen un comportamiento más bien caótico a lo largo del tiempo. En cambio cuando las frecuencias son menores o iguales a 0.5, se observa una periodicidad más clara en las visitas, en especial cuando la frecuencia es cercana a 0.5. Por último, cuando la frecuencia es cercana a 1, se observa que el comportamiento es más bien estable, ya que se visita al local en prácticamente todos los periodos.

Estas observaciones son válidas esencialmente cuando se aleatorizan las demandas inciales. Si es que se trabaja directametne con las instancias de Coelho, se observa que:

1) La correlación entre la frecuencia de visitas y la proporción de demanda inicial respecto a la capacidad de la bodega se mantiene alta, inlcuso más cercana a 1.

2) Se separan notoriamente 2 grupos: los cuya frecuencia ronda 0.4 y los cuya frecuencia ronda 0.6. En el primer caso, se trata de los locales que inicialmente tienen una demanda igual a un tercio de la capacidad de la bodega. En el segundo caso, se trata de los locales que inicialmente tienen una demanda igual a la mitad de la capacidad de la bodega.

3) Si bien hay grupos de locales con frecuencias cercanas, no se observa una periodicidad más uniforme entre las visitas en la visualización, se podría interpretar como que deben compatibilizarse numerosos locales con frecuencias similares más no tan exactas. Sin embargo, la diferencia entre grupos de frecuencias se puede apreciar en la visualización.

Se trabajó con 2 instancias N=10, con una instancia N=25 y con una instancia N=30.
