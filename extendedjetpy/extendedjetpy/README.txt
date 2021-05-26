La forma más simple de compilar el programa es utilizando "make" (tiene que estar instalado).

Desde la terminal ejecutamos:
> make

Al compilar genera un ejecutable "run.x". Lo corremos como:
> ./run.x

Para borrar todos los archivos que crea al compilar:
> make clean
(esto no borra los archivos .f90 ni los gráficos)

----------------------------------------------------------------------------

Explicación de los archivos:

* Makefile -> Se fija en la carpeta cuáles son todos los archivos que terminan en .f90 y busca en Makefile.depend cuáles son las dependencias entre ellos. Es el que se encarga de la compilación básicamente.

* Makefile.depend -> Acá se especifica de qué otros módulos depende cada módulo, para que defina el orden de compilación. Si se crean módulos nuevos (para calcular otros procesos radiativos, por ejemplo) es necesario actualizar este archivo.

* global.f90 -> Tiene metidas todas las constantes, por lo que prácticamente todos los módulos lo usan. Si se cambia un parámetro de acá, hay que volver a compilar el código (con make) y correr el progrgama.

* main.f90 -> Programa principal que llama a las distintas subrutinas/módulos.

* controles.f90 -> Se especifican qué módulos vamos a querer correr (puede ser útil cuando se está probando un nuevo módulo y no se quiere correr el programa entero para no perder tiempo). En este caso son sólo 2 módulos y corre rápido, así que no hace falta tocarlo.

* dist.f90 -> Calcula la distribución de partículas y tiempos de enfriamiento.

* plot_diste.sc y plot_tiempose.sc -> scripts de gnuplot, se pueden llamar como:
> gnuplot plot_diste.sc (sin entrar en gnuplot primero)
gnuplot> load plot_diste.sc (desde gnuplot)

Además de graficar las curvas en sí, ajustan algunas rectas haciendo uso de un filtrado con awk. Puede ser útil para hacer algunos chequeos.

*fixbb -> darle permiso de ejecución; lo llaman los scripts de gnuplot para recortarle los bordes blancos a la imagen (es muy útil si luego hay que incorporarlas en un artículo, tesis, diapositiva...).
