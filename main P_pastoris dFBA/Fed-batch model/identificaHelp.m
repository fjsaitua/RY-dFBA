% IDENTIFICA calcula la matriz de correlación (Mc) con el fin de determinar
% si los parámetros de un modelo son 'localmente identificables' (a priori)
% o no, según el vector tiempo (tpo) dado. Parámetros que son 'localmente
% identificables' poseen correlaciones entre -1 y 1 con todos los otros
% parámetros. Parámetros que no son 'localmente identificables' poseen
% correlaciones iguales a +1 o -1 con, al menos, otro parámetro.
%
% IDENTIFICA requiere los siguientes 'inputs':
%
%                 identifica(NumEqs,k,X0,tpo,NombreModelo,umbral)
%
% NumEqs = ecuaciones diferenciales que cuentan con información experimental
%   o que son observables. En algunos casos, no todos los estados
%   del modelo han sido son medidos u observados. Por ejemplo, podría darse
%   que en un modelo con cinco ecuaciones diferenciales, se desea estudiar
%   un subconjunto de ellas; en ese caso, tendremos que NumEqs = [1 3 4].
%   Por otra parte, si las cinco ecuaciones fueran estudiadas, entonces
%   NumEqs = [1:5]. El orden de las ecuaciones debe ser ascendente.
%
% k = valor nominal del vector de parámetros a estudiar. En torno a estos
%   valores se realizará el análisis de sensibilidad. Si no son conocidos
%   a priori, pueden estimarse mediante mínimos cuadrados (e.g., LSQCURVEFIT).
%
% X0 = condición inicial del modelo diferencial completo.
%
% tpo = vector de valores tiempo que se utilizarán en el estudio de
%   identificabilidad. Por ejemplo, [0:2:142]. Si solo se especifican dos
%   valores, IDENTIFICA asumirá que se ingresó el tiempo inicial y final de
%   integración. En caso que se desee estudiar valores promedios de correlación,
%   pueden utilizarse arreglos aleatorios de tiempo creciente.
%
% NombreModelo = nombre del archivo *.m que contiene al modelo diferencial.
%   Deberá ingresarse sin extensión y con un signo arroba (@) antes del nombre
%   del archivo. El modelo deberá tener como entrada: la variable
%   independiente, la variable dependiente y el vector de parámetros.
%   Por ejemplo, modelo(t,x,k).
%
% Umbral = valor mínimo, a partir del cual es aceptable la correlación
%   entre parámetros. Si IDENTIFICA encuentra parámetros que poseen correlaciones
%   mayores a UMBRAL, los mostrará al final del análisis y avisará que podrían
%   existir problemas de identificabilidad. Por defecto, puede usar un valor
%   de 0.95. Si IDENTIFICA encuentra valores mayores a 0.99, los indicará
%   automáticamente; en ese caso existen serios problemas de identificabilidad.
%
% La matriz de correlación, Mc, es guardada automáticamente en el archivo
% Mc.txt y es graficada usanda un mapa de colores para indicar los distintos
% valores de correlación.
%
% ----------------------------------------
% El siguiente ejemplo analiza el modelo diferencial 'bioreactor':
%
%   NumEqs = [1:8];
%   k      = [20.99 0.00078 1.15 2.6 0.00061 2.95 0.894 0.0913];
%   X0     = [0.004641 3.5/1000 0 0.004641 0 0 0 182.55/1000];
%   tpo    = [0 142];
%   NombreModelo = @bioreactor;
%   umbral = 0.95;
%
%   identifica(NumEqs,k,X0,tpo,NombreModelo,umbral)
%
% En el ejemplo se estudian 8 ecuaciones diferenciales, 8 parámetros, seguidos
% de las condiciones iniciales de integración. Para el tiempo de integración se
% ingresó únicamente el tiempo inicial (t = 0) y final (t = 142). El archivo
% con el modelo diferencial se denomina bioreactor. El umbral utilizado fue 0.95.
% ----------------------------------------
% Referencias:
%
% 1. Jacquez J.A. and Greif P. (1985). Numerical parameter identifiability and
% estimability: Integrating identifiability, estimability, and optimal
% sampling design. Math. Biosci. 77:201-227.
%
% 2. Jacquez J.A. and Perry T. (1990). Parameter estimation: local
% identifiability of parameters. Am. J. Physiol. 258:727-736.
%
% NOTA: IDENTIFICA NO SE EJECUTARÁ SI SU NOMBRE SE MODIFICA.
%
% por Claudio Gelmi, Ph.D. (www.iiq.cl)
% Departamento de Ingeniería Química y Bioprocesos,
% Pontificia Universidad Católica de Chile.
%
% Última actualización: 13/01/2012.