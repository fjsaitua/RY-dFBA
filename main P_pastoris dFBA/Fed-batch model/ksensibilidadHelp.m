% KSENSIBILIDAD calcula la sensibilidad de los estados xi con respecto a
% los parámetros del modelo, es decir, calcula numéricamente dxj/dki.
%
% KSENSIBILIDAD requiere los siguientes 'inputs':
%
%                 ksensibilidad(NumEqs,k,X0,tpo,NombreModelo,flag)
%
% NumEqs = las ecuaciones diferenciales que serán estudiadas. Es posible
%   escoger un subconjunto de ellas. Por ejemplo, podría darse
%   que en un modelo con cinco ecuaciones diferenciales, se desea estudiar
%   por ejemplo, la primera, tercera y cuarta, de tal forma que NumEqs = [1 3 4].
%   Por otra parte, si las cinco ecuaciones fueran las estudiadas, entonces
%   NumEqs = [1:5]. El orden de las ecuaciones debe ser ascendente.
%
% k = valor nominal del vector de parámetros a estudiar. En torno a estos
%   valores se realizará el análisis de sensibilidad. Si no son conocidos a
%   priori, pueden estimarse mediante mínimos cuadrados (e.g., LSQCURVEFIT).
%
% X0 = condición inicial del modelo diferencial completo.
%
% tpo = vector de valores tiempo que se utilizarán en el estudio de
%   sensibilidad. Si solo se especifican dos valores, KSENSIBILIDAD asumirá
%   que se ingresó el tiempo inicial y final de integración.
%
% NombreModelo = nombre del archivo *.m que contiene al modelo diferencial.
%   Deberá ingresarse sin extensión y con un símbolo arroba (@) antes del
%   nombre del archivo. El modelo deberá tener como entrada: la variable
%   independiente, la variable dependiente y luego el vector de parámetros.
%   Por ejemplo, modelo(t,x,k).
%
% Flag = 1, entrega la sensibilidad absoluta dxi/dkj.
% Flag = 2, entrega la sensibilidad kj*dxi/dkj.
% Flag = 3, entrega la sensibilidad relativa (kj/xi)*dxi/dkj. En este caso
%   se entregan gráficos resumen, en donde se muestra el efecto de los
%   parámetros sobre las variable de estado en un mismo gráfico.
%
% Los gráficos que no exhiban sensibilidad con respecto a algún parámetro
% no se desplegarán.
%
% Una vez finalizado el análisis, KSENSIBILIDAD guardará dos figuras:
% La primera, Sensib1_HHMMSS, resume todas las sensibilidades en un solo
% gráfico. Si el número de ecuaciones o parámetros es muy grande, esta
% figura no se mostrará. La segunda figura es Sensib2_HHMMSS, la cual muestra
% una comparación relativa del impacto de cada parámetro en las distitas
% ecuaciones de estado. El método para calcular la comparación relativa es
% una modificación de "S" presentado en:
%
%   Hao et al. Modeling the VPAC2-Activated cAMP/PKA Signaling Pathway: From
%   Receptor to Circadian Clock Gene Induction. Biophysical Journal,
%   Vol 90:1560–1571 (2006).
%
% en cambio, KSENSIBILIDAD calcula "S" como (usa valor absoluto de las
% diferencias del numerador):
%
%        S = sum(abs(y(ti)-ybase(ti)))/sum(abs(ybase(ti)))*((delta k)/k)^-1
%
% NOTA 1: HHMMSS se refiere a la hora, minuto y segundo en que fueron
% guardadas las figuras Sensib1 y Sensib2.
%
%--------------------------------
% El siguiente ejemplo analiza el modelo diferencial 'bioreactor':
%
%   NumEqs = [1:3];
%   k      = [20.99 0.00078 1.15 2.6 0.00061 2.95 0.894 0.0913];
%   X0     = [0.004641 3.5/1000 0 0.004641 0 0 0 182.55/1000];
%   tpo    = [0 142];
%   NombreModelo = @bioreactor;
%   flag = 2;
%
%   ksensibilidad(NumEqs,k,X0,tpo,NombreModelo,flag)
%
% En este ejemplo se estudian las tres primeras ecuaciones diferenciales.
% Luego se ingresan los valores nominales de los 8 parámetros, seguidos de
% las condiciones iniciales de integración. Para el tiempo de integración se
% ingresó únicamente el tiempo inicial (t = 0) y final (t = 142). El archivo
% con el modelo diferencial se denomina bioreactor.
%--------------------------------
%
% NOTA 2: KSENSIBILIDAD NO SE EJECUTARÁ SI SU NOMBRE SE MODIFICA.
%
% por <a href = "http://www.iiq.cl" >Claudio Gelmi, Ph.D.</a>
% Departamento de Ingeniería Química y Bioprocesos,
% Pontificia Universidad Católica de Chile.
%
% Última actualización: 19/09/2013.