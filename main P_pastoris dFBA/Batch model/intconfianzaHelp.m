% INTCONFIANZA calcula intervalos de confianza cuando se ajustan par�metros
% a ecuaciones diferenciales o modelos no lineales usando la funci�n LSQCURVEFIT.
% Esta funci�n no requiere del toolbox de estad�stica para ser ejecutada y
% es un reemplazo a la funci�n NLPARCI disponible en el toolbox de
% estad�stica. Sin embargo, se requiere el Jacabiano ,J, calculado por la
% rutina LSQCURVEFIT, disponible en el toolbox de optimizaci�n.
%
% INTCONFIANZA(J,RESID,PARAMS,ALPHA) entrega intervalos de 100*(1-alpha)%
% de confianza para el ajuste de par�metros de modelos no lineales, sean
% �stos din�micos o no. Los argumentos de entrada a la funci�n INTCONFIANZA
% son:
%
%   J = Jacobiano calculado por la funci�n LSQCURVEFIT del toolbox de
%       optimizaci�n.
%
%   RESID = diferencia entre valores predichos por el modelo y los datos
%       experimentales. Este vector es calculado por la funci�n
%       LSQCURVEFIT.
%
%   PARAMS = valor de los par�metros encontrados por LSQCURVEFIT.
%
%   ALPHA = nivel de significancia, usualmente 0.05 (i.e., IC de 95%).
%
% Las salidas de la funci�n INTCONFIANZA son DELTA y CV. DELTA corresponde
% al intervalo +/- 100*(1-alpha) de confianza de cada par�metro ajustado.
% La segunda salida (opcional) CV corresponde al coeficiente de variaci�n,
% definido aqu� para cada par�metro "i" como 100*sigma_i/param_i (sigma_i = 
% desviaci�n est�ndar del par�metro i). Este valor corresponde al PORCENTAJE
% que representa la incertidumbre con respecto al tama�o del par�metro.
%
% Mediante el uso del Jacobiano, el algoritmo hace una estimaci�n de la
% matriz de covarianza de los par�metros, de tal forma que:
%
%       cov = inv(J'*J)*sum(resid.^2)/(n-nparam)
%
% El valor inverso de la distribuci�n acumulada de Student se program�
% usando la integral de la distribuci�n t, la cual tambi�n fue
% incluida en esta rutina:
%
%       t-student(x,v) =
%       gamma((v+1)/2)/((pi*v)^0.5*gamma(v/2)*(1+x^2/v)^((v+1)/2))
%
% donde v son los grados de libertad de la distribuci�n.
%
% Recomendaciones:
%
% 1. El c�lculo del jacobiano se ve afectado por el tipo de integrador y
% su precisi�n utilizada para integrar el sistema diferencial. Para obtener
% resultados confiables se recomienda utilizar el integrador ODE113 o ODE15s,
% este �ltimo con valores 'RelTol' y 'AbsTol' del orden de 1e-7. Si el error
% de integraci�n no es controlado adecuadamente, los valores del intervalo
% de confianza no ser�n confiables.
%
% 2. Si se desea usar LSQCURVEFIT solo con el fin del c�lculo del
% jacobiano, se recomienda usar optimset('MaxFunEvals',1,'MaxIter',1) y
% definir como punto de partida inicial la soluci�n �ptima encontrada (X*).
% De esta forma, el optimizador solo calcular� el jacobiano en el punto �ptimo,
% sin modificar la soluci�n X*.
%
% Referencias:
%
%   P. Englezos y N. Kalogerakis, Applied Parameter Estimation for Chemical
%   Engineers. Marcel Dekker Inc, New York, 2001.
%
%   K. Krishnamoorthy, Handbook of Statistical Distributions with Applications.
%   Chapman & Hall/CRC, Boca Raton, 2006.
%
% NOTA: INTCONFIANZA NO SE EJECUTAR� SI SU NOMBRE SE MODIFICA.
%
% por Claudio Gelmi, Ph.D. (www.iiq.cl)
% Departamento de Ingenier�a Qu�mica y Bioprocesos,
% Pontificia Universidad Cat�lica de Chile.
%
% �ltima actualizaci�n: 13/01/2012.