% INTCONFIANZA calcula intervalos de confianza cuando se ajustan parámetros
% a ecuaciones diferenciales o modelos no lineales usando la función LSQCURVEFIT.
% Esta función no requiere del toolbox de estadística para ser ejecutada y
% es un reemplazo a la función NLPARCI disponible en el toolbox de
% estadística. Sin embargo, se requiere el Jacabiano ,J, calculado por la
% rutina LSQCURVEFIT, disponible en el toolbox de optimización.
%
% INTCONFIANZA(J,RESID,PARAMS,ALPHA) entrega intervalos de 100*(1-alpha)%
% de confianza para el ajuste de parámetros de modelos no lineales, sean
% éstos dinámicos o no. Los argumentos de entrada a la función INTCONFIANZA
% son:
%
%   J = Jacobiano calculado por la función LSQCURVEFIT del toolbox de
%       optimización.
%
%   RESID = diferencia entre valores predichos por el modelo y los datos
%       experimentales. Este vector es calculado por la función
%       LSQCURVEFIT.
%
%   PARAMS = valor de los parámetros encontrados por LSQCURVEFIT.
%
%   ALPHA = nivel de significancia, usualmente 0.05 (i.e., IC de 95%).
%
% Las salidas de la función INTCONFIANZA son DELTA y CV. DELTA corresponde
% al intervalo +/- 100*(1-alpha) de confianza de cada parámetro ajustado.
% La segunda salida (opcional) CV corresponde al coeficiente de variación,
% definido aquí para cada parámetro "i" como 100*sigma_i/param_i (sigma_i = 
% desviación estándar del parámetro i). Este valor corresponde al PORCENTAJE
% que representa la incertidumbre con respecto al tamaño del parámetro.
%
% Mediante el uso del Jacobiano, el algoritmo hace una estimación de la
% matriz de covarianza de los parámetros, de tal forma que:
%
%       cov = inv(J'*J)*sum(resid.^2)/(n-nparam)
%
% El valor inverso de la distribución acumulada de Student se programó
% usando la integral de la distribución t, la cual también fue
% incluida en esta rutina:
%
%       t-student(x,v) =
%       gamma((v+1)/2)/((pi*v)^0.5*gamma(v/2)*(1+x^2/v)^((v+1)/2))
%
% donde v son los grados de libertad de la distribución.
%
% Recomendaciones:
%
% 1. El cálculo del jacobiano se ve afectado por el tipo de integrador y
% su precisión utilizada para integrar el sistema diferencial. Para obtener
% resultados confiables se recomienda utilizar el integrador ODE113 o ODE15s,
% este último con valores 'RelTol' y 'AbsTol' del orden de 1e-7. Si el error
% de integración no es controlado adecuadamente, los valores del intervalo
% de confianza no serán confiables.
%
% 2. Si se desea usar LSQCURVEFIT solo con el fin del cálculo del
% jacobiano, se recomienda usar optimset('MaxFunEvals',1,'MaxIter',1) y
% definir como punto de partida inicial la solución óptima encontrada (X*).
% De esta forma, el optimizador solo calculará el jacobiano en el punto óptimo,
% sin modificar la solución X*.
%
% Referencias:
%
%   P. Englezos y N. Kalogerakis, Applied Parameter Estimation for Chemical
%   Engineers. Marcel Dekker Inc, New York, 2001.
%
%   K. Krishnamoorthy, Handbook of Statistical Distributions with Applications.
%   Chapman & Hall/CRC, Boca Raton, 2006.
%
% NOTA: INTCONFIANZA NO SE EJECUTARÁ SI SU NOMBRE SE MODIFICA.
%
% por Claudio Gelmi, Ph.D. (www.iiq.cl)
% Departamento de Ingeniería Química y Bioprocesos,
% Pontificia Universidad Católica de Chile.
%
% Última actualización: 13/01/2012.