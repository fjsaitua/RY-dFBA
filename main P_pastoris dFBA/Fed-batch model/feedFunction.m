%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F = feedFunction(t)
% Returns feed of fedbatch fermentation (or 0 if fermentation is batch)
%
% INPUT:
% t             Time of simulation [h]
% 
% OUTPUT:
% F             feed [L/h]
%
% All rights reserved to Benjamín J. Sánchez
%
% Last Update: 2015-02-06 - Franciso Saitua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = feedFunction(t)

    
    ti      = evalin('base','t_feed');
    A = 0.07;
    B = 0.03;
    C = 0.07;
    mu_set = A+B*exp(-C*(t-ti));

if t < ti
    F = 0;
else
    % Rest of fed-batch parameters
    Sin     = 500;      % [g/L]
    X0_gL   = evalin('base','X0_FB');  % [g/L]
    V0      = evalin('base','V0_FB'); % Ojo a que no se considera la disminución de volumen en la modelación
    Ysx     = 0.75913; % gramos biomasa / gramos glucosa
    
    % Parámetros exponencial normal
    a = mu_set*(X0_gL*V0)/(Sin*Ysx);
    b = mu_set;
    F = a*exp(b*(t-ti));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%