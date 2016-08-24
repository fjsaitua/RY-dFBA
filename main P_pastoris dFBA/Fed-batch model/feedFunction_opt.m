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

function F = feedFunction_opt(t,mu_set)

%d = evalin('base','dataset'); %Dataset number
d = 5;

if d == 5 % 9A fermentation
    
    ti      = 30.63;
    mu_set  = mu_set;
    Sin     = 534;
    X0_OD   = 32.8;
    ODDW    = 0.8;
    X0_gL   = X0_OD*ODDW;
    V0      = 0.386; % Ojo a que no se considera la disminución de volumen en la modelación
    Ysx     = 1*ODDW; % gramos biomasa / gramos glucosa
    
    % Parámetros exponencial normal
    a = mu_set*(X0_gL*V0)/(Sin*Ysx);
    b = mu_set;
    
else
    ti = 2000;
    a = 0.001;
    b = 0.0664;
end

if t <= ti
    F = 0;
else
    F = a*exp(b*(t-ti));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%