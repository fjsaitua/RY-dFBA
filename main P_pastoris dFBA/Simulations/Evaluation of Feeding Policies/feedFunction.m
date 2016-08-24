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

t_feed_carb = evalin('base','t_feed_carb');

%d = evalin('base','dataset'); %Dataset number
d = 5;

if d == 5 % 9A fermentation
    
    ti      = evalin('base','t_feed');

    % Decreasing growth rate
    A = evalin('base','mu_min');
    B = evalin('base','mu_max')-A;
    C = evalin('base','mu_rate');

    if t>=ti
        mu_set = A+B*exp(-C*(t-ti));
    else
        mu_set = 0;
    end
    
    assignin('base','instant_mu_set',mu_set)
   
else
    ti = 2000;
    a = 0.001;
    b = 0.0664;
end

if t <= ti
    F = 0;
else
    % Rest of fed-batch parameters
    Sin     = 534;      % [g/L]
    X0_gL   = evalin('base','X0_FB');  % [g/L]
    V0      = evalin('base','V0_FB'); % Ojo a que no se considera la disminución de volumen en la modelación
    Ysx     = 0.6549    ; % gramos biomasa / gramos glucosa
    
    % Parámetros exponencial normal
    a = mu_set*(X0_gL*V0)/(Sin*Ysx);
    b = mu_set;
    if t_feed_carb < ti
        F = a*exp(b*(t-ti));
    else 
        F = a*exp(b*(t-t_feed_carb));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%