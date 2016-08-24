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

function F = feedFunction_PATIL(t)

% Initial glucose concentration of the batch culture
ti      = evalin('base','t_feed');
mu_set  = evalin('base','mu_set');
superUmbral = evalin('base','superUmbral');
t_feed_carb = evalin('base','t_feed_carb');

if t <= ti && superUmbral == 0;
    
    F = 0;
    
elseif t > ti && superUmbral == 0
    
    % Rest of fed-batch parameters
    Sin     = 534;      % [g/L]
    X0_gL   = evalin('base','X0_FB');  % [g/L]
    V0      = evalin('base','V0_FB'); % Ojo a que no se considera la disminución de volumen en la modelación
    Ysx     = 0.6549; % gramos biomasa / gramos glucosa

    % Parámetros exponencial normal
    a = mu_set*(X0_gL*V0)/(Sin*Ysx);
    b = mu_set;
    if t_feed_carb < ti
        F = a*exp(b*(t-ti));
    else 
        F = a*exp(b*(t-t_feed_carb));
    end
    
elseif t > ti && superUmbral == 1
    k       = evalin('base','k_post');
    V0      = evalin('base','V0_FB'); % Ojo a que no se considera la disminución de volumen en la modelación
    t_start = evalin('base','t_start');
    
    F = k*V0*exp(k*(t-t_start));
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%