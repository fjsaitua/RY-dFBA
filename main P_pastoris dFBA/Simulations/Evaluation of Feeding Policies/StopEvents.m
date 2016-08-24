function [value,isterminal,direction] = StopEvents(t,y,p)
% Declare important variables
O2_instant_demand = evalin('base','O2_instant_demand');

kLa = 200; % h-1


% Restrictions
value(1) = y(1)-1; % Volume less than 1L
value(2) = O2_instant_demand-kLa*(5*7.8-2.8); % Oxygen demand more than 2 L/min
isterminal = ones(size(value));
direction = zeros(size(value));

end