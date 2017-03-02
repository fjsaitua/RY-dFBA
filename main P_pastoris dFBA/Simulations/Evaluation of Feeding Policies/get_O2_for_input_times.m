function O2 = get_O2_for_input_times(fluxDistrib,Exp_Times)
% Using the flux distribution during the culture and the sampling time, 
% this script determines the instantaneous oxygen uptake rate for the 
% different experimental times

t_O2 = fluxDistrib(end,:);
r_O2 = fluxDistrib(1318,:);

O2 = zeros(length(Exp_Times),1);
O2(1) = r_O2(1);
O2(end) = r_O2(end);


for i=2:length(Exp_Times)-1
    
    % lower bound
    LB = find(t_O2<Exp_Times(i));
    LB = LB(end);
    
    % upper bound
    UB = find(t_O2>Exp_Times(i));
    UB = UB(1);
    
    % Interpolation
    m = (r_O2(UB)-r_O2(LB))/(t_O2(UB)-t_O2(LB));
    
    O2(i) = r_O2(LB) + m*(Exp_Times(i)-t_O2(LB));
end    


end