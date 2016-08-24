%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dx = pseudoSteadyState(t,x,k)
% Optimizes biomass using FBA, under the pseudo-steady state assumption
%
% INPUTS:
% t             Time of simulation [h]
% x             Concentrations of each variable in the last iteration 
%               ([L] or [g/L])
% k             Kinetic parameters for glucose consumption
% 
% OUTPUT:
% dx            Derivatives of concentrations [g/Lh]
%
% Benjamín J. Sánchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = pseudoSteadyState_PATIL(t,x,k) % cambié k por p ya que no estoy fijando nada 

model       = evalin('base','model');
excMet      = evalin('base','excMet');
excRxn      = evalin('base','excRxn');
PM          = evalin('base','PM');
feed        = evalin('base','feed');
kfixed      = evalin('base','kfixed');
post_starv  = evalin('base','post_starv');
t_feed      = evalin('base','t_feed');
kLa         = evalin('base','kLa');
Csat        = evalin('base','Csat');
CsetPoint   = evalin('base','CsetPoint');
superUmbral = evalin('base','superUmbral');
ateEtOH     = evalin('base','ateEtOH');
atePyr      = evalin('base','atePyr');
ateArab     = evalin('base','ateArab');

%Construct p from the fixed values of kfixed and the estimated values of k.
p = zeros(length(kfixed),1);
j = 1;
for i = 1:length(kfixed)
    if isnan(kfixed(i))
        p(i) = k(j); % ¿?
        j    = j+1;
    else
        p(i) = kfixed(i);
    end    
end


% Fixed modifications of the model
model = fixedModifications(model,p);

%Kinetic constraints
model = kineticConstraints_PATIL(model,t,x,excRxn,p);

%FBA
FBAsol = solveFBA_PATIL(model,t,excRxn,p);

if numel(FBAsol.x) == 0
    FBAsol.x = zeros(length(model.rxns),1);
end

% Save the data from the end of the batch
if t >= (t_feed-1) && post_starv == 0
    assignin('base','X0_FB',x(2))
    assignin('base','V0_FB',x(1))
    post_starv = 1;
    assignin('base','post_starv',post_starv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oxygen Supply Demand Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
O2 = FBAsol.x(findRxnIDs(model,'EX_o2(e)'));
O_supply = kLa*(Csat-CsetPoint); % mg/Lh;
O2_demand = abs(O2*x(2)*32); % mg/Lh


%Integration:
N  = length(excMet);
dx = zeros(N,1);

if ateEtOH + atePyr + ateArab < 3
    
    if O2_demand < O_supply
        %Volume:
        Fin   = feedFunction_PATIL(t);
        dx(1) = Fin;               %[L/h]

        % Biomass:
        mu    = FBAsol.x(excRxn(2,1));  %[1/h]
        dx(2) = mu*x(2)-x(2)/x(1)*Fin;  %[gDW/Lh]

        % Glucose
        i     = excMet(3);
        j     = excRxn(3,1);
        v     = -FBAsol.x(j)*model.S(i,j);              %[mmol/gDWh]
        dx(3) = (feed(3)-x(3))*Fin/x(1) + v*PM(3)*x(2); %[g/Lh]
    
    else
        if superUmbral == 0
            mu      = FBAsol.x(excRxn(2,1));  % [1/h]
            rs      = FBAsol.x(excRxn(3,1))*PM(3);  % [g/g*h]
            Ysx     = abs(mu/rs);
            q_max   = mu*x(2);
            Sf      = 534; % g/L
            k_post       = q_max/(Ysx*Sf);
            t_start = t;
            assignin('base','X0_FB',x(2))
            assignin('base','V0_FB',x(1))
            assignin('base','q_max',q_max);
            assignin('base','mu_set',k_post)
            superUmbral = 1;
            assignin('base','superUmbral',superUmbral);
            assignin('base','k_post',k_post);
            assignin('base','t_start',t_start);
        end

        q_max = evalin('base','q_max');
        k_post = evalin('base','k_post');

        % Volume:
        Fin   = feedFunction_PATIL(t);
        dx(1) = Fin;               %[L/h]

        % Biomass:
        mu    = FBAsol.x(excRxn(2,1));  %[1/h]
%        dx(2) = q_max-k_post*x(2);  %[gDW/Lh]
        dx(2) = mu*x(2)-k_post*x(2);  %[gDW/Lh]

        % Glucose
        i     = excMet(3);
        j     = excRxn(3,1);
        v     = -FBAsol.x(j)*model.S(i,j);              %[mmol/gDWh]
        dx(3) = (feed(3)-x(3))*k_post + v*PM(3)*x(2); %[g/Lh]
    end

else
    
    if O2_demand < O_supply
        %Volume:
        Fin   = feedFunction_PATIL(t);
        dx(1) = Fin;               %[L/h]

        % Biomass:
        mu    = FBAsol.x(excRxn(2,1));  %[1/h]
        dx(2) = mu*x(2)-x(2)/x(1)*Fin;  %[gDW/Lh]

        % Glucose
        i     = excMet(3);
        j     = excRxn(3,1);
        v     = -FBAsol.x(j)*model.S(i,j);              %[mmol/gDWh]
        dx(3) = (feed(3)-x(3))*Fin/x(1) + v*PM(3)*x(2); %[g/Lh]

    else
        
        if superUmbral == 0
            mu      = FBAsol.x(excRxn(2,1));  % [1/h]
            rs      = FBAsol.x(excRxn(3,1))*PM(3);  % [g/g*h]
            Ysx     = abs(mu/rs);
            q_max   = mu*x(2);
            Sf      = 534; % g/L
            k_post       = q_max/(Ysx*Sf);
            t_start = t;
            assignin('base','X0_FB',x(2))
            assignin('base','V0_FB',x(1))
            assignin('base','q_max',q_max);
            assignin('base','mu_set',k_post)
            superUmbral = 1;
            assignin('base','superUmbral',superUmbral);
            assignin('base','k_post',k_post);
            assignin('base','t_start',t_start);
        end

        q_max = evalin('base','q_max');
        k_post = evalin('base','k_post');

        % Volume:
        Fin   = feedFunction_PATIL(t);
        dx(1) = Fin;               %[L/h]

        % Biomass:
        mu    = FBAsol.x(excRxn(2,1));  %[1/h]
%        dx(2) = q_max-k_post*x(2);  %[gDW/Lh]
        dx(2) = mu*x(2)-k_post*x(2);  %[gDW/Lh]

        % Glucose
        i     = excMet(3);
        j     = excRxn(3,1);
        v     = -FBAsol.x(j)*model.S(i,j);              %[mmol/gDWh]
        dx(3) = (feed(3)-x(3))*k_post + v*PM(3)*x(2); %[g/Lh]

    end
end

%Extracellular Metabolites:
for q = 4:N
    i     = excMet(q);
    j     = excRxn(q,1);
    v     = -FBAsol.x(j)*model.S(i,j);              %[mmol/gDWh]
    dx(q) = (feed(q)-x(q))*Fin/x(1) + v*PM(q)*x(2); %[g/Lh]
end


skip_delays = evalin('base','skip_delays');
persistent odetime
maxtime = 300;

if skip_delays
    if isempty(odetime)
        odetime = tic;
    elseif toc(odetime) > maxtime
        odetime = [];
        error('Slow integration');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%