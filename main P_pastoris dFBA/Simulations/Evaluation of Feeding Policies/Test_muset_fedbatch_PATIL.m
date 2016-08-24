clear all % Sacar si quiero probar un k_best

%% Parameter values

k_param = [   2.3102; ...      % Vmax
        0.9083;...       % Kg
        0.5255;...     % pEtOH
        0.0703;...    % pPyr
        0.1174;...    % pArab
        1.7082e-5;...    % pCit
        1e-7;...       % pTau
        -.8056;...   % pEtOH FED -.8056
        -.4399;...   % pPyr FED -.4399
        -.1233;...   % pArab FED -.1233
        0;...       % pCit FED
        1e-7;...       % pTau FED
        1.1604e-5;...         % w F.O. Batch
        2.0217e-6;...         % w F.O. FedBatch
        2];             % m_atp

%% Initial conditions
x0 = [  0.5; ...    % Volume
        0.1; ...    % Biomass
        50;  ...    % Glucose
        0.1; ...    % Ethanol
        0.001; ...  % Pyruvate
        0.1; ...    % Arabitol
        0.0001; ... % Citrate
        0];         % Thaumatine


%% Simulation
Time1 = tic;
odeTime = 0:2:120;

% Oxygen Data
kLa = 150; % [h-1]
Csat = 7.8*5; % mg/L
CsetPoint = 2.8; % mg/L
mu_set = 0.1;
superUmbral = 0;

% Start feed
t_feed = 35;
t_feed_carb = 0;
post_starv = 0;
assignin('base','t_feed',t_feed)
assignin('base','post_starv',post_starv)
assignin('base','kLa',kLa)
assignin('base','Csat',Csat)
assignin('base','CsetPoint',CsetPoint)
assignin('base','mu_set',mu_set);

% Checks that Pichia ate the remaining ethanol, arabitol and pyruvate
ateEtOH = 0;
atePyr  = 0;
ateArab = 0;  
assignin('base','ateEtOH',ateEtOH)
assignin('base','atePyr',atePyr)
assignin('base','ateArab',ateArab)
        
[tWT,xWT] = run_PPdFBA_PATIL(k_param,x0,odeTime,[],[],[]);
figure(1)
subplot(2,3,1)
plot(tWT,xWT(:,1),'b-')
title('Volume')
subplot(2,3,2)
plot(tWT,xWT(:,2),'r-')
title('Biomass')
subplot(2,3,3)
plot(tWT,xWT(:,3),'k-')
title('Glucose')
subplot(2,3,4)
plot(tWT,xWT(:,4),'c-')
title('Ethanol')
subplot(2,3,5)
plot(tWT,xWT(:,5),'g-')
title('Pyruvate')
subplot(2,3,6)
plot(tWT,xWT(:,6),'m-')
title('Arabitol')

% Plot oxygen demand
load('metMovieWT.mat')
O2 = get_O2_for_input_times(fluxDistrib,tWT);
O2_demand = abs(O2).*xWT(:,2)*32; % mg/Lh

figure(2)
plot(tWT,O2_demand,'b-')
title('Oxygen demand')
ylabel('[mg/Lh]')
xlabel('Time [h]')
optTime =toc(Time1);

save('time mu_set_01_ox _final_2.mat','tWT')
save('StateVariables mu_set_01_ox_final_2.mat','xWT')