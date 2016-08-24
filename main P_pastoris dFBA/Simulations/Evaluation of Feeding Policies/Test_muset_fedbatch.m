clear all % Sacar si quiero probar un k_best

%% Parameter values

k_param = [   2.3102; ...      % Vmax
        0.9083;...       % Kg
        0.5255;...     % pEtOH
        0.0703;...    % pPyr
        0.1174;...    % pArab
        1.7082e-5;...    % pCit
        1e-3;...       % pTau
        -1.2;...   % pEtOH FED
        -.4399;...   % pPyr FED
        -.1233;...   % pArab FED
        0;...       % pCit FED
        5e-4;...       % pTau FED
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
odeTime = 0:2:100;
mu_maxs = 0.15; % h^-1
mu_rates = [0.07 0.5];
mu_mins = 0.03;
conditions = length(mu_maxs)*length(mu_rates)*length(mu_mins);
simTimes = cell(conditions,1);
simVariables = cell(conditions,1);
O2_demand = cell(conditions,1);
O2_rates = cell(conditions,1);
Biom = cell(conditions,1);

% Checks that Pichia ate the remaining ethanol, arabitol and pyruvate
ateEtOH = 0;
atePyr  = 0;
ateArab = 0;
t_feed_carb = 0;
assignin('base','ateEtOH',ateEtOH)
assignin('base','atePyr',atePyr)
assignin('base','ateArab',ateArab)

% Oxygen transfer rate limitation
kLa = 150; 
max_OTR = kLa*(7.8*5-2.8);

% Start feed
t_feed = 35;
post_starv = 0;
assignin('base','t_feed',t_feed)
assignin('base','post_starv',post_starv)
l = 1;

for i=1:length(mu_maxs)
    for j=1:length(mu_rates);
        for k=1:length(mu_mins)
        
            delete('metMovieWT.mat')
        
            % Change parameters in the feed function
            assignin('base','mu_max',mu_maxs(i));
            assignin('base','mu_rate',mu_rates(j));
            assignin('base','mu_min',mu_mins(k));
       
            % Solve ODE
            [tWT,xWT] = run_PPdFBA(k_param,x0,odeTime,[],[],[]);
            simTimes{l} = tWT;
            simVariables{l} = xWT;
    
            load('metMovieWT.mat')
    
            O2 = get_O2_for_input_times(fluxDistrib,tWT);
            O2_rates{l} = O2;
            Biom{l} = xWT(:,2);
            O2_demand{l} = abs(O2.*xWT(:,2)*32);
            
            display(['check: ' num2str(l)])
                        
            %% Plots
            figure(l)
            labels = {  'Volume' 'Biomass' 'Glucose' 'Ethanol' 'Pyruvate'...
                        'Arabitol' 'Citrate' 'Thaumatin'};
            for m=1:8
                subplot(3,3,m)
                plot(tWT,xWT(:,m))
                xlabel('[h]')
                ylabel('[g/L]')
                title(labels(m))
            end
            suptitle([  '\mu_M_A_X: ' num2str(mu_maxs(i))...
                        ', \mu_M_I_N: ' num2str(mu_mins(k))...
                        ', rate: ' num2str(mu_rates(j))])

            l = l+1;
        end
    end
end

optTime =toc(Time1);

%% Calculation of performance indicators
Final_Biomass   = zeros(conditions,1);
Final_Thaumatin = zeros(conditions,1);
Volumetric_Prod = zeros(conditions,1);
Biomass_Prod    = zeros(conditions,1);

for i=1:conditions
    Final_Biomass(i)    = simVariables{i,1}(end,2);
    Final_Thaumatin(i)  = simVariables{i,1}(end,8);
    Volumetric_Prod(i)  = Final_Thaumatin(i)/simTimes{i,1}(end);
    Biomass_Prod(i)     = Final_Biomass(i)/simTimes{i,1}(end);
end

save('Final_Biomass.mat','Final_Biomass')
save('Final_Thaumatin.mat','Final_Thaumatin')
save('Volumetric_Prod.mat','Volumetric_Prod')
save('Biomass_Prod.mat','Biomass_Prod')


break
figure(13)
mu_value = cell(1,conditions);
lineshapes = {  'k-' 'b-' 'g-' 'r-' 'c-' 'y-' ...
                'r--' 'b--' 'g--' 'k--' 'c--' 'y--'...
                'r.-' 'b.-' 'g.-' 'k.-' 'c.-' 'y.-'};
for i=1:conditions
    plot(simTimes{i},O2_demand{i},lineshapes{i})
    hold on
    % mu_value{i} = ['\mu ' num2str(mu_sets(i)) ' h^-^1'];
end
    plot([0 odeTime(end)],[max_OTR max_OTR],'r-')
%legend(mu_value,'Location','Best')
ylabel('Oxygen Uptake [mg/Lh]')
xlabel('Time [h]')

save('time mu_set_01_ox.mat','tWT')
save('StateVariables mu_set_01_ox.mat','xWT')