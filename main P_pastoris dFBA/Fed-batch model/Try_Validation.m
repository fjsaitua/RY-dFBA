clear all
n_struct = 33;

% Structure 1: Vmax, Ks, fpyr, fcit, v_etoh_FB, v_pyr_FB, v_arab_FB, v_cit_FB, alpha_batch, m_ATP, T_cons.
%struct = [1 2 4 6 8 9 10 11 13 15 16]; 

% Structure 3: Vmax, Ks, fpyr, farab, fcit, v_pyr_FB, alpha_batch, alpha_FB, m_atp, Tcons;
struct = [1 2 4 5 6 9 13 14 15 16];

% Structure 4: Vmax, fetoh, fpyr, farab, fcit, v_arab_FB, alpha_batch, alpha_FB, m_atp, Tcons;
%struct = [1 3 4 5 6 10 13 14 15 16];

% Structure 3 + v_etoh,b
%struct = [1 2 3 4 5 6 9 13 14 15 16];

% Prededined values
kSS = [ 6;      % Vmax
        2.7e-3; % Ks
        0;      % v_etoh,batch
        0;      % v_pyr,batch
        0;      % v_arab,batch
        0;      % v_cit,batch 
        0;      % v_tau,batch
        1.21;   % v_etoh,FB
        0.14;   % v_pyr,FB
        0.15;   % v_arab,FB
        0;      % v_cit,FB
        0;      % v_tau,FB
        0;      % alpha batch
        0;      % alpha FB
        2.18;   % m_atp
        22];    % Tiempo inicio consump de metabolitos secundarios


kSS(struct) = NaN;

% Filename and datasheet
filename = 'Datos Calibración Fed-batch.xlsx';
datasets = 4;
CrossCal_times = zeros(size(datasets));

for i = 1:length(datasets)
    time1 = tic;
    iteration_complete_valid(filename,datasets(i),kSS,n_struct);
    t1 = toc(time1);
    CrossCal_times(i) = t1;
end
