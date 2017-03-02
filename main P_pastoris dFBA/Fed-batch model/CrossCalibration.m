clear all
n_struct = 33;

% Structure 3: Vmax, Ks, fpyr, farab, fcit, v_pyr_FB, alpha_batch, alpha_FB, m_atp, Tcons;
struct = [1 2 4 5 6 8 11 12 13 14];

% Predefined values to fix parameter not included in the adjustable set
kSS = [ 6;      % Vmax
        2.7e-3; % Ks
        0;      % v_etoh,batch
        0;      % v_pyr,batch
        0;      % v_arab,batch
        0;      % v_cit,batch 
        1.21;   % v_etoh,FB
        0.14;   % v_pyr,FB
        0.15;   % v_arab,FB
        0;      % v_cit,FB
        0;      % alpha batch
        0;      % alpha FB
        2.18;   % m_atp
        22];    % initial time


kSS(struct) = NaN;

% Filename and datasheet
filename = 'Additional File 11 - Fed-batch model calibration data.xlsx';
datasets = 4;
CrossCal_times = zeros(size(datasets));

for i = 1:length(datasets)
    time1 = tic;
    iteration_complete_valid(filename,datasets(i),kSS,n_struct);
    t1 = toc(time1);
    CrossCal_times(i) = t1;
end
