clear all
n_struct = 9;
% struct = [1 2 8]; % Structure 1: Vmax, Ks, alpha
% struct = [1 2 6 8]; % Struture 2: Vmax, Ks, fCit, alpha
%struct = [1 3 4 5]; % Structure 5: Vmax,fet,fpyr,farab
%struct = [1 6 8 9]; % Structure 7: Vmax, fcit, alpha y matp
struct = [1 3 4 6 8]; % Structure 9: Vmax,fetoh, fpyr, fcit, alpha.
kSS = [5 2.7e-3 0 0 0 0 0 0 2.18];
kSS(struct) = NaN;

% Filename and datasheet
filename = 'Calibration data Picha pastoris dFBA.xlsx';
dataset = 9;
time1 = tic;
iteration_complete_valid(filename,dataset,kSS,n_struct);
t1 = toc(time1);
