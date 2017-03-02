clear all

% Indicate the structure to be assessed
n_struct = 9;

% Indicate the ajustable parameters of the structure
struct = [1 3 4 5 6 8]; % Structure 1: Vmax,fetoh, fpyr, farab, fcit, alpha.
kSS = [5 2.7e-3 0 0 0 0 0 2.18];
kSS(struct) = NaN;

% Filename and datasheet
filename = 'Additional File 10 - Batch model calibration data.xlsx';

% Indicate the dataset to be calibrated with the reduced model structure
dataset = 9;
time1 = tic;
iteration_complete_valid(filename,dataset,kSS,n_struct);
t1 = toc(time1);
