function Robustness_info = Test_Robust_Structures(PreRegResults,RobustStruct)

% Load the Results of the calibration
load(PreRegResults);
load(RobustStruct);

% EXPERIMENTAL DATA SHOULD BE AVAILABLE IN WORKSPACE

n = length(cmp_group.CC);

% Assing the pre-regression analysis
k_SS = it_results.k_SS;

Robustness_info = cell(n,1);

    for i=1:n
        % Upload the structure to test
        Struct = cmp_group.sols(i,:);
        assignin('base','kfixed',Struct);
    
        % Isolate indexes of unfixed parameters
        index = find(isnan(Struct));
    
        % Parameters to analyze
        k = k_SS(index);
        
        [AICc,CI,CC,Mc,Ms,diff] = reg_analysis_complete(k);
        
        Robustness_info{i}.kfixed     = Struct;
        Robustness_info{i}.k_SS       = k;
        Robustness_info{i}.AICc       = AICc;
        Robustness_info{i}.CI         = CI;
        Robustness_info{i}.CC         = CC;
        Robustness_info{i}.Mc         = Mc;
        Robustness_info{i}.Ms         = Ms;
        Robustness_info{i}.diff       = diff;
        %Robustness_info{i}.ess_report = load('ess_report.mat');
        

    
        close all
    end

end