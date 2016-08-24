%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ktofix = decision(Mc,Ms)
% Constructs ktofix; a vector indicating which parameters should be fixed
% according to sensitivity and identifiability analysis. Returns a vector
% with the parameters to be fixed.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ktofix = decision(Mc,Ms)

kfixed = evalin('base','kfixed');
m      = length(kfixed);
[~,n]  = size(Ms);
Mc2    = zeros(m,m);
Ms2    = ones(m,n);

%Construct kf (row vector with 1 if parameter is fixed and 0 otherwise).
kf = zeros(1,m);
for i = 1:m
    if ~isnan(kfixed(i))
        kf(i)=1;
    end    
end

%Modify correlation matrix (Mc) and sensitivity matrix (Ms) to include
%zeros and ones in the fixed parameters position, respectively.
for i = 1:m
    if kf(i) == 0
        for j = 1:m
            if kf(j) == 0
                Mc2(i,j) = Mc(i-sum(kf(1:i)),j-sum(kf(1:j)));
            end
        end
        Ms2(i,:) = Ms(i-sum(kf(1:i)),:);
    end    
end

%Construct ksensit and kidentif, which are vectors with 1 if the
%corresponding parameter has the corresponding problem, and 0 otherwise.
ksensit  = zeros(1,m);
kidentif = zeros(1,m);
for i = 1:m
    %Sensitivity:
    sensitive = false;
    for j = 1:n
        if Ms2(i,j) >= 0.01;
            sensitive = true;
        end
    end
    if ~sensitive
        ksensit(i) = 1;
    end
    
    %Identifiability:
    for j = 1:m
        if abs(Mc2(i,j)) >= 0.95 && i ~= j
            kidentif(i) = 1;
            kidentif(j) = 1;
        end
    end
end

%ksum is the sum of all problems for each parameter:
ksum = ksensit + kidentif;

%Detect if there is at least one parameter with both problems
%(problem2 = true) or with 1 problem (problem1 = true). The algorithm
%prioritizes parameters that have 2 problems, and after parameters with
%only 1 problem.
problem1 = false;
problem2 = false;
if sum(ismember(ksum,2))>=1
    problem2 = true;
elseif sum(ismember(ksum,1))>=1
    problem1 = true;
end

%Construct ktofix (1 if the analized parameter should be fixed and 0
%otherwise). Assign 1s to ktofix only if the number of problems of the 
%corresponding parameter is of the corresponding priority. 
ktofix = zeros(1,m);
for i = 1:m
    if problem2 && ksum(i) == 2
        ktofix(i) = 1;
    elseif problem1 && ksum(i) == 1
        ktofix(i) = 1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%