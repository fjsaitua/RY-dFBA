%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CC = reg_analysis_onlyCC(fixed_values,kfixed)
% Calculates the CC values for a given estimation.
%
% INPUTS:
% fixed_values      Values from the first iteration
% kfixed            Vector indicating which parameters are fixed
%                   Complete the rest of the parameters with NaN
%
% OUTPUT:
% CC                CC values from significance analysis
%
% Benjamín J. Sánchez
% Last Update: 2014-11-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CC = reg_analysis_onlyCC(fixed_values,kfixed)

assignin('base','kfixed',kfixed);

k = fixed_values;
k(find(~isnan(kfixed))) = [];

%Load experimental data:
texp    = evalin('base','texp');
yexp    = evalin('base','ydata');

%Eliminate rows with NaN in any measurement:
[N,n]    = size(yexp);
NaN_rows = zeros(size(texp))';
for i = 1:N
    for j = 1:n
        if isnan(yexp(i,j))
            NaN_rows(i) = 1;
        end
    end
end
NaN_rows = find(NaN_rows);
texp(NaN_rows,:) = [];
yexp(NaN_rows,:) = [];


%Jacobian and residual calculations using lsqcurvefit:
kL = k.*0.9999999999;
kU = k.*1.0000000001 + ones(1,length(k))*1e-10;
for i = 1:length(k)
    if k(i) < 0
        kL(i) = k(i)*1.0000000001;
        kU(i) = k(i)*0.9999999999;
    end
end
options = optimset('MaxIter',1,'MaxFunEvals',1,'TolFun',1,'TolX',1);

[k2,~,res,~,~,~,Jac] = lsqcurvefit(@minSquares2,k,texp,yexp,kL,kU,options);
%Confidence intervals:
[~,CC_SS]  = intconfianza(Jac,res,k2,0.05);
CC_SS      = CC_SS./100;
m          = length(k2);
for i = 1:m
    if k2(i) == 0
        CC_SS(i) = 0;
    end
end

%Adjust the results for easier integration with Excel:
m  = length(kfixed);
CC = zeros(1,m);
j  = 1;
for i = 1:m
    if isnan(kfixed(i))
        CC(i) = CC_SS(j);
        j     = j+1;
    end    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%