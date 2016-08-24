%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [AICc,CI,CC,Mc,Ms,diff] = reg_analysis(k)
% Performs al the regression analysis for a parameter set: significance,
% sensitivity and identifiability analysis, and AICc calculation.
% Returns all the different regression analysis outputs.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AICc,CI,CC,Mc,Ms,diff] = reg_analysis_complete(k)

assignin('base','skip_delays',false);

%Load experimental data:
texp    = evalin('base','texp');
yexp    = evalin('base','ydata');
weights = evalin('base','weights');

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

simTime = texp(length(texp));

% %Normalize data with maximum measures:
% for i = 1:8
%     yexp(:,i) = yexp(:,i)./(max(yexp(:,i))*weights(:,i));
% end

%Jacobian and residual calculations using lsqcurvefit:
kL = k.*0.9999999999;
kU = k.*1.0000000001 + ones(1,length(k))*1e-15; % ojo cambiar
for i = 1:length(k)
    if k(i) < 0
        kL(i) = k(i)*1.0000000001;
        kU(i) = k(i)*0.9999999999;
    end
end

options = optimset('MaxIter',1,'MaxFunEvals',1,'TolFun',1,'TolX',1);

try
    [k2,~,res,~,~,~,Jac] = lsqcurvefit(@minSquares2,k,texp,yexp,kL,kU,options);
    diff = sum((k-k2).^2);
    %Confidence intervals:
    [CI,CC] = intconfianza(Jac,res,k2,0.05);
    CC      = CC./100;
    m       = length(k2);
    disp('Unsignificant Parameters:');
    for i = 1:m
        if k2(i) == 0
            CC(i) = 0;
        elseif CC(i) >= 2
            disp(num2str(i));
        end
    end
    
    %AIC calculation:
    SS   = sum(sum(res.^2));
    AIC  = N*log(SS)+2*(m+1);
    AICc = AIC + 2*(m+1)*(m+2)/(N-m-2);
    
catch exception
    disp(exception.identifier);
    CI   = zeros(size(k));
    CC   = zeros(size(k));
    diff = 0;
    AICc = 0;
end

%Identifiability analysis:
x0 = evalin('base','x0');
post_starv = 0;
assignin('base','post_starv',post_starv)
Mc = identificaBSB(1:7,k,x0,[0 simTime],@pseudoSteadyState,0.95,1);

%Sensitivity analysis:
close all
post_starv = 0;
assignin('base','post_starv',post_starv)
ksensibilidadBSB(1:7,k,x0,[0 simTime],@pseudoSteadyState,3,1);
Ms = load('GraficoBarra.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%