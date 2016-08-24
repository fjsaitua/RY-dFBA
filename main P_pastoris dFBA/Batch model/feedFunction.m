%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F = feedFunction(t)
% Returns feed of fedbatch fermentation (or 0 if fermentation is batch)
%
% INPUT:
% t             Time of simulation [h]
% 
% OUTPUT:
% F             feed [L/h]
%
% All rights reserved to Benjamín J. Sánchez
%
% Last Update: 2015-02-06 - Franciso Saitua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = feedFunction(t)

% d = evalin('base','dataset'); %Dataset number

%if d == 1           %12-12-05 A
%    ti = 18.03333;
%    a  = 0.0010;
%    b  = 0.0620;
        
%elseif d == 2       %12-12-12 A
%    ti = 14.36667;
%    a  = 0.0011;
%    b  = 0.0806;
%    
%elseif d == 3       %12-12-17 A
%    ti = 13.33333;
%    a  = 0.0011;
%    b  = 0.0743;
%    
%elseif d == 5       %13-01-02 B
%    ti = 12.86667;
%    a  = 0.0009;
%    b  = 0.0620;
%    
%elseif d == 7       %13-01-17 A
%    ti = 13.28333;
%    a  = 0.0023;
%    b  = 0.0743;
%    
%elseif d == 9       %13-01-23 A
%    ti = 12.06667;
%    a  = 0.0022;
%    b  = 0.0664;
%    
%elseif d == 10      %13-01-23 B
%    ti = 12.10000;
%    a  = 0.0027;
%    b  = 0.0743;
%    
%elseif d == 12      %13-01-28 B
%    ti = 13.21667;
%    a  = 0.0020;
%    b  = 0.0664;
%    
%else                %Any anaerobic batch
    ti = 2000;
    a = 0.001;
    b = 0.0664;
%end

if t <= ti
    F = 0;
else
    F = a*exp(b*(t-ti));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%