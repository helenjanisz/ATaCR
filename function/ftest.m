function [ P ] = ftest( res1,parms1,res2,parms2 )
%[ P ] = ftest( res1,parms1,res2,parms2 )
%
%
% function to perform a formal F-test on two sets of residuals of fits to
% data (res1, res2) that have respectively been fitted using parms1,parms2
% (where these are integers equal to the number of parameters used in each
% fit)
% IMPORTANTLY: 
%   This will test whether the second fit (yielding res2) is statistically
%   superior to the first fit, but used more parms. I.e.:
%          sum(abs(res2)) < sum(abs(res1)) 
%   and            parms2 > parms1
%
% the degrees of freedom for each fit are therefore:
% 1) N - parms1
% 2) N - parms2
% where N is the number of data, or length(res*)
% The residuals are just equal to dobs-dpred, so we square and sum these to
% get the chi^2 values (Ea)
% 
% Z. Eilon

N1 = length(res1);
N2 = length(res2);
% generally N1 and N2 will be equal

% calculate degrees of freedom
v1 = N1-parms1;
v2 = N2-parms2;

% calculate chi^2 sums
Ea_1 = sum(res1.^2);
Ea_2 = sum(res2.^2);

Fobs = (Ea_1/v1)/(Ea_2/v2);

P = 1 - ( fcdf(Fobs,v1,v2) - fcdf(1/Fobs,v1,v2) );

end

