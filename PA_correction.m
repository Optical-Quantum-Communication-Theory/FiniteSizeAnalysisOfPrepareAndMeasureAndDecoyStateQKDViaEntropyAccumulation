function pac = PA_correction (nu,epsilon_PA,n)
%Privacy Amplification correction
%Inputs:
%       nu:             alpha-1 with alpha being the renyi parameter
%       epsilon_PA:     privacy amplification epsilon
%       n:              number of signals
first_correction = (nu+1/(nu))*log2(1/epsilon_PA);
pac = (first_correction-2)/n;
end