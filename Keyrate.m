function l = Keyrate (nu,crossover,sv,dimA)
%Calculate keyrate using the single round quantity
%Inputs:
%       nu:             alpha-1 with alpha being the renyi parameter
%       crossover:      crossover min-tradeoff function
%       sv:             single round quantity
%       dimA:           dimension of Alice's quantum system

%calculate the maximum of crossover function
[~,max_g]=bounds(crossover,"all");

coef2 = ((nu)/(1-nu))^2;
ka_coef_num = (1-nu)^3;
ka_coef_den = 6*log(2)*((3-2*(nu+1))^3);
ka_coef = ka_coef_num/ka_coef_den;
ka_exponent = log2(dimA)+max_g;
frac = (nu)/(1-nu);
k_alpha1 = ka_coef*2^(frac*ka_exponent);
k_alpha2 = log((2^ka_exponent)+exp(2));
k_alpha = k_alpha1*(k_alpha2)^3;
l = (sv-coef2*k_alpha);
end