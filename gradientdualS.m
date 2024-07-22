function fval = gradientdualS(tplus,dimR,nu,testprob)
%Function Calculates the gradient of h* after relaxing it
%Input: tplus:              vector corresponding to the violations in the
%                           observations
%       dimR:               dimension of Alice's key register
%       nu:                 alpha-1, where alpha is the Renyi parameter
%       testprob:           testing probability
%Output: fval:              gradient of s

    %phi_0
    phi0 = 1/testprob*(nu)/(1-nu)*log(2)/2;

    %phi_1 
    phi1 = 1/sqrt(testprob)*(nu)/(1-nu)*log(2)/2*log2(2*dimR+1);

    if sum(tplus) >= 2*phi1
        fval = (sum(tplus) - 2*phi1)/(8*phi0)*ones(numel(tplus),1);
    else 
        fval = zeros(numel(tplus),1);
    end
end