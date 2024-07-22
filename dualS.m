function fval = dualS(tplus,dimR,nu,testprob)
%Function Calculates the value of h* after relaxing it
%Input: tplus:              vector corresponding to the violations in the
%                           observations
%       dimR:               dimension of Alice's key register
%       nu:                 alpha-1, where alpha is the Renyi parameter
%       testprob:           testing probability
%Output: fval:              value of s

    %phi_0
    phi0 = 1/testprob*(nu)/(1-nu)*log(2)/2;

    %phi_1 
    phi1 = 1/sqrt(testprob)*(nu)/(1-nu)*log(2)/2*log2(2*dimR+1);

    if sum(tplus) >= 2*phi1
        fval = ((sum(tplus) - 2*phi1)^2)/(16*phi0);
    else 
        fval = 0;
    end
end