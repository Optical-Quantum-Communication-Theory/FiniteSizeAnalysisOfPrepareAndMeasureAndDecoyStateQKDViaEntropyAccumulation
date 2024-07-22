function ec = Error_corr (testprob,dimA,dimB,eta,depol,epsilon_sec,epsilon_PA,f_EC,n)
%Calculate the error correction 
%Inputs:
%       testprob:           initial choice of testing probability
%       dimA:               dimension of Alice's quantum system
%       dimB:               dimension of Bob's quantum system
%       eta:                channel loss
%       depol:              depolarisation
%       epsilon_sec:        security epsilon
%       epsilon_PA:         privacy amplification epsilon
%       f_EC:               error correction efficiency
%       n:                  Number of signals

%Generating the joint expectations
[expectationsJoint,~,~,~,~,~,~] = EATStatGen(testprob,dimA,dimB,eta,depol);
    %Error validation epsilon
    epsilon_EV = epsilon_sec - epsilon_PA;
    %calculate expectation values of generation rounds
    expectationsGen= expectationsJoint(1:2,1:2)*1/sum(expectationsJoint(1:2,1:2),"all");
    %calculate the joint entropy
    entjoint = Entropy(expectationsGen);
    %calculate the reduced state in B register
    B = zeros(2,2);
    B(1,1) = sum(expectationsGen(:,1));
    B(2,2) = sum(expectationsGen(:,2));
    %calculate the entropy of Bob
    entB = Entropy(B);
    %calculate conditional entropy
    entCond = entjoint - entB;
    %probability of detection
    Pdet = 1 - sum(expectationsJoint(:,5),"all");
    %other corrections
    pre = Pdet*(1-testprob)^2;
    cor = (1/n)*log2(2/epsilon_EV);
    %calculate error correction 
    ec = f_EC*pre*entCond + cor;
end