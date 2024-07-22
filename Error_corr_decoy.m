function ecd = Error_corr_decoy (testprob,dimA,dimB,eta,pd,depol,probsDecoy,decoys,epsilon_sec,epsilon_PA,f_EC,n)
%Calculate the error correction 
%Inputs:
%       testprob:           initial choice of testing probability
%       dimA:               dimension of Alice's quantum system
%       dimB:               dimension of Bob's quantum system
%       eta:                channel loss
%       pd:                 dark count probability
%       depol:              depolarisation
%       probsDecoy          decoys probabilities
%       decoys:             decoys intensities
%       epsilon_sec:        security epsilon
%       epsilon_PA:         privacy amplification epsilon
%       f_EC:               error correction prefactor
%       n:                  Number of signals

%Regenerating EAT statistics
    [~,~,~,~,~,~,expectationsJoint] = EATStatGenDecoy(testprob,dimA,dimB,eta,pd,depol,probsDecoy,decoys);
%epsilon evaluation
    epsilon_EV = epsilon_sec - epsilon_PA;
%correction prefactors
    prefactors = (1-testprob)^2*probsDecoy(1);
%expectation value of generation rounds
    expectationsGen=expectationsJoint(1:2,1:2)/sum(expectationsJoint(1:2,1:2),'all');
%calculate the joint entropy
    entjoint = Entropy(expectationsGen);
%calculate the reduced state on B register
    B = zeros(2,2);
    B(1,1) = sum(expectationsGen(:,1));
    B(2,2) = sum(expectationsGen(:,2));
%entropy of Bob
    entB = Entropy(B);
%calculate conditional entropy
    entCond = entjoint - entB;
%detection probability    
    Pdet = 1 - sum(expectationsJoint(:,5),"all");
%other corrections
    pre = Pdet*(1-testprob)^2;
    cor = (1/n)*log2(2/epsilon_EV);
%calculate error correction 
    ecd = f_EC*pre*entCond + cor;
end
