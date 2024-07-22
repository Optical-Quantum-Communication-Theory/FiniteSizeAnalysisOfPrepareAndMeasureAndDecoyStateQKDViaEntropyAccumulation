function [pre_keyrate,optdeltacomp_temp] = Connectors (testprob,dimA,dimB,dimAprime,dimR,nu,eta,depol,epsilon_sec,f_EC,n)
%Given the required inputs, calculates the single round quantity in EAT
%Inputs:
%       testprob:           initial choice of testing probability
%       dimA:               dimension of Alice's quantum system
%       dimB:               dimension of Bob's quantum system
%       dimAprime:          dimension of quantum system sent to Bob
%       dimR:               dimension of Alice's key register
%       nu:                 alpha-1 with alpha being the renyi parameter
%       eta:                channel loss
%       depol:              depolarisation
%       epsilon_sec:        security epsilon
%       f_EC:               error correction efficiency
%       n:                  Number of signals

%Generating EAT statistics
        [~,stateTestRounds,stateGenRounds,observablesEAT,acceptfrequency,krausOp,keyMap] = EATStatGen(testprob,dimA,dimB,eta,depol);
%Finding the optimal crossover min-tradeoff function
        [~,crossover,~,~,~] = FWminimiseRelEntChoi(stateGenRounds,stateTestRounds,acceptfrequency,observablesEAT,keyMap,krausOp,dimA,dimAprime,dimB,dimR,nu,testprob);
%Finding the optimal value of the FW
        [optvalue,~,~] = ImprovedSecondOrderOpt(crossover,stateGenRounds,stateTestRounds,acceptfrequency,observablesEAT,keyMap,krausOp,dimA,dimAprime,dimB,dimR,nu,testprob);
%Calculating the optimal epsilon_PA
        epsilon_PA = (epsilon_sec*(nu+1))/(1+2*nu);
%Privacy amplification Correction
        Privacy_amp_corr = PA_correction(nu,epsilon_PA,n);
%Error correction correction
        er_cor = Error_corr (testprob,dimA,dimB,eta,depol,epsilon_sec,epsilon_PA,f_EC,n);
%Completeness parameter
        if abs(sum(acceptfrequency(:))-1) > 1e-10 || min(sum(acceptfrequency(:))) < 0
            pre_keyrate = -2;
            optdeltacomp_temp = -1;
        else
        [~,optdeltacomp_temp , ~ , ~ , ~ , ~ , ~] = optdeltacom(n,testprob,acceptfrequency,crossover,1e-3,100);
%Keyrate calculation with all the corrections
        pre_keyrate = Keyrate (nu,crossover,optvalue,dimA)-Privacy_amp_corr-er_cor-optdeltacomp_temp;
        end

end