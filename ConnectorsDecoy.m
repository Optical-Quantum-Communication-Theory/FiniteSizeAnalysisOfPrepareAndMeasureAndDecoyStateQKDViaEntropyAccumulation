function [pre_keyrate,flag,crossover,acceptfrequency,optvalue,fst_opt] = ConnectorsDecoy (testprob,dimA,dimB,dimAprime,dimR,nu,eta,pd,depol,prob1,prob2,prob3,decoy1,decoy2,decoy3,epsilon_sec,f_EC,n,n_photon)
%Given the required inputs, calculates the single round quantity in EAT
%Inputs:
%       testprob:           initial choice of testing probability
%       dimA:               dimension of Alice's quantum system
%       dimB:               dimension of Bob's quantum system
%       dimAprime:          dimension of quantum system sent to Bob
%       dimR:               dimension of Alice's key register
%       nu:                 alpha-1 with alpha being the renyi parameter
%       eta:                channel loss
%       pd:                 dark count probability
%       depol:              depolarisation
%       prob1:              first decoy probability
%       prob2:              second decoy probability
%       prob3:              third decoy probability
%       decoy1:             first decoy intensity
%       decoy2:             second decoy intensity
%       decoy3:             third decoy intensity
%       epsilon_sec:        security epsilon
%       f_EC:               error correction prefactor
%       n:                  Number of signals
%       n_photon:           photon number cut-off
%Outputs:
%       pre_keyrate:        keyrate excluding the completeness correction
%       flag:               flag indicating final solver status
%       crossover:          gradient of crossover fmin
%       acceptfrequency:    expected honest (normalized) distribution on test-outcome alphabet

%generating vectors for decoy intensities and their probabilities
        decoys = [decoy1, decoy2, decoy3];
        probsDecoy = [prob1,prob2, prob3];
%generating EAT statistics
        [stateTestRounds,stateGenRounds,observablesEAT,acceptfrequency,krausOp,keyMap,~] = EATStatGenDecoy(testprob,dimA,dimB,eta,pd,depol,probsDecoy,decoys);
%finding optimal crossover
        [fst_opt,crossover,~,~,~] = FWminimiseRelEntChoi_DecoyV2(stateGenRounds,stateTestRounds,acceptfrequency,observablesEAT,keyMap,krausOp,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon);
%finding the optimal single round quantity
        [optvalue,~,flag] = ImprovedSecondOrderOptDecoy(crossover,stateGenRounds,stateTestRounds,acceptfrequency,observablesEAT,keyMap,krausOp,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon);
%optimal privacy amplification epsilon
        epsilon_PA = (epsilon_sec*(nu+1))/(1+2*nu);
%privacy amplification correction
        Privacy_amp_corr = PA_correction(nu,epsilon_PA,n);
%error correction
            er_cor = Error_corr_decoy (testprob,dimA,dimB,eta,pd,depol,probsDecoy,decoys,epsilon_sec,epsilon_PA,f_EC,n);
%keyrate calculation with all the corrections except completeness
            pre_keyrate = Keyrate(nu,crossover,optvalue,dimA)-Privacy_amp_corr-er_cor;
end