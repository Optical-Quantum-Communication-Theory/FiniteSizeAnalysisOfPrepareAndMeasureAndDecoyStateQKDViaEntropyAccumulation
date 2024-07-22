function [stateTestRounds,stateGenRounds,observablesEAT,acceptfrequency,krausOp,keyMap,expectationsCondGen] = EATStatGenDecoy(testprob,dimA,dimB,eta,pd,depol,probsDecoy,decoys)
%Function generates statistics and Kraus operators for EAT single round channel model
%Input:
%       testprob:           testing probability
%       dimA:               dimension of Alice's quantum system
%       dimB:               dimension of Bob's quantum system
%       eta:                channel loss
%       depol:              depolarisation

    pz = 1- testprob;
    px = testprob;

    %Assume perfect detectors
    etad = 1;

    %Only works for active!
    active = 1;

    %%%%%%%%%%%%%%%%%%%%% Construct POVM elements, Kraus ops and Key map %%%%%%%%%%%%%%%%%%%%%%%%%
    dimPB = 5;
    ketP = 1/sqrt(2)*[1;1];
    ketM = 1/sqrt(2)*[1;-1];
    
    %Alice's POVMs
    POVMsAGen = {diag([1,0]),diag([0,1])};
    POVMsATest = {ketP*ketP', ketM*ketM'};
    
    %POVMs for channel
    POVMsA = {pz*diag([1,0]),pz*diag([0,1]),(1-pz)*(ketP*ketP'),(1-pz)*(ketM*ketM')};
    
    %Bob's POVMs
    % include vacuum
    % block diagonal structure and order: 2x2 qubit, 1x1 vac
    POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};
    
    observablesJointTest = cell(numel(POVMsATest),numel(POVMsB));
    for indexA = 1:numel(POVMsATest)
        for indexB = 1:numel(POVMsB)
            observablesJointTest{indexA,indexB} = kron(POVMsATest{indexA},POVMsB{indexB});
        end
    end
    
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)    
    krausOpZ = kron((kron(zket(2,1),diag([1,0])) + kron(zket(2,2),diag([0,1])) ),diag([1,1,0])) ; % for Z basis
    
    %We use the Z basis only
    krausOp = {krausOpZ};
     
    %% key map
    proj0 = kron(diag([1,0]),eye(dimB*dimA));
    proj1 = kron(diag([0,1]),eye(dimB*dimA));
    keyMap = {proj0,proj1};
    
    % compute expectations
    probList = [pz/2,pz/2, px/2, px/2];

    for indexDecoy =1:numel(decoys)
        expectationsJoint(:,:,indexDecoy) = probsDecoy(indexDecoy)*diag(probList)*ExpectationsWCP(active,decoys(indexDecoy),eta,etad,pd,depol,pz);
        expectationsJointConditional(:,:,indexDecoy) = diag(probList)*ExpectationsWCP(active,decoys(indexDecoy),eta,etad,pd,depol,pz);
    end
    
    %Expectations conditioned on test rounds
    expectationsTest = 1/testprob*expectationsJoint(3:4,:,:);
  
    
    %%%%%%%%%%%%%%%%%%%%% States in Test and Gen rounds  and observables %%%%%%%%%%%%%%%%%%%%%%%%%
    
    vecGenRounds = 1/sqrt(2)*(kron([1; 0],[1;0]) + kron([0; 1],[0;1])) ;
    stateGenRounds = perturbationChannel(vecGenRounds*vecGenRounds');
    
    vecTestRounds = 1/sqrt(2)*(kron(ketP,ketP) + kron(ketM,ketM)) ;
    stateTestRounds = perturbationChannel(vecTestRounds*vecTestRounds');
    
    observablesEAT = observablesJointTest;
    acceptfrequency = expectationsTest;
    expectationsCondGen = expectationsJointConditional(:,:,1); %P(a,b|musig)
end
