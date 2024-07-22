function [expectationsJoint,stateTestRounds,stateGenRounds,observablesEAT,acceptfrequency,krausOp,keyMap] = EATStatGen(testprob,dimA,dimB,eta,depol)
%Function generates statistics and Kraus operators for EAT single round channel model
%Input:
%       testprob:           testing probability
%       dimA:               dimension of Alice's quantum system
%       dimB:               dimension of Bob's quantum system
%       eta:                channel loss
%       depol:              depolarisation

    pz = 1- testprob;
    px = testprob;

    %%%%%%%%%%%%%%%%%%%%% Construct POVM elements, Kraus ops and Key map %%%%%%%%%%%%%%%%%%%%%%%%%
    dimPB = 5;
    ketP = 1/sqrt(2)*[1;1];
    ketM = 1/sqrt(2)*[1;-1];
    
    %Alice's POVMs
    POVMsAGen = {diag([1,0]), diag([0,1])};
    POVMsATest = {ketP*ketP', ketM*ketM'};
    
    %POVMs for channel
    POVMsA = {pz*diag([1,0]),pz*diag([0,1]),(1-pz)*(ketP*ketP'),(1-pz)*(ketM*ketM')};
    
    %Bob's POVMs
    % includes vacuum
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
    
    
    % generate rho and apply quantum channel
    
    observablesJointChannel = cell(numel(POVMsA),numel(POVMsB));
    for indexA = 1:numel(POVMsA)
        for indexB = 1:numel(POVMsB)
            observablesJointChannel{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
        end
    end
    
    
    % generate the maximally entangled density matrix for these 
    rhoAB = 0;
    for index = 1:dimA
        rhoAB = rhoAB + kron(zket(dimA,index),zket(dimA,index));
    end
    
    rhoAB = (rhoAB*rhoAB')/dimA;
    rhoAB = depolarizationChannel(rhoAB,depol);
    
    rhoAB = lossChannel(rhoAB, eta);
    
    % compute expectations
    % load observables from description to do so
    expectationsJoint = zeros(size(observablesJointChannel));
    
    for index = 1:numel(observablesJointChannel)
        expectationsJoint(index) = trace(observablesJointChannel{index}'*rhoAB);
    end
    
    expectationsTest = 1/testprob*expectationsJoint(3:4,:);
    
    %%%%%%%%%%%%%%%%%%%%% States in Test and Gen rounds  and observables %%%%%%%%%%%%%%%%%%%%%%%%%
    
    vecGenRounds = 1/sqrt(2)*(kron([1; 0],[1;0]) + kron([0; 1],[0;1])) ;
    stateGenRounds = perturbationChannel(vecGenRounds*vecGenRounds');
    
    vecTestRounds = 1/sqrt(2)*(kron(ketP,ketP) + kron(ketM,ketM)) ;
    stateTestRounds = perturbationChannel(vecTestRounds*vecTestRounds');
    
    observablesEAT = observablesJointTest;
    acceptfrequency = expectationsTest;
end

%%%%%%%%%%%%%%%%%%%%% Lossy Channel %%%%%%%%%%%%%%%%%%%%%%%%%
% Lossy channel takes in a 2x2 qubit and outputs a 3x3 density matrix formed by qubit + perp
function rhoPrime = lossChannel(rho,eta)
    loss = 1-eta;
    krausOps = cell(1,3);
    % mapping to loss
    krausOps{1} = sqrt(loss)*zket(3,3)*zket(2,1)';
    krausOps{2} = sqrt(loss)*zket(3,3)*zket(2,2)';
    % remainder is not lost
    V = [1,0;0,1;0,0];
    krausOps{3} = sqrt(1-loss)*V;
    %Add alice’s system
    for i = 1:length(krausOps)
        krausOps{i} = kron(eye(2),krausOps{i});
    end
    rhoPrime = krausFunc(rho,krausOps);
end

%%%%%%%%%%%%%%%%%%%%% Depolarisation Channel %%%%%%%%%%%%%%%%%%%%%%%%%

% Applies a depolarization channel to rho
function rhoPrime = depolarizationChannel(rho, depol) 
    % Nielsen & Chuang p. 379
    krausOps = cell(1,4);
    pauliMatrices = {eye(2), [0,1;1,0], [0,1i;-1i,0], [1,0;0,-1]};
    krausOps{1} = kron(eye(2), sqrt(1-3/4*depol)*eye(2));
    for i = 2 : length(pauliMatrices)
        krausOps{i} = kron(eye(2), sqrt(depol)*pauliMatrices{i}/2);
    end
    rhoPrime = krausFunc(rho, krausOps);
end
