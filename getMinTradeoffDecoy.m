function crossover = getMinTradeoffDecoy(stateTestRounds,observables,acceptfrequency,gradW,grads,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon)
    %Function finds the gradient of the crossover min trade-off function
    %Input: 
    %       stateGenRounds:     State in generation rounds
    %       stateTestRounds:    State in test rounds
    %       acceptfrequency:    Accept frequency in test rounds
    %       observables:        Observables of test rounds
    %       keyMap:             Kraus operators used for pinching channel Z
    %       krausOperators:     Kraus operators used for G map
    %       dimA:               dimension of Alice's quantum system
    %       dimAprime:          dimension of quantum system sent to Bob
    %       dimB:               dimension of Bob's quantum system
    %       dimR:               dimension of Alice's key register
    %       nu:                 alpha-1, where alpha is the Renyi parameter
    %       testprob:           testing probability
    %       decoys:             decoy intensities entered as a list [mu_signal, mu_2, mu_3,...]
    %       probsDecoy:         probabilities for decoys
    %       n_photon:           photon number cut-off
    %       gradW:              Gradient rel entropy
    %       gradS:              Gradient correction term
    %       optChoi:            optimal Choi matrix
    %
    %Output: 
    %        crossover:         gradient of the crossover min trade-off function

        %Dimension of Choi matrix
    dim = dimAprime*dimB;

    %number of observables
    nobs = numel(observables);

    %tolerance
    choi_tolerance = 1e-12;
    decoy_tolerance = 1e-12;

    %number of intensities used
    n_decoy=size(decoys,2);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(acceptfrequency,1);
    n = size(acceptfrequency,2);

    %Matrices to pick n=0,1,...-photon components
    M = {1,n_photon+1};
    for i = 1:n_photon+1
        Mi = zeros(m, m*(n_photon + 1));
        indr = [1:m]';
        indc = [i:n_photon+1:m*(n_photon + 1)]';
        indx = sub2ind(size(Mi),indr,indc);
        Mi(indx) = 1;
        M{i}=Mi;
    end

    %Calculate Alice's signal probabilities given a test round
    % uses equal probabilities of Alice's signal for each intensity!
    probsAliceTest = sum(acceptfrequency,[2,3]);

    cvx_solver mosek
    cvx_precision high
    cvx_begin sdp quiet
        variable Y(m*(n_photon+1),n) 
        variable J(dim,dim) hermitian semidefinite
        variable xplus(nobs*n_decoy) nonnegative
        variable lambda(m,n,n_decoy)
        variable YieldCutOff(m,n,n_decoy)
        dual variables duallambda{m,n,n_decoy}

        minimize real(trace(gradW'*J)) + transpose(grads)*xplus
        subject to

            %Density matrix in test rounds given Alice sent a single photon
            rhoTest = PartialTrace(kron(eye(dimA),J)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
            
            % 0 <= Y <= 1
            vec(Y) >= vec(zeros(m*(n_photon+1),n));
            vec(Y) <= vec(ones(m*(n_photon+1),n));

            % 1-photon component treated seperately with add.
            %constraints
            Y1 = M{2}*Y;

            % Require that yields sum to prob of Alice sending the signal
            for indexPhoton = 0: n_photon
                Yi = M{indexPhoton+1}*Y;
                sum(Yi,2) == probsAliceTest;
            end

            %0-photon yields
            Y0 = M{1}*Y;

            %0-photon error rate e_0=1/2
            Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));

            %sum(lammda) = 0
            ones(1,m*n*n_decoy)*vec(lambda) == 0; %decoy_tolerance;
            
            %Usual decoy bounds rewritten as matrix times vector
            
            for kdecoy = 1:n_decoy
                %Initialize probabilities
                pmu = Poisson(decoys(kdecoy),0:1:n_photon);
                Pmu = kron(eye(m),pmu);
                ptot = sum(pmu);

                %Compute sum_n p_mu(n) Y_n
                PmuY = Pmu*Y;
                
                %Contraints on delta^mu
                % 0<= delta^mu <= 1 - p_tot(mu)
                vec(YieldCutOff(:,:,kdecoy)) >= zeros(m*n,1);
                vec(YieldCutOff(:,:,kdecoy)) <= (1-ptot)*ones(m*n,1) + decoy_tolerance;
                
                %Remaining constraints on yields from observations and Choi
                %matrix
                for indexrow = 1:m
                    for indexcol = 1:n
                        Y1(indexrow,indexcol) - trace(observables{indexrow,indexcol}'*rhoTest) == 0;% decoy_tolerance;
                        duallambda{indexrow,indexcol,kdecoy} :0 == acceptfrequency(indexrow,indexcol,kdecoy) - probsDecoy(kdecoy)*(PmuY(indexrow,indexcol) + YieldCutOff(indexrow,indexcol,kdecoy)) - lambda(indexrow,indexcol,kdecoy);
                    end
                end
            end          
            
            %xplus >= lambda
            vec(xplus) >= abs(vec(lambda));

            %Partial trace constraint on J
            PartialTrace(J,[2],[dimAprime,dimB]) - eye(dimAprime) == zeros(dimAprime);

    cvx_end
    if strcmp(cvx_status, 'Infeasible') % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, subproblem status: %s ****\n",cvx_status);
    end
    
    %Dual variables are gradient of crossover min trade-off function.
    % Convert dual variables stored as cell to matrix.
    crossover = cell2mat(duallambda);
end

%%%%%%%%%%%%%%%%%%%%% Possonian distribution %%%%%%%%%%%%%%%%%%%%%%%%%
function prob = Poisson(mu,n)
    prob = exp(-mu).*mu.^n./factorial(n);
end