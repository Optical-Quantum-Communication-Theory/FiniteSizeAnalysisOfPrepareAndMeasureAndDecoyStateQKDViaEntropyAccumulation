function [lowerbnd_final,upperbnd,flag] = ImprovedSecondOrderOptDecoy(crossover,stateGenRounds,stateTestRounds,acceptfrequency,observables,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon)
%Function minimises relative entropy with penalty term h* and finds the
%gradient of the crossover min trade-off function
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
%
%Output: lowerbnd_final:    final lower bound on optimal value of optimisation
%        upperbnd:          final upper bound on optimal value of optimisation
%        flag:              flag indicating final solver status
    
    %%%%%%%%%%%%%%%%%%%%% Parameters for FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%
    %Enable printing
    verbose = 1;

    %Store status of FW steps
    status = [];
       
    %Maximum iteratins
    maxiter = 150;
    
    %Tolerance in checking vertices are in active set
    tolvert = 1e-10;
    
    %search precision for line step
    linesearchprecision = 1e-20;
    
    %Tolerance in solution
    maxgap = 1e-6;

    %%%%%%%%%%%%%%%%%%%%% Initial guess %%%%%%%%%%%%%%%%%%%%%%%%%

    %Initial guess using Choi Matrix
    J0 = eye(dimAprime*dimB)/dimAprime;

    [Jinit, Yieldinit, YieldCutOffinit] = closestChoiMatrixDecoy(J0,stateTestRounds,observables,acceptfrequency,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon);  
        
    %%%%%%%%%%%%%%%%%%%%% Prepare FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

    %Initialise initial guess for FW iteration
    Jmat = Jinit;
    Yield = Yieldinit;
    YieldCutOff = YieldCutOffinit;
    
    %Coefficients of vectors in active set
    alpha_t = {1};
    
    %Indices indicating active vectors in S_t (Although S_t is the active
    %set, vectors are not removed only their indices in I_active)
    I_active = [1];
    
    %Initialize active set
    S_t = {{Jmat,Yield,YieldCutOff}};

    %Keep track of gammas
    gamma_t_list = [];
    
    %%%%%%%%%%%%%%%%%%%% Run FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Perform iteration
    if verbose == 1
        fprintf("Performing FW Iterations \n")
    end

    for t = 1:maxiter
        if verbose == 1        
            fprintf('FW iteration:%d',t)
        end
        tstart_FW=tic;
    
        is_away = false;
        
        %Gradients
        gradW = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*gradientRelEntChoi(Jmat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
        gradCorYield = gradientCorYield(Yield, YieldCutOff, crossover,testprob,nu,dimR,probsDecoy,decoys,n_photon);
        gradCorYieldCutOff = gradientCorYCutOff(Yield, YieldCutOff, crossover,testprob,nu,dimR,probsDecoy,decoys,n_photon);
        
        %Find step direction for FW step
        %solve for FW step
        [sJt, sYieldt, sYieldCutOfft, status] = subproblem(stateTestRounds,observables,acceptfrequency,gradW,gradCorYield,gradCorYieldCutOff,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon,status);
              
        %step direction in FW case
        dJt_FW = sJt - Jmat;
        dYieldt_FW = sYieldt - Yield;
        dYieldCutOfft_FW = sYieldCutOfft - YieldCutOff;
    
    
        %Find step direction for Away step
        %Evaluate <grad(f), v> for all v in S_t
        vtcell = cellfun(@(x) trace(gradW'*x{1}) + evaluateGradYdotY(gradCorYield,x{2}) + transpose(gradCorYieldCutOff)*vec(x{3}), S_t, 'Uniform', 0);
    
        %convert to array
        vtmat = cell2mat(vtcell);
    
        %maximise only over active vectors
        [~,id] = max(vtmat(I_active));
    
        %save index of maximum
        index_away = I_active(id);
        
        %vector reaching maximum of <grad(f), v> for all v in S_t
        vJt = S_t{index_away}{1};
        vYieldt = S_t{index_away}{2};
        vYieldCutOfft = S_t{index_away}{3};

        %Pair direction
        dJt = sJt - vJt;
        dYieldt = sYieldt - vYieldt;
        dYieldCutOfft = sYieldCutOfft - vYieldCutOfft;
        gamma_max = alpha_t{index_away};
        
    
        % perform an exact line search for step size
        optimoptions = optimset('TolX',linesearchprecision,'Tolfun',linesearchprecision);
        minfunc = @(delta) primalf(Jmat + delta*dJt, Yield + delta*dYieldt, YieldCutOff + delta*dYieldCutOfft, crossover, acceptfrequency,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon);
        %step size given by solution
        gamma_t = fminbnd(minfunc,0,gamma_max,optimoptions);     
    

        %Perform updates of active set S_t, alpha_t and set of active
        %indices I_active
    
        % Away part:
        if abs(gamma_t - gamma_max) < 1e-10
                if verbose == 1 
                    fprintf("  Drop step")
                end
    
                %If a drop step happens the coeficient alpha is set to 0
                alpha_t{index_away} = 0;
                %and the according vector is set inactive by removing it
                %from the set of active indices
                I_active(I_active == index_away) = [];
        else
            alpha_t{index_away} = alpha_t{index_away} - gamma_t;
        end
    
        %Towards part
    
        % Is the step st a new vector?
        %transform S_t into cell array of vectors
        S_t_cellvec = cellfun(@(x) stackentriesCell(x),S_t,'Uniform', 0);
        FWindex = find(ismembertol([S_t_cellvec{:}]',[reshape(real(sJt),[],1); reshape(imag(sJt),[],1) ; reshape(sYieldt,[],1); reshape(sYieldCutOfft,[],1)]',tolvert,'ByRows',true,'DataScale',1));
    
        %If st is not in S_t
        if isempty(FWindex) 
            %add st to S_t
            S_t{end+1} = {sJt, sYieldt, sYieldCutOfft};
    
            %add step size gamma_t as another component to alpha_t
            alpha_t{end+1} = gamma_t;
    
            %Add new vector as an active vector to the set of active
            %vectors (will be last entry in S_t, thus numel(S_t) added
            %to I_active)
            I_active = [I_active, numel(S_t)];
        
        %If st is in S_t
        else
            alphaMatFW = [alpha_t{FWindex}];
            
            if isempty(alphaMatFW(alphaMatFW < 1e-13)) ~= 1 %vector was in St but not active
                %Activate index corresponding to vector st
                I_active = [I_active, FWindex];
            end
    
            %vector was in S_t and the coefficient corresponding to st is updated
            %accordingly by only adding gamma_t since the
            %multiplication was done before
            alpha_t(FWindex) = cellfun(@(x) x + gamma_t ,alpha_t(FWindex),'Uniform', 0); %alpha_t{FWindex} + gamma_t;
        end
        
        %Exceptional case if gamma is 1, the set S_t collapses to 
        % S_t = {st}, here done by removing all indices from set of
        % active indices and replacing them by index of st
        if gamma_t > 1- (1e-13)
            I_active = [FWindex];
        end
        
        %Calculate next step of optimisation variable
        Jnew = Jmat + gamma_t*dJt;
        YieldNew = Yield + gamma_t*dYieldt;
        YieldCutOffNew = YieldCutOff + gamma_t*dYieldCutOfft;
        
        %Calculate duality gap        
        gap = -trace(gradW' * dJt_FW) - evaluateGradYdotY(gradCorYield, dYieldt_FW) - transpose(gradCorYieldCutOff) * vec(dYieldCutOfft_FW);
    
        %Calculate projection value        
        proj_val = trace(gradW' * Jnew) + evaluateGradYdotY(gradCorYield, YieldNew) + transpose(gradCorYieldCutOff) * vec(YieldCutOffNew);
    
        %Calculate function value of new and old optimisation variable 
        fNew = primalf(Jnew, YieldNew, YieldCutOffNew, crossover, acceptfrequency,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon);
        fOld = primalf(Jmat, Yield, YieldCutOff, crossover, acceptfrequency,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon);
    
        %Calculate lower bound on minimum
        lowerbnd = fOld + trace(gradW' * dJt_FW) + evaluateGradYdotY(gradCorYield, dYieldt_FW) + transpose(gradCorYieldCutOff) * vec(dYieldCutOfft_FW);
    
        %print progress
        t_FW=toc(tstart_FW);
        fprintf('\n    FW iteration time: %f\n',t_FW);
        fprintf('lower bound:%f     projection value:%f    gap:%f    fvalue:%f    fgap:%f\n',lowerbnd,proj_val,gap,fNew,fOld-fNew)
        
        %Stopping criterion if gap is small enough
        if  ( abs(gap) < maxgap)
            Jmat = Jmat + gamma_t*dJt;
            Yield = Yield + gamma_t*dYieldt;
            YieldCutOff = YieldCutOff + gamma_t*dYieldCutOfft;
            break;
        end
        
        %Update optimisation variables
        Jmat = Jmat + gamma_t*dJt;
        Yield = Yield + gamma_t*dYieldt;
        YieldCutOff = YieldCutOff + gamma_t*dYieldCutOfft;
    
        %Display warning if algorithm hasn't converged and reached the
        %maximum iterations
        if t == maxiter
            disp('**** Warning: Maximum iteration limit reached. ****');
        end
    end

    %Fix Choi matrix to be hermitian
    J_opt = (Jmat +Jmat')/2;
    Yieldopt = Yield;
    YieldCutOffopt = YieldCutOff;

    %Calculate rho from optimal Choi matrix
    rho_opt = PartialTrace(kron(eye(dimA),J_opt)*(kron(PartialTranspose(stateGenRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
    
    % check size of rho
    if isempty(krausOperators)
        dprime = size(rho_opt, 1);
    else
        dprime = size(krausFunc(rho_opt, krausOperators), 1);
    end

    %Calculate perturbation needed for rho
    perturbation = computePerturbationEpsilon(rho_opt, krausOperators, keyMap);
    
    if perturbation > 1/(exp(1)*(dprime-1))
        ME = MException('FW2StepSolver:epsilon too large','Theorem cannot be applied. Please have a better rho to start with');
        throw(ME);
    end

    %Calculate gradient with fixed perturbation
    gradW_eps = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*gradientRelEntChoi_eps(perturbation,J_opt,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
    gradCorYield = gradientCorYield(Yield, YieldCutOff, crossover,testprob,nu,dimR,probsDecoy,decoys,n_photon);
    gradCorYieldCutOff = gradientCorYCutOff(Yield, YieldCutOff, crossover,testprob,nu,dimR,probsDecoy,decoys,n_photon);
    
    %Run one more FW iteration at optimal Choi matrix and perturbation for
    %gradient dJt_final
    statusfinal = [] ;
    [sJt_final, sYieldt_final, sYieldCutOfft_final, statusfinal] = subproblem(stateTestRounds,observables,acceptfrequency,gradW_eps,gradCorYield,gradCorYieldCutOff,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon,statusfinal);
    
    dJt_final = sJt_final - J_opt;
    dYieldt_final = sYieldt_final - Yieldopt;
    dYieldCutOfft_final = sYieldCutOfft_final - YieldCutOffopt;
    
    %Function value at correctly perturbed optimal points
    f_final = primalf_eps(perturbation,J_opt,Yieldopt, YieldCutOffopt,crossover,acceptfrequency,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon);
              
    %Calculate correction term due to perturbation
    if perturbation == 0
        zetaEp = 0;
    else
        zetaEp = 2*perturbation*(dprime-1)*log(dprime/(perturbation*(dprime-1)));
    end
    
    upperbnd = f_final;
    lowerbnd_final = f_final + trace(gradW_eps' * dJt_final) + evaluateGradYdotY(gradCorYield, dYieldt_final) + transpose(gradCorYieldCutOff) * vec(dYieldCutOfft_final) - zetaEp;
    flag = statusfinal;
end

%%%%%%%%%%%%%%%%%%%%% Objective functions of FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

%Objective function in FW iteration with perturbation calculated at each
%iteration

function fval = primalf(ChoiMat,Yield,YieldCutOff,crossover,acceptfrequency,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon)
    m = size(crossover,1);
    n = size(crossover,2);
    n_decoy = size(crossover,3);

    genfrequency = zeros(m,n,n_decoy);

    for indexDecoy = 1:n_decoy
        pmu = Poisson(decoys(indexDecoy),0:1:n_photon);
        Pmu = kron(eye(m),pmu);
        PmuY = Pmu*Yield;

        genfrequency(:,:,indexDecoy) = probsDecoy(indexDecoy)*(PmuY + YieldCutOff(:,:,indexDecoy));
    end
    fval = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*relEntropyChoi(ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB) + sum(crossover.*(acceptfrequency - genfrequency),'all') - nu/(1-nu)*log(2)/2*Varq(genfrequency,crossover,testprob,dimR);
end

%Objective function with fixed perturbation for final calculation

function fval = primalf_eps(perturbation,ChoiMat,Yield,YieldCutOff,crossover,acceptfrequency,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon)
    m = size(crossover,1);
    n = size(crossover,2);
    n_decoys = size(crossover,3);

    genfrequency = zeros(m,n,n_decoys);

    for indexDecoy = 1:n_decoys
        pmu = Poisson(decoys(indexDecoy),0:1:n_photon);
        Pmu = kron(eye(m),pmu);
        PmuY = Pmu*Yield;

        genfrequency(:,:,indexDecoy) = probsDecoy(indexDecoy)*(PmuY + YieldCutOff(:,:,indexDecoy));
    end
    
    fval = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*relEntropyChoi_eps(perturbation,ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB) + sum(crossover.*(acceptfrequency - genfrequency),'all') - nu/(1-nu)*log(2)/2*Varq(genfrequency,crossover,testprob,dimR);
end


function fval = Varq(genfrequency,crossover,testprob,dimR)
    varp = 1/testprob*sum(genfrequency.*(max(crossover, [],'all') - crossover).^2,'all') - (max(crossover, [],'all') - sum(crossover.*genfrequency,'all'))^2;
    fval = (log2(1+2*dimR) + sqrt(2 + varp))^2;
end

%%%%%%%%%%%%%%%%%%%%% Subproblem calculating deltaJ %%%%%%%%%%%%%%%%%%%%%%%%%

function [Jtilde, Yield, YCutOff, status] = subproblem(stateTestRounds,observables,acceptfrequency,gradW,gradCorY,gradCorYCutOff,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon,status)
    %Dimension of Choi matrix
    dim = dimAprime*dimB;

    %number of observables
    nobs = numel(observables);

    %tolerance
    decoy_tolerance = 0;

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
        variable YieldCutOff(m,n,n_decoy) 

        sumYieldObjective = 0;
        
        %Construct sum of all yield components used in objective funstion
        for indexPhoton = 0:1:n_photon
            %Yields for n photons as vector
            vecYn = vec(M{indexPhoton+1}*Y);
            %gradient for yields of n photons as vector
            vecgradCorYn = vec(gradCorY{indexPhoton+1});
    
            %add gradYn*Yn to objective
            sumYieldObjective = sumYieldObjective + transpose(vecgradCorYn)*vecYn;
        end

        minimize real(trace(gradW'*J)) + transpose(gradCorYCutOff)*vec(YieldCutOff) + sumYieldObjective
        subject to

            %Density matrix in test rounds given Alice sent a single photon
            rhoTest = PartialTrace(kron(eye(dimA),J)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
            
            % 0 <= Y <=1
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

            %Usual decoy bounds rewritten as matrix times vector
            
            for k = 1:n_decoy
                %Initialize probabilities
                pmu = Poisson(decoys(k),0:1:n_photon);
                Pmu = kron(eye(m),pmu);

                %Compute sum_n p_mu(n) Y_n
                PmuY = Pmu*Y;

                %Contraints on delta^mu
                % 0 <= delta^mu <= 1 - p_tot(mu)
                vec(YieldCutOff(:,:,k)) >= zeros(m*n,1);
                vec(YieldCutOff(:,:,k)) <= (1-sum(pmu))*ones(m*n,1) + decoy_tolerance;
                
                %Remaining constraints on yields from observations and Choi
                %matrix
                for p = 1:m
                    for q = 1:n
                        %Single Photon
                        Y1(p,q) - trace(observables{p,q}'*rhoTest) == 0;
                    end
                end
            end           

            %Partial trace constraint on J
            PartialTrace(J,[2],[dimAprime,dimB]) - eye(dimAprime) == zeros(dimAprime);

    cvx_end
    if strcmp(cvx_status, 'Infeasible') % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, subproblem status: %s ****\n",cvx_status);
    end

    %record status for debugging
    status = [status, string(cvx_status)];
    
    %Output results from optimisation
    Jtilde = (full(J)+full(J)')/2;
    Yield = full(Y);
    YCutOff = full(YieldCutOff);
end

%%%%%%%%%%%%%%%%%%%%% Closest Choi matrix %%%%%%%%%%%%%%%%%%%%%%%%%

function [Jtilde, Yield, YCutOff] = closestChoiMatrixDecoy(J0,stateTestRounds,observables,acceptfrequency,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon)
    %Dimension of Choi matrix
    dim = dimAprime*dimB;

    %number of observables
    nobs = numel(observables);

    %tolerance
    decoy_tolerance = 0;

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
        variable YieldCutOff(m,n,n_decoy)

        minimize norm(J - J0)
        subject to

            %Density matrix in test rounds given Alice sent a single photon
            rhoTest = PartialTrace(kron(eye(dimA),J)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
            
            % 0 <= Y <=1
            vec(Y) >= vec(zeros(m*(n_photon+1),n));
            vec(Y) <= vec(ones(m*(n_photon+1),n));

            %0- and 1-photon component treated seperately with add.
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
  
            %Usual decoy bounds rewritten as matrix times vector
            
            for k = 1:n_decoy
                %Initialize probabilities
                pmu = Poisson(decoys(k),0:1:n_photon);
                Pmu = kron(eye(m),pmu);

                %Compute sum_n p_mu(n) Y_n
                PmuY = Pmu*Y;

                %Contraints on delta^mu
                % 0 <= delta^mu <= 1 - p_tot(mu)
                vec(YieldCutOff(:,:,k)) >= zeros(m*n,1);
                vec(YieldCutOff(:,:,k)) <= (1-sum(pmu))*ones(m*n,1) + decoy_tolerance;
                
                %Remaining constraints on yields from observations and Choi
                %matrix
                for p = 1:m
                    for q = 1:n
                        %Single Photon
                        Y1(p,q) - trace(observables{p,q}'*rhoTest) == 0;
                    end
                end
            end

            %Additional constraints for 0- and 1-photon components
            %in terms of Choi matrices
           

            %Partial trace constraint on J
            PartialTrace(J,[2],[dimAprime,dimB]) - eye(dimAprime) == zeros(dimAprime);

    cvx_end
    if strcmp(cvx_status, 'Infeasible') % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, subproblem status: %s ****\n",cvx_status);
    end
    
    %Output results from optimisation
    Jtilde = (full(J)+full(J)')/2;
    Yield = full(Y);
    YCutOff = full(YieldCutOff);
end

%%%%%%%%%%%%%%%%%%%%% Compute needed perturbation %%%%%%%%%%%%%%%%%%%%%%%%%

function epsilon = computePerturbationEpsilon(rho, krausOps, keyMap)
% This function uses perturbationChannelEpsilon to compare the epsilons
% obtained from G(rho) and Z(G(rho)) and takes the maximum of the two.
    % first check G(rho)
    gRho = krausFunc(rho, krausOps);
    epsilonG = perturbationChannelEpsilon(gRho);

    zRho = krausFunc(gRho, keyMap);
    epsilonZ = perturbationChannelEpsilon(zRho);

    epsilon = max([epsilonG, epsilonZ]);
end

%%%%%%%%%%%%%%%%%%%%% Possonian distribution %%%%%%%%%%%%%%%%%%%%%%%%%
function prob = Poisson(mu,n)
    prob = exp(-mu).*mu.^n./factorial(n);
end

%%%%%%%%%%%%%%%%%%%%% Helper function evaluating gradY.Y %%%%%%%%%%%%%%%%%%%%%%%%%

function fval = evaluateGradYdotY(gradCorY,Y)
    n_photon = size(gradCorY,2) -1;
    m = size(gradCorY{1},1);

    M = {1,n_photon+1};
    for i = 1:n_photon+1
        Mi = zeros(m, m*(n_photon + 1));
        indr = [1:m]';
        indc = [i:n_photon+1:m*(n_photon + 1)]';
        indx = sub2ind(size(Mi),indr,indc);
        Mi(indx) = 1;
        M{i}=Mi;
    end

    sumYieldObjective = 0;
        
    %Construct sum of all yield components used in objective funstion
    for indexPhoton = 0:1:n_photon
        %Yields for n photons as vector
        vecYn = vec(M{indexPhoton+1}*Y);
        %gradient for yields of n photons as vector
        vecgradCorYn = vec(gradCorY{indexPhoton+1});

        %add gradYn*Yn to objective
        sumYieldObjective = sumYieldObjective + transpose(vecgradCorYn)*vecYn;
    end
    fval = sumYieldObjective;
end

%%%%%%%%%%%%%%%%%%%%% Helper function for FW iteration converting matrix in cell array to vector for searching %%%%%%%%%%%%%%%%%%%%%%%%%


function cellout = stackentriesCell(cellin)
    stackedvec = [reshape(real(cellin{1}),[],1); reshape(imag(cellin{1}),[],1); reshape(cellin{2},[],1); reshape(cellin{3},[],1)];
    cellout = stackedvec;
end