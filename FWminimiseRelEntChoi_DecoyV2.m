function [optvalue,crossover,optChoi,rhoGen,rhoTest] = FWminimiseRelEntChoi_DecoyV2(stateGenRounds,stateTestRounds,acceptfrequency,observables,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob,decoys,probsDecoy,n_photon)
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
%       decoys:             decoy intensities entered as a list [mu_signal, mu_2, mu_3,...]
%       probsDecoy:         probabilities for decoys
%       n_photon:           photon number cut-off
%
%Output: optvalue:          optimal value of optimisation
%        crossover:         gradient of the crossover min trade-off function
%        optChoi:           optimal Choi matrix
%        rhoGen:            resulting rho in generation rounds after optimal attack with optChoi
%        rhoTest:           resulting rho in test rounds after optimal attack with optChoi

    
    %%%%%%%%%%%%%%%%%%%%% Parameters for FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%
    %Enable printing
    verbose = 1;

    %Store status of FW steps
    status = [];
       
    %Maximum iteratins
    maxiter = 300;
    
    %Tolerance in checking vertices are in active set
    tolvert = 1e-8;
    
    %search precision for line step
    linesearchprecision = 1e-20;
    
    %Tolerance in solution
    maxgap = 1e-5;

    %%%%%%%%%%%%%%%%%%%%% Initial guess %%%%%%%%%%%%%%%%%%%%%%%%%

    %Initial guess using Choi Matrix
    J0 = eye(dimAprime*dimB)/dimAprime;

    [Jinit, xplusinit, lambdainit, Yinit] = closestChoiMatrixDecoy(J0,stateTestRounds,observables,acceptfrequency,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon);
        
    %%%%%%%%%%%%%%%%%%%%% Prepare FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

    %Initialise initial guess for FW iteration
    Jmat = Jinit;
    xplus = xplusinit;
    
    %Coefficients of vectors in active set
    alpha_t = {1};
    
    %Indices indicating active vectors in S_t (Although S_t is the active
    %set, vectors are not removed only their indices in I_active)
    I_active = [1];
    
    %Initialize active set
    S_t = {{Jmat,xplus}};

    %Keep track of gammas
    gamma_t_list = [];
    
    %%%%%%%%%%%%%%%%%%%% Run FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Perform iterations
    if verbose == 1
        fprintf("Performing FW Iterations \n")
    end

    for t = 1:maxiter
    
        if verbose == 1        
            fprintf('FW iteration:%d',t)
        end
        
        tstart_FW=tic;
   
        
        %Gradients
        gradW = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*gradientRelEntChoi(Jmat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
        grads = gradientdualS(xplus,dimR,nu,testprob);
        
        %Find step direction for FW step
        %solve for FW step
        [sJt, stxplus, Yield, lambda,status] = subproblem(stateTestRounds,observables,acceptfrequency,gradW,grads,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon,status);
        
        %step direction in FW case
        dJt_FW = sJt - Jmat;
        dxt_FW = stxplus - xplus;
    
        %Find step direction for Away step
        %Evaluate <grad(f), v> for all v in S_t
        vtcell = cellfun(@(x) trace(gradW'*x{1}) + grads'*x{2},S_t,'Uniform', 0);
    
        %convert to array
        vtmat = cell2mat(vtcell);
    
        %maximise only over active vectors
        [~,id] = max(vtmat(I_active));
    
        %save index of maximum
        index_away = I_active(id);
        
        %vector reaching maximum of <grad(f), v> for all v in S_t 
        vJt = S_t{index_away}{1};
        vxplust = S_t{index_away}{2};
    
        %step direction in Away step case
        dJt_A = Jmat - vJt;
        dxt_A = xplus - vxplust;
    
        %Pair direction
        dJt = sJt - vJt;
        dxt = stxplus - vxplust;
        gamma_max = alpha_t{index_away};
        
    
        % perform an exact line search for step size
        optimoptions = optimset('TolX',linesearchprecision,'Tolfun',linesearchprecision);
        minfunc = @(delta) primalf(Jmat + delta*dJt, xplus + delta*dxt ,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,decoys,nu,testprob);
        %step size given by solution
        gamma_t = fminbnd(minfunc,0,gamma_max,optimoptions);
        
        %Perform updates of active set S_t, alpha_t and set of active
        %indices I_active
    
        % Away part:
        if ismembertol(gamma_max, gamma_t,1e-4) == 1 %|| gamma_t <= 1e-12
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
        FWindex = find(ismembertol([S_t_cellvec{:}]',[reshape(real(sJt),[],1); reshape(imag(sJt),[],1) ;stxplus]',tolvert,'ByRows',true,'DataScale',1));
    
        %If st is not in S_t
        if isempty(FWindex) 
            %add st to S_t
            S_t{end+1} = {sJt, stxplus};
    
            %add step size gamma_t as another component to alpha_t
            alpha_t{end+1} = gamma_t;
    
            %Add new vector as an active vector to the set of active
            %vectors (will be last entry in S_t, thus numel(S_t) added
            %to I_active)
            I_active = [I_active, numel(S_t)];
        
        %If st is in S_t
        else
            alphaMatFW = [alpha_t{FWindex}];
            
            if isempty(alphaMatFW(abs(alphaMatFW) < 1e-13)) ~= 1 %vector was in St but not active
                %Activate index corresponding to vector st (.' for case that
                %two vales are in FWindex)
                I_active = [I_active, FWindex.'];
            end
    
            %vector was in S_t and the coefficient corresponding to st is updated
            %accordingly by only adding gamma_t since the
            %multiplication was done before
            alpha_t(FWindex) = cellfun(@(x) x + gamma_t ,alpha_t(FWindex),'Uniform', 0); 
        end
        
        %Exceptional case if gamma is 1, the set S_t collapses to 
        % S_t = {st}, here done by removing all indices from set of
        % active indices and replacing them by index of st
        if ismembertol(gamma_t, 1, 1e-13) == 1
            I_active = [FWindex];
        end
        
        
        %Calculate next step of optimisation variable
        Jnew = Jmat + gamma_t*dJt;
        xplusnew = xplus + gamma_t*dxt;
        
        %Calculate gap
        gap = - trace(gradW'*dJt_FW) - grads'*dxt_FW;
    
        %Calculate projection value
        proj_val = trace(gradW'*Jnew) + grads'*xplusnew;
    
        %Calculate function value of new and old optimisation variable 
        fnew = primalf(Jnew, xplusnew,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,decoys,nu,testprob);
        fold = primalf(Jmat,xplus,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,decoys,nu,testprob);
    
        %Calculate lower bound on minimum
        lowerbnd = fold + trace(gradW'*dJt_FW) + grads'*dxt_FW;
        
        %time needed for FW step
        t_FW=toc(tstart_FW);
    
        %print progress
        if verbose == 1        
            fprintf('\n    FW iteration time: %f\n',t_FW);
            fprintf('lower bound:%f     projection value:%f    gap:%f    fvalue:%f    fgap:%f\n',lowerbnd,proj_val,gap,fnew,fnew-fold)
        end
    
        %Stopping criterion if gap is small enough
        if  abs(gap) <= maxgap %abs(gap/lowerbnd) <= maxgap || abs(gap) <= 1e-6
            Jmat = Jmat + gamma_t*dJt;
            xplus = xplus + gamma_t*dxt;
            break;
        end
        
        %Update optimisation variable
        Jmat = Jmat + gamma_t*dJt;
        xplus = xplus + gamma_t*dxt;
        
        %Display warning if algorithm hasn't converged and reached the
        %maximum iterations
        if t == maxiter
            disp('**** Warning: Maximum iteration limit reached. ****');
        end
    end
    
    %Assign optimal values
    optChoi = Jmat;
    optxplus = xplus;
    optvalue = primalf(Jmat,xplus,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,decoys,nu,testprob);
    
    %Assign optimal crossover min trade-off function with same shape as the
    %accept frequency
    gradW = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*gradientRelEntChoi(optChoi,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
    grads = gradientdualS(optxplus,dimR,nu,testprob);
    crossover = getMinTradeoffDecoy(stateTestRounds,observables,acceptfrequency,gradW,grads,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon);

    %Calculate rho in test rounds
    rhoTest = PartialTrace(kron(eye(dimA),optChoi)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);

    %Calculate rho in generation rounds
    rhoGen = PartialTrace(kron(eye(dimA),optChoi)*(kron(PartialTranspose(stateGenRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
end

%%%%%%%%%%%%%%%%%%%%% Objective function of FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

function fval = primalf(ChoiMat,tplus,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,decoys,nu,testprob)
    fval = Poisson(decoys(1),1)*(1-testprob)^2/log(2)*relEntropyChoi(ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB) + dualS(tplus,dimR,nu,testprob);
end

%%%%%%%%%%%%%%%%%%%%% Closest Choi matrix %%%%%%%%%%%%%%%%%%%%%%%%%
function [Jtilde,tplustilde, lambdatilde,Yield] = closestChoiMatrixDecoy(J0,stateTestRounds,observables,acceptfrequency,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon)
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
        variable xplus(nobs*n_decoy) nonnegative
        variable lambda(m,n,n_decoy)
        variable YieldCutOff(m,n,n_decoy)

        minimize norm(J - J0) + norm(xplus)
        subject to

            %Density matrix in test rounds given Alice sent a single photon
            rhoTest = PartialTrace(kron(eye(dimA),J)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
            
            % 0 <= Y 
            vec(Y) >= vec(zeros(m*(n_photon+1),n));
            vec(Y) <= vec(ones(m*(n_photon+1),n));

            % 1-photon component treated seperately with add.
            %constraints
            Y1 = M{2}*Y;

            %0-photon yields
            Y0 = M{1}*Y;

            %0-photon error rate e_0=1/2
            Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));

            % Require that yields sum to prob of Alice sending the signal
            for indexPhoton = 0: n_photon
                Yi = M{indexPhoton+1}*Y;
                sum(Yi,2) == probsAliceTest;
            end

            %sum(lammda) = 0
            ones(1,m*n*n_decoy)*vec(lambda) == 0; %decoy_tolerance;
            
            %Usual decoy bounds rewritten as matrix times vector
            
            for kdecoy = 1:n_decoy
                %Initialize probabilities
                pmu = Poisson(decoys(kdecoy),0:1:n_photon);
                Pmu = kron(eye(m),pmu);

                %Compute sum_n p_mu(n) Y_n
                PmuY = Pmu*Y;
                
                %Contraints on delta^mu
                % 0<= delta^mu <= 1 - p_tot(mu)
                vec(YieldCutOff(:,:,kdecoy)) >= zeros(m*n,1);
                vec(YieldCutOff(:,:,kdecoy)) <= (1-sum(pmu))*ones(m*n,1) + decoy_tolerance;
                
                %Remaining constraints on yields from observations
                for indexrow = 1:m
                    for indexcol = 1:n
                        Y1(indexrow,indexcol) - trace(observables{indexrow,indexcol}'*rhoTest) == 0;% decoy_tolerance;
                        0 == acceptfrequency(indexrow,indexcol,kdecoy) - probsDecoy(kdecoy)*(PmuY(indexrow,indexcol) + YieldCutOff(indexrow,indexcol,kdecoy)) - lambda(indexrow,indexcol,kdecoy);
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
    
    
    Jtilde = (full(J)+full(J)')/2;
    tplustilde = xplus;
    lambdatilde = lambda;
    Yield = full(Y);
end

%%%%%%%%%%%%%%%%%%%%% Subproblem calculating deltaJ %%%%%%%%%%%%%%%%%%%%%%%%%
function [Jtilde,tplustilde, Yield,lambdatilde,status] = subproblem(stateTestRounds,observables,acceptfrequency,gradW,grads,dimA,dimAprime,dimB,probsDecoy,decoys,n_photon,status)
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
                        0 == acceptfrequency(indexrow,indexcol,kdecoy)/probsDecoy(kdecoy) - (PmuY(indexrow,indexcol) + YieldCutOff(indexrow,indexcol,kdecoy)) - lambda(indexrow,indexcol,kdecoy)/probsDecoy(kdecoy);
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
    
    %record status for debugging
    status = [status, string(cvx_status)];
    
    Jtilde = (full(J)+full(J)')/2;
    tplustilde = xplus;
    lambdatilde = lambda;
    Yield = full(Y);
end

%%%%%%%%%%%%%%%%%%%%% Helper function for FW iteration converting matrix in cell array to vector for searching %%%%%%%%%%%%%%%%%%%%%%%%%

function cellout = stackentriesCell(cellin)
    stackedvec = [reshape(real(cellin{1}),[],1);reshape(imag(cellin{1}),[],1);reshape(cellin{2},[],1)];
    cellout = stackedvec;
end

%%%%%%%%%%%%%%%%%%%%% Possonian distribution %%%%%%%%%%%%%%%%%%%%%%%%%
function prob = Poisson(mu,n)
    prob = exp(-mu).*mu.^n./factorial(n);
end