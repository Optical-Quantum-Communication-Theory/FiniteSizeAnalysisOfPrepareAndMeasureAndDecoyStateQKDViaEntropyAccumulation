function [lowerbnd_final,upperbnd,flag] = ImprovedSecondOrderOpt(crossover,stateGenRounds,stateTestRounds,acceptfrequency,observables,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob)
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
    
    %Initial guess using Choi Matrix
    J0 = eye(dimAprime*dimB)/dimB;
    Jinit= closestChoiMatrix(J0,stateTestRounds,acceptfrequency,observables,dimA,dimAprime,dimB);
    
    %%%%%%%%%%%%%%%%%%%%% FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%
    
    status = [];
    
    %Maximum iteratins
    maxiter = 50;
    
    %Tolerance in checking vertices are in active set
    tolvert = 1e-10;
    
    %Initialise initial guess x0
    Jmat = Jinit;
    
    %Coefficients of vectors in active set
    alpha_t = {1};
    
    %Indices indicating active vectors in S_t (Although S_t is the active
    %set, vectors are not removed only their indices in I_active)
    I_active = [1];
    
    %Initialize active set
    S_t = {Jmat};
    
    %search precision for line step
    linesearchprecision = 1e-20;
    
    %Tolerance in solution
    maxgap = 1e-6;
    
    %Perform iteration
    for t = 1:maxiter
        
        fprintf('FW iteration:%d',t)
        tstart_FW=tic;
    
        is_away = false;
        
        %Gradients
        gradW = (1-testprob)^2/log(2)*gradientRelEntChoi(Jmat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
        gradCor = gradientCorChoi(Jmat,crossover,stateTestRounds,observables,testprob,nu,dimA,dimAprime,dimB,dimR);
        
        %Find step direction for FW step
        %solve for FW step
        [sJt, status] = subproblem(gradW,gradCor,dimAprime,dimB,status); 

        %step direction in FW case
        dJt_FW = sJt - Jmat;
    
    
        %Find step direction for Away step
        %Evaluate <grad(f), v> for all v in S_t
        vtcell = cellfun(@(x) trace(gradW'*x) + trace(gradCor'*x),S_t,'Uniform', 0);
    
        %convert to array
        vtmat = cell2mat(vtcell);
    
        %maximise only over active vectors
        [~,id] = max(vtmat(I_active));
    
        %save index of maximum
        index_away = I_active(id);
        
        %vector reaching maximum of <grad(f), v> for all v in S_t 
        vJt = S_t{index_away};
    
        %step direction in Away step case
        dJt_A = Jmat - vJt;
    
        %Decide if FW or away step is taken
        if - trace(gradW'*dJt_FW) - trace(gradCor'*dJt_FW) >= - trace(gradW'*dJt_A) - trace(gradCor'*dJt_A)
            %save identifier for away step
            is_away = false;
    
            %Step direction is taken to be a FW step
            dJt = dJt_FW;
    
            %Maximum step size set to 1
            gamma_max = 1;
        else
            %save identifier for away step
            is_away = true;
    
            %Step direction is taken to be an Away step
            dJt = dJt_A;
    
            %Maximum step size set to alpha_t(index maximum in
            %S_t)/(1-alpha_t(index maximum in S_t))
            gamma_max = min(alpha_t{index_away}/(1-alpha_t{index_away}),1);
        end
    
        % perform an exact line search for step size
        optimoptions = optimset('TolX',linesearchprecision);
        minfunc = @(delta) primalf(Jmat + delta*dJt,crossover,acceptfrequency,observables,stateGenRounds,stateTestRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
            
        %step size given by solution
        gamma_t = fminbnd(minfunc,0,gamma_max,optimoptions);
        
    
        %Perform updates of active set S_t, alpha_t and set of active
        %indices I_active
    
        % Away step:
        if is_away == true
            %All Away steps:
            disp(" ")
            fprintf("Away step")
    
            %For all away steps all alpha_t are multiplied by (1+gamma_t)
            alpha_t = cellfun(@(x) (1+gamma_t)*x,alpha_t,'Uniform', 0);
                        
            if abs(gamma_t - gamma_max) < 10*eps(1) %Only Drop steps
                fprintf(" and Drop step")
    
                %If a drop step happens the coeficient alpha is set to 0
                alpha_t{index_away} = 0;
                %and the according vector is set inactive by removing it
                %from the set of active indices
                I_active(I_active == index_away) = [];
    
            else %no drop step
    
                %If it wasn't a drop step S_t stays the same and the alpha
                %corresponding to the maximum v in S_t is updated
                alpha_t{index_away} = alpha_t{index_away} - gamma_t;
            end
    
        %FW step:    
        else
            %For all FW steps all alpha_t are multiplied by (1-gamma_t)
            alpha_t = cellfun(@(x) (1-gamma_t)*x,alpha_t,'Uniform', 0);
            
            % Is the step st a new vector?
            %transform S_t into cell array of vectors
            S_t_cellvec = cellfun(@(x) stackentriesCell(x),S_t,'Uniform', 0);
            FWindex = find(ismembertol([S_t_cellvec{:}]',[reshape(real(sJt),[],1); reshape(imag(sJt),[],1)]',tolvert,'ByRows',true));
                                                         
            %If st is not in S_t
            if isempty(FWindex) 
                %add st to S_t
                S_t{end+1} = sJt;
    
                %add step size gamma_t as another component to alpha_t
                alpha_t{end+1} = gamma_t;
    
                %Add new vector as an active vector to the set of active
                %vectors (will be last entry in S_t, thus numel(S_t) added
                %to I_active)
                I_active = [I_active, numel(S_t)];
            
            %If st is in S_t
            else 
                if alpha_t{FWindex} < eps %vector was in St but not active
                    %Activate index corresponding to vector st
                    I_active = [I_active, FWindex];
                end
                %vector was in S_t and the coefficient corresponding to st is updated
                %accordingly by only adding gamma_t since the
                %multiplication was done before
                alpha_t(FWindex) = cellfun(@(x) x + gamma_t, alpha_t(FWindex),'Uniform', 0);
            end
            
            %Exceptional case if gamma is 1, the set S_t collapses to 
            % S_t = {st}, here done by removing all indices from set of
            % active indices and replacing them by index of st
            if gamma_t > 1-10*eps(1) 
                I_active = [FWindex];
            end
        end
        
        %Calculate next step of optimisation variable
        Jnew = Jmat + gamma_t*dJt;
        
        %Calculate duality gap
        gap = - trace(gradW'*dJt_FW) - trace(gradCor'*dJt_FW);
    
        %Calculate projection value
        proj_val = trace(gradW'*Jnew) + trace(gradCor'*Jnew);
    
        %Calculate function value of new and old optimisation variable 
        fNew = primalf(Jnew,crossover,acceptfrequency,observables,stateGenRounds,stateTestRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
        fOld = primalf(Jmat,crossover,acceptfrequency,observables,stateGenRounds,stateTestRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
    
        %Calculate lower bound on minimum
        lowerbnd = fOld + trace(gradW'*dJt_FW) + trace(gradCor'*dJt_FW);
    
        %print progress
        t_FW=toc(tstart_FW);
        fprintf('\n    FW iteration time: %f\n',t_FW);
        fprintf('lower bound:%f     projection value:%f    gap:%f    fvalue:%f    fgap:%f\n',lowerbnd,proj_val,gap,fNew,fOld-fNew)
        
        %Stopping criterion if gap is small enough
        if  ( abs(gap) < maxgap)
            Jmat = Jmat + gamma_t*dJt;
            break;
        end
        
        %Update optimisation variable
        Jmat = Jmat + gamma_t*dJt;
    
        %Display warning if algorithm hasn't converged and reached the
        %maximum iterations
        if t == maxiter
            disp('**** Warning: Maximum iteration limit reached. ****');
        end
    end

    %Fix Choi matrix to be hermitian
    J_opt = (Jmat +Jmat')/2;

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
    gradW_eps = (1-testprob)^2/log(2)*gradientRelEntChoi_eps(perturbation,J_opt,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
    gradCor = gradientCorChoi(J_opt,crossover,stateTestRounds,observables,testprob,nu,dimA,dimAprime,dimB,dimR);
    
    %Run one more FW iteration at optimal Choi matrix and perturbation for
    %gradient dJt_final
    statusfinal = [];
    [sJt_final, statusfinal] = subproblem(gradW,gradCor,dimAprime,dimB,status);

    dJt_final = sJt_final - J_opt;
    
    %Function value at correctly perturbed optimal points
    f_final = primalf_eps(perturbation,Jmat,crossover,acceptfrequency,observables,stateGenRounds,stateTestRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
    
    %Calculate correction term due to perturbation
    if perturbation == 0
        zetaEp = 0;
    else
        zetaEp = 2*perturbation*(dprime-1)*log(dprime/(perturbation*(dprime-1)));
    end
    
    upperbnd = f_final;
    lowerbnd_final = f_final + trace(gradW_eps'*dJt_final) + trace(gradCor'*dJt_final) - zetaEp;
    flag = statusfinal;
end

%%%%%%%%%%%%%%%%%%%%% Objective functions of FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

%Objective function in FW iteration with perturbation calculated at each
%iteration

function fval = primalf(ChoiMat,crossover,acceptfrequency,observables,stateGenRounds,stateTestRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob)
    rhoTest = PartialTrace(kron(eye(dimA),ChoiMat)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
    
    m = size(crossover,1);
    n = size(crossover,2);

    genfrequency = zeros(m,n);

    for indexrow = 1:m
        for indexcolumn = 1:n
            genfrequency(indexrow,indexcolumn) = trace(observables{indexrow,indexcolumn}'*rhoTest);
        end
    end
    
    fval = (1-testprob)^2/log(2)*relEntropyChoi(ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB) + sum(crossover.*(acceptfrequency-genfrequency),'all') - nu/(1-nu)*log(2)/2*Varq(genfrequency,crossover,testprob,dimR);
end

%Objective function with fixed perturbation for final calculation

function fval = primalf_eps(perturbation,ChoiMat,crossover,acceptfrequency,observables,stateGenRounds,stateTestRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob)
    rhoTest = PartialTrace(kron(eye(dimA),ChoiMat)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
    
    m = size(crossover,1);
    n = size(crossover,2);

    genfrequency = zeros(m,n);

    for indexrow = 1:m
        for indexcolumn = 1:n
            genfrequency(indexrow,indexcolumn) = trace(observables{indexrow,indexcolumn}'*rhoTest);
        end
    end
    
    fval = (1-testprob)^2/log(2)*relEntropyChoi_eps(perturbation,ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB) + sum(crossover.*(acceptfrequency-genfrequency),'all') - nu/(1-nu)*log(2)/2*Varq(genfrequency,crossover,testprob,dimR);
end


function fval = Varq(genfrequency,crossover,testprob,dimR)
    varp = 1/testprob*sum(genfrequency.*(max(crossover, [],'all') - crossover).^2,'all') - (max(crossover, [],'all') - sum(crossover.*genfrequency,'all'))^2;
    fval = (log2(1+2*dimR) + sqrt(2 + varp))^2;
end

%%%%%%%%%%%%%%%%%%%%% Subproblem calculating deltaJ %%%%%%%%%%%%%%%%%%%%%%%%%

function [Jtilde,status] = subproblem(gradW,gradCor,dimAprime,dimB,status)
    %Dimension of Choi matrix
    dim = dimAprime*dimB;
    
    cvx_solver Mosek
    cvx_precision high

    cvx_begin sdp quiet
        variable J(dim,dim) hermitian semidefinite
        minimize real(trace(gradW'*J)) + real(trace(gradCor'*J))
        subject to
                        
            %Partial trace constraint on J
            PartialTrace(J,[2],[dimAprime,dimB]) - eye(dimAprime) == zeros(dimAprime);

    cvx_end
    if strcmp(cvx_status, 'Infeasible') % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, subproblem status: %s ****\n",cvx_status);
    end
    
    %record status for debugging
    status = [status, string(cvx_status)];
    
    Jtilde = full(J);
end

%%%%%%%%%%%%%%%%%%%%% Closest Choi matrix %%%%%%%%%%%%%%%%%%%%%%%%%

function Choi = closestChoiMatrix(J0,stateTestRounds,expectations,observables,dimA,dimAprime,dimB)
    %solving for closest Choi matrix compatible with our observations

    %disp("solving for closest Choi matrix")

    dim = dimAprime*dimB;

    cvx_solver mosek
    cvx_precision high
    
    cvx_begin sdp quiet
        variable J(dim,dim) hermitian semidefinite
        minimize norm(J0-J)
        subject to
            
            %Partial trace constraint Tr_B[J] = 1_A'
            PartialTrace(J,[2],[dimAprime,dimB]) - eye(dimAprime) == zeros(dimAprime);

            %Calculate rho in test rounds
            rho = PartialTrace(kron(eye(dimA),J)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
            for i=1:numel(observables)
                trace(observables{i}'*rho) - expectations(i) == 0;
            end
    cvx_end
    Choi = (J+J')/2;
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

%%%%%%%%%%%%%%%%%%%%% Helper function for FW iteration converting matrix in cell array to vector for searching %%%%%%%%%%%%%%%%%%%%%%%%%


function cellout = stackentriesCell(cellin)
    stackedvec = [reshape(real(cellin),[],1);reshape(imag(cellin),[],1)];
    cellout = stackedvec;
end