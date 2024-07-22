function [optvalue,crossover,optChoi,rhoGen,rhoTest] = FWminimiseRelEntChoi(stateGenRounds,stateTestRounds,acceptfrequency,observables,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob)
%Function minimises relative entropy with penalty term T* and finds the
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
%Output: optvalue:          optimal value of optimisation
%        crossover:         gradient of the crossover min trade-off function
%        optChoi:           optimal Choi matrix
%        rhoGen:            resulting rho in generation rounds after optimal attack with optChoi
%        rhoTest:           resulting rho in test rounds after optimal attack with optChoi

    %%%%%%%%%%%%%%%%%%%%% Parameters for FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

    status = [];
    
    %Maximum iteratins
    maxiter = 80;
    
    %Tolerance in checking vertices are in active set
    tolvert = 1e-10;

    %Tolerance in solution
    maxgap = 1e-6;

    %search precision for line step
    linesearchprecision = 1e-20;

    %%%%%%%%%%%%%%%%%%%%% Initial guess %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Initial guess using Choi Matrix
    J0 = eye(dimAprime*dimB)/dimB;
    Jinit= closestChoiMatrix(J0,stateTestRounds,acceptfrequency,observables,dimA,dimAprime,dimB);
    
    %Calculate corresponding density matrix to check validity
    rhoInitChoi = PartialTrace(kron(eye(dimA),Jinit)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
    
    %%%%%%%%%%%%%%%%%%%%% FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Initialise initial guess x0
    Jmat = Jinit;
    xplus = zeros(numel(observables),1);
    
    %Coefficients of vectors in active set
    alpha_t = {1};
    
    %Indices indicating active vectors in S_t (Although S_t is the active
    %set, vectors are not removed only their indices in I_active)
    I_active = [1];
    
    %Initialize active set
    S_t = {{Jmat,xplus}};
   
        
    %Perform iteration
    for t = 1:maxiter
        
        fprintf('FW iteration:%d',t)
        tstart_FW=tic;
        
        %Gradients
        gradW = (1-testprob)^2/log(2)*gradientRelEntChoi(Jmat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB);
        grads = gradientdualS(xplus,dimR,nu,testprob);
        
        %Find step direction for FW step
        %solve for FW step
        [sJt, stxplus, lambda,duallambda,status] = subproblem(stateTestRounds,observables,acceptfrequency,gradW,grads,dimA,dimAprime,dimB,testprob,status);
        
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
        optimoptions = optimset('TolX',linesearchprecision);
        minfunc = @(delta) primalf(Jmat + delta*dJt, xplus + delta*dxt ,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
        %step size given by solution
        gamma_t = fminbnd(minfunc,0,gamma_max,optimoptions);
        
        %Perform updates of active set S_t, alpha_t and set of active
        %indices I_active
    
        %Away part:
    
        if abs(gamma_t - gamma_max) < 1e-13
            fprintf("Drop step")
    
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
            
            if isempty(alphaMatFW(alphaMatFW < 1e-13)) ~= 1 %vector was in St but not active
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
        if gamma_t > 1- (1e-13)
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
        fnew = primalf(Jnew, xplusnew,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
        fold = primalf(Jmat,xplus,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
    
        %Calculate lower bound on minimum
        lowerbnd = fold + trace(gradW'*dJt_FW) + grads'*dxt_FW;
    
        %print progress
        t_FW=toc(tstart_FW);
        fprintf('\n    FW iteration time: %f\n',t_FW);
        fprintf('lower bound:%f     projection value:%f    gap:%f    fvalue:%f    fgap:%f\n',lowerbnd,proj_val,gap,fnew,fold-fnew)
        
        %Stopping criterion if gap is small enough
        if  ( abs(gap) < maxgap)
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
    optvalue = primalf(Jmat,xplus,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob);
    
    %Assign optimal crossover min trade-off function with same shape as the
    %accept frequency
    crossover = cell2mat(reshape(duallambda,size(acceptfrequency)));

    %Calculate rho in test rounds
    rhoTest = PartialTrace(kron(eye(dimA),optChoi)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);

    %Calculate rho in generation rounds
    rhoGen = PartialTrace(kron(eye(dimA),optChoi)*(kron(PartialTranspose(stateGenRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
end

%%%%%%%%%%%%%%%%%%%%% Objective function of FW iteration %%%%%%%%%%%%%%%%%%%%%%%%%

function fval = primalf(ChoiMat,tplus,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB,dimR,nu,testprob)
    fval = (1-testprob)^2/log(2)*relEntropyChoi(ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB) + dualS(tplus,dimR,nu,testprob);
end

%%%%%%%%%%%%%%%%%%%%% Subproblem calculating deltaJ %%%%%%%%%%%%%%%%%%%%%%%%%
function [Jtilde,tplustilde,lambdatilde,duallambdatilde,status] = subproblem(stateTestRounds,observables,acceptfrequency,gradW,grads,dimA,dimAprime,dimB,testprob,status)
    %Dimension of Choi matrix
    dim = dimAprime*dimB;

    %number of observables
    nobs = numel(observables);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(acceptfrequency,1);
    n = size(acceptfrequency,2);
    
    %Select solver
    cvx_solver Mosek
    cvx_precision high
    
    %Setup SDP
    cvx_begin sdp quiet
        variable J(dim,dim) hermitian semidefinite
        variable tplus(nobs) nonnegative
        variables lambda(m*n)
        dual variables duallambda{m*n}
        minimize real(trace(gradW'*J)) + transpose(grads)*tplus
        subject to
            
            %Calculate rho in test rounds
            rho = PartialTrace(kron(eye(dimA),J)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
            
            %Lambda is equal to difference between accepted freq and Tr[Gamma rho] 
            for i = 1:numel(observables)
                duallambda{i} : 0 == acceptfrequency(i) - trace(observables{i}'*rho) - lambda(i) ; %assigns dual variable to constraints
            end

            %xplus >= lambda
            vec(tplus) >= abs(vec(lambda));

            %Partial trace constraint on J
            PartialTrace(J,[2],[dimAprime,dimB]) - eye(dimAprime) == zeros(dimAprime);

    cvx_end
    if strcmp(cvx_status, 'Infeasible') % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, subproblem status: %s ****\n",cvx_status);
    end
    
    %record status for debugging
    status = [status, string(cvx_status)];
    
    Jtilde = full(J);
    tplustilde = tplus;
    lambdatilde = lambda;
    duallambdatilde = duallambda;
end

%%%%%%%%%%%%%%%%%%%%% Closest Choi matrix %%%%%%%%%%%%%%%%%%%%%%%%%

function Choi = closestChoiMatrix(J0,stateTestRounds,expectations,observables,dimA,dimAprime,dimB)
    %solving for closest Choi matrix compatible with our observations

    %disp("solving for closest Choi matrix")
    
    %Dimension of Choi matrix
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
    Choi = J;
end

%%%%%%%%%%%%%%%%%%%%% Helper function for FW iteration converting matrix in cell array to vector for searching %%%%%%%%%%%%%%%%%%%%%%%%%

function cellout = stackentriesCell(cellin)
    stackedvec = [reshape(real(cellin{1}),[],1);reshape(imag(cellin{1}),[],1);reshape(cellin{2},[],1)];
    cellout = stackedvec;
end