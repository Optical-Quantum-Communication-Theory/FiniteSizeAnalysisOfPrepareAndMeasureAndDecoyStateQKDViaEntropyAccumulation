function [iqw,optdeltacom , uppbndstest , lowbndstest , uppbndgen , lowbndgen , pabortATbnd] = optdeltacom(N,testprob,acceptfrequency,crossover,epscomAT,numtrials)
    % Computes a value of the "completeness keyrate penalty" Delta_com for the specified parameters, using fmincon to try optimizing the acceptance set to minimize Delta_com
    % Requires that acceptfrequency and crossover must be arrays of the same dimension, but otherwise their dimensions can be arbitrary 
    % N: number of rounds
    % testprob: test probability
    % acceptfrequency: expected honest (normalized) distribution on test-outcome alphabet
    % crossover: gradient of crossover fmin
    % epscomAT: desired completeness parameter for acceptance test
    % numtrials: number of fmincon trials to run (if set to 0, will just use a specific valid choice instead of optimizing)

    % Output values are defined as follows:
    % optdeltacom: best value of Delta_com found
    % uppbndstest: corresponding choice of upper-bound t values for the test-round terms
    % lowbndstest: corresponding choice of lower-bound t values for the test-round terms (returned as negative values)
    % uppbndgen: corresponding choice of upper-bound t value for the generation-round term
    % lowbndgen: corresponding choice of lower-bound t value for the generation-round term (returned as negative value)
    % pabortATbnd: actual computed bound on P[AT aborts]; should be smaller than epscomAT
    % iqw: a flag when the LP is infeasbile

    % Some basic input validation tests
    if size(acceptfrequency) ~= size(crossover)
        error('Inconsistent sizes of acceptfrequency and crossover')
    elseif abs(sum(acceptfrequency(:))-1) > 1e-10 || min(sum(acceptfrequency(:))) < 0
        sum(acceptfrequency(:))
        error('acceptfrequency is not a probability distribution')
    elseif numtrials < 0
        error('Invalid numtrials value')
    end
    
    % Compute gradient of actual (not crossover) min-tradeoff function
    gradientf = [crossover(:)', max(crossover(:))]/testprob;
    % Compute expected honest distribution on full "parameter-estimation register"
    acceptfrequencyfull = [testprob*acceptfrequency(:)', 1-testprob];
    % Dimension of "parameter-estimation register", which we shall call C
    dimC = length(gradientf);
    % Tolerance by which to relax the LP constraints (this slightly loosens the bound on Delta_com, but in the reliable direction)
    extratol = 1e-11;

    % We now set up inputs to fmincon for finding the choice of t bounds that minimizes deltacom
    % We use the convention that all t values are positive, and only negate the lower-bound values when calling entrywiseLP
    % First we define the objective function, based on computing a closed-form solution to the completeness LP
    % (We'll use the helper_minus function to handle fmincon violating non-negativity constraints; this will be accounted for when constructing best solution later)
    deltacom = @(tbnds)entrywiseLPtbnds(-gradientf,helper_minus(tbnds),extratol);
    % No linear inequality/equality constraints needed, apart from simple upper/lower bounds 
    % Note we will not constrain t values to be such that the bounds on p_acc are within [0,1]; such a constraint creates problems in fmincon when there are small entries in acceptfrequencyfull. Instead we will rely on the P_single_fail function being appropriately defined in such regimes.
    Aineq = []; bineq = []; Aeq = []; beq = []; lb = 1e-14*ones(1,2*dimC); ub = ones(1,2*dimC);
    % Define the nonlinear constraint P[abort] < epscomAT (with a rescaling parameter to increase "penalty" for constraint violation)
    nonlcon = @(tbnds)epscomAT_constraint(N,acceptfrequencyfull,helper_minus(tbnds),epscomAT,1e5);
    % Set fmincon options
    options = optimset('Display','off'); 

    % Prepare variables to store results from multiple trials of optimizing tbnds
    trialsols = zeros(numtrials+1,2*dimC); trialvals = zeros(numtrials+1,1); 
    % First note that a simple first choice for tbnds would be one where the upper/lower t values are symmetric and such that epscomAT is "uniformly distributed" across contributions
    % Let us compute that case and use it as our first viable solution
    simpletbnds = simple_t_from_eps(N,acceptfrequencyfull,epscomAT);
    trialsols(1,:) = simpletbnds; trialvals(1) = deltacom(simpletbnds);
    % Now let's run fmincon for specified number of trials to gather other possible solutions
    % Note that the initial points chosen in these trials (other than the first) are not necessarily initially feasible wrt epsAT constraint
    for j=2:numtrials+1
        if j==2 % On the first trial, use simpletbnds as initial point
            inipt = simpletbnds;
        elseif mod(j,2)==0 % Otherwise, on "half" the remaining trials use a random perturbation of that above choice 
            randomscales = .9 + 0.2*rand(1,2*dimC);
            inipt = randomscales.*simpletbnds;
        else % On all other trials, just use random rescaling of acceptfrequencyfull 
            inipt = rand(1,2*dimC).*[acceptfrequencyfull , acceptfrequencyfull];
        end
        [trialsols(j,:) , trialvals(j)] = fmincon(deltacom,inipt,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
    end

    % We want the best solution out of the above trials, but since fmincon does not always obey constraints, we'll have to iterate until an actually valid solution is found
    for attempt=1:numtrials+1

        % Identify the current best claimed solution
        [bestval, bestvalpos] = min(trialvals);
        besttbnds = helper_minus(trialsols(bestvalpos,:)); % Using helper_minus preserves consistency with earlier inputs to fmincon; this means that in fact our final result besttbnds might not technically not be equal to any of the trialsols, but we'll verify it's valid
        % Compute its actual P[abort] bound by slightly abusing epscomAT_constraint function
        pabortATbnd = epscomAT_constraint(N,acceptfrequencyfull,besttbnds,0,1);
    
        if pabortATbnd < (1+1e-8)*epscomAT % If specified completeness bound was indeed satisfied (to relative error 1e-8), exit loop with this as best solution
            break
        else % Otherwise remove this claimed solution and proceed with loop to find next-best attempt
            trialsols(bestvalpos,:) = []; trialvals(bestvalpos) = [];
        end

        if attempt==numtrials+1
            error('No feasible solutions found!')
          %  iqw = 1

        end

    end
    %%%
    iqw = 0;
    % Sanity check: run the LP once using CVX instead, to check against the analytic solution computed above
    [certsol , certval] = entrywiseLPcert(-gradientf,[besttbnds(1:dimC);-besttbnds(dimC+1:2*dimC)],extratol);
    if certval == -Inf
        format long
        [besttbnds(1:dimC) ; -besttbnds(dimC+1:2*dimC)]'
        format short
        disp('CVX finds completeness LP infeasible, check above values')
        iqw = 1;
    elseif abs(bestval-certval)>1e-5 || norm(sum(certsol))>1e-8
        format long
        [bestval , certval , bestval-certval]
        [-gradientf ; besttbnds(1:dimC) ; -besttbnds(dimC+1:2*dimC); certsol]'
        sum(certsol)
        format short
        disp('Incompatibility between analytic/certified solution, check above values')
    end
    
    % Set final output values
    optdeltacom = bestval;
    uppbnds = besttbnds(1:dimC) + extratol;
    uppbndstest = reshape(uppbnds(1:dimC-1),size(acceptfrequency)); uppbndgen = uppbnds(dimC);
    lowbnds = -besttbnds(dimC+1:2*dimC) - extratol;
    lowbndstest = reshape(lowbnds(1:dimC-1),size(acceptfrequency)); lowbndgen = lowbnds(dimC);
    
end

function [constraintviolation , dummy] = epscomAT_constraint(N,acceptfrequencyfull,tbnds,epscomAT,scalefactor)
    % Imposes the nonlinear constraint P[abort] < epscomAT in required format for fmincon 
    % scalefactor is a factor by which to scale the constraint violation magnitude, to increase the penalty in fmincon for constraint violations
    % constraintviolation = scalefactor*(pabortATbnd - epscomAT) where pabortATbnd is the upper bound on AT aborting, i.e. we want this to be <=0 in the constrained optimization
    % dummy is simply a dummy variable to fulfill fmincon requirements on the nonlinear constraint

    dimC = length(acceptfrequencyfull);
    pabortATbnd = 0;
    for j = 1:dimC
        pabortATbnd = pabortATbnd + P_single_fail(N,acceptfrequencyfull(j),tbnds(j),tbnds(dimC+j));
    end
    constraintviolation = scalefactor*(pabortATbnd - epscomAT);
    dummy = 0;

end

function simpletbnds = simple_t_from_eps(N,acceptfrequencyfull,epscomAT)
    % Computes a simple choice of t values for entrywise accept conditions, in the following sense:
    % - Each entry's condition is a "symmetric" interval about acceptfrequencyfull(j)
    % - Each entry's condition has the same probability of failing in AT (apart from special-case handling for probabilities close to zero)

    dimC = length(acceptfrequencyfull);
    tol = 1e-14; % Tolerance-to-0 parameter: probabilities smaller than this will be treated separately
    num_nonzeros = sum(acceptfrequencyfull>tol); % Number of nonzero values (up to tolerance)

    % Compute desired t values (in terms of the magnitudes only, since the intervals are "symmetric")
    t_magnitudes = zeros(1,dimC);
    for j = 1:dimC
        if acceptfrequencyfull(j) < tol 
            % For probabilities < tol, we'll assume that setting the t value equal to tol is sufficient to make the P[abort] contribution from that term effectively zero
            t_magnitudes(j) = tol; 
        else
            % For other probabilities, compute t value such that the probability of that condition failing is 0.95*epscomAT/num_nonzeros 
            % (The 0.95 prefactor is to cope with numerical imprecisions that might otherwise cause the probabilities to sum to more than epscomAT)
            eps_from_t = @(t)(0.95*epscomAT/num_nonzeros - P_single_fail(N,acceptfrequencyfull(j),helper_minus(t),helper_minus(t))); % Using helper_minus to avoid problems with negative t; this should be a monotone increasing function of t
            if eps_from_t(tol) > 0 % Special-case handling if distribution is extreme enough that even t=tol yields Pfail < 0.95*epscomAT/num_nonzeros for that term
                t_magnitudes(j) = tol; 
            else % Assume that if setting t=tol does not make Pfail sufficiently small, then fzero can find some value (necessarily > tol) to make Pfail sufficiently small
                t_magnitudes(j) = helper_minus(fzero(eps_from_t,1e-3*acceptfrequencyfull(j))); % Using helper_minus consistently with above formulation
            end
        end
    end
    
    % Produce final tbnds value (recalling we focus on "symmetric" intervals)
    simpletbnds = [t_magnitudes , t_magnitudes];
    
end

function optval = entrywiseLPtbnds(objcoeffs,tbnds,extratol)
    % Minor variant of the entrywiseLP function (see below), for compatibility with the way we're using helper_minus when defining fmincon objective

    dimC = length(tbnds)/2;
    optval = entrywiseLP(objcoeffs,[tbnds(1:dimC);-tbnds(dimC+1:2*dimC)],extratol);

end

function optval = entrywiseLP(objcoeffs,upplowbnds,extratol)
    % Computes an analytic solution to an LP with entrywise constraints AND the constraint that the optimization variables sum to zero
    % objcoeffs: coefficients of objective function
    % upplowbnds: matrix with two rows, where the first row specifies upper bounds on the optimization variables and second row specifies lower bounds
    % extratol: extra relaxation parameter for the upper & lower bounds (note that this will only increase the LP value and hence still yield a reliable keyrate lower bound in the end)
    % (More explicitly: solves max(objcoeffs*x') s.t. sum(x) == 0 AND upplowbnds(2,j)-extratol <= x(j) <= upplowbnds(1,j)+extratol for all j)
    % Warning: assumes (without verifying) that the supplied entrywise constraints are feasible (after extratol relaxation) even together with the sum-to-zero constraint

    numterms = size(objcoeffs,2);
    uppbnds = upplowbnds(1,:) + extratol; lowbnds = upplowbnds(2,:) - extratol;
    
    % Basic sanity/feasibility checks: note however that we do not check if the LP is still feasible with the sum-to-zero constraint, which is a necessary assumption in order for the solution here to be valid.
    if size(upplowbnds) ~= [2 numterms]
        size(upplowbnds)
        size(objcoeffs)
        error('Inconsistent sizing!')
    elseif min(uppbnds - lowbnds) < 0
        upplowbnds
        min(uppbnds - lowbnds)
        error('Infeasible upper-lower bound pair!')
    end
    
    % Identify positive and negative coefficients in objective function
    posinds = find(objcoeffs>=0); neginds = find(objcoeffs<0); %It's basically arbitrary which group we include the coefficient-zero terms into
    if sort([posinds neginds]) ~= 1:numterms % Sanity check
        [posinds neginds]
        error('Problem with partitioning')
    end
    
    % Create initial "candidate solution" that maximizes the objective by saturating the entrywise constraints but does not satisfy sum-to-zero constraint
    % It will be iteratively changed over the upcoming process until it sums to zero
    candsol = zeros(1,numterms);
    for j=posinds
        candsol(j) = uppbnds(j);
    end
    for j=neginds
        candsol(j) = lowbnds(j);
    end
    
    % Creating copies of the index lists, which will be iteratively changed over the upcoming process
    % (Basically these will be lists of indices on which we have not yet changed the value of candsol; we will remove indices from these lists after each loop iteration)
    posindschange = posinds; negindschange = neginds; 

    % Iterative updating of candidate solution until it becomes feasible
    for changenum=1:numterms % Number of iterations should not exceed numterms if LP is feasible
        
        if sum(candsol) > 0 % If candidate solution currently sums to value larger than zero

            if isempty(posindschange) % Sanity check (this should never happen)
                error('Problem with iterative solution adjustment')
            end
            % Find smallest positive coefficient whose corresponding candsol term was not previously changed (this will be the term we change in this iteration)
            [mincoeff , minindsubset] = min(objcoeffs(posindschange));
            changeind = posindschange(minindsubset); % Index of that term
            if mincoeff<0 % Sanity check (this should never happen)
                error('Serious problem!')
            end
            % Remove that index from list of terms that could be updated in the future
            posindschange = posindschange(posindschange~=changeind);
            % Compute maximum amount we can decrease that candsol term while still satisfying constraint
            maxsolchange = uppbnds(changeind)-lowbnds(changeind);
            if maxsolchange > sum(candsol) % Implies that changing that candsol term suffices to make sum(candsol)==0
                % Changes that candsol term by the required amount to make sum(candsol)==0, then exit loop
                candsol(changeind) = candsol(changeind)-sum(candsol);
                break
            else 
                % Decrease that candsol term by the maximum amount and allow loop to continue to next term
                candsol(changeind) = lowbnds(changeind);
            end
    
        elseif sum(candsol)<0 % If candidate solution currently sums to value smaller than zero
    
            if isempty(negindschange) % Sanity check (this should never happen)
                error('Problem with iterative solution adjustment')
            end
            % Find largest negative coefficient whose corresponding candsol term was not previously changed (this will be the term we change in this iteration)
            [maxcoeff , maxindsubset] = max(objcoeffs(negindschange));
            changeind = negindschange(maxindsubset); % Index of that term
            if maxcoeff>=0 % Sanity check (this should never happen)
                error('Serious problem!')
            end
            % Remove that index from list of terms that could be updated in the future
            negindschange = negindschange(negindschange~=changeind);
            % Compute maximum amount we can increase that candsol term while still satisfying constraint
            maxsolchange = uppbnds(changeind)-lowbnds(changeind);
            if maxsolchange > abs(sum(candsol)) % Implies that changing that candsol term suffices to make sum(candsol)==0
                % Changes that candsol term by the required amount to make sum(candsol)==0, then exit loop
                candsol(changeind) = candsol(changeind)+abs(sum(candsol));
                break
            else
                % Increase that candsol term by the maximum amount and allow loop to continue to next term
                candsol(changeind) = uppbnds(changeind);
            end
    
        elseif sum(candsol)==0 % Implies candidate solution is already feasible
            break
        end
    
        if changenum==numterms %Should not be reached as long as the LP is feasible
            error('Too many solution adjustment iterations')
        end
    
    end
    
    %Sanity check for feasibility
    if abs(sum(candsol))>1e-10 || min(candsol-lowbnds)<0 || min(uppbnds-candsol)<0
        candsol
        sum(candsol)
        upplowbnds
        error('Created infeasible solution!')
    end
    
    optval = objcoeffs*candsol';

end

function [optsol , optval] = entrywiseLPcert(objcoeffs,upplowbnds,extratol)
    % Uses CVX to solve an LP with entrywise constraints AND the constraint that the optimization variables sum to zero
    % objcoeffs: coefficients of objective function
    % upplowbnds: matrix with two rows, where the first row specifies upper bounds on the optimization variables and second row specifies lower bounds
    % extratol: extra relaxation parameter for the upper & lower bounds (note that this will only increase the LP value and hence still yield a reliable keyrate lower bound in the end)
    % (More explicitly: solves max(objcoeffs*x') s.t. sum(x) == 0 AND upplowbnds(2,j)-extratol <= x(j) <= upplowbnds(1,j)+extratol for all j)

    numterms = size(objcoeffs,2);
    uppbnds = upplowbnds(1,:) + extratol; lowbnds = upplowbnds(2,:) - extratol;
    
    cvx_begin quiet
    cvx_precision high
    
    variable candsol(1,numterms)
    
    maximize(objcoeffs*candsol')
    subject to
    
    sum(candsol) == 0
    
    lowbnds <= candsol <= uppbnds
    
    cvx_end
    
    optsol = candsol;
    optval = cvx_optval;

end

function value = helper_minus(x)
    %helper function to make array only have positive entries
    value = zeros(size(x));
    for index = 1:length(x)
        if x(index) < 0
            value(index) = 0;
        else
            value(index) = x(index);
        end
    end
end

function prob = P_single_fail(N,p,uppert,lowert)
    %Calculates reject probability Pr(X >= N(p+uppert) or X <= N(p-lowert))
    % given a probability p for input t and N

    %N: Total number of rounds
    %p: probability
    %uppert, lowert: Acceptance parameters

    %tolerance to 0
    tol = 1e-14;

    %Upper range of the accepted parameters
    m1 = min([N , floor(N*(p+uppert))]);

    %Lower range of the accepted parameters
    m2 = max([0 , floor(N*(p-lowert))]);

    if p <= tol
        prob = 0;
    elseif abs(p-1) <= tol
        prob = 0;
    else
        % "Right tail" of distribution
        prob = 1 - betainc(1-p, helper_minus(N-m1), m1+1);
        % Add "left tail" of distribution
        prob = prob + betainc(1-p, helper_minus(N-m2+1), m2);
    end
end
