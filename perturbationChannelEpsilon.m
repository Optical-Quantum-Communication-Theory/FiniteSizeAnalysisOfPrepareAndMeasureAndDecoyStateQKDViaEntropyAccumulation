function epsilon = perturbationChannelEpsilon(rho)
% PERTURBATIONCHANNELEPSILON Computes the necessary epsilon by which to
% perturb according to Theorem 2 of Reliable numerical key rates paper by
% Winick, Lutkenhaus, and Coles. 
%
% Inputs:
% * rho: The density matrix to be perturbed
%
% Outputs:
% * epsilon: The magnitude of the perturbation
%
% See also FW2StepSolver

    defaultPerturbation = 1e-14;
    dim = size(rho,1);
    eigMin = lambda_min(rho);
    epsilon=0;
    if eigMin<=0
        if real(trace(rho))  > 0 
            epsilon = (eigMin  * dim)/(eigMin*dim - real(trace(rho))) + defaultPerturbation;
             % check if epsilon is in the range where the theorem can apply
          	if epsilon < 0 || epsilon > 1/(exp(1)*(dim-1))
                    ME = MException('perturbationChannelEpsilon:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                    fprintf("**** Error: Perturbation calculation failed. The smallest eigenvalue is less than zero. ****\n")
                    throw(ME);
           	end % end of if for validity of realEpsilon
        else
             ME = MException('perturbationChannelEpsilon:badRho','Trace of rho should be positive');
             throw(ME);
       	end % end of check perturbed rho
    end

end