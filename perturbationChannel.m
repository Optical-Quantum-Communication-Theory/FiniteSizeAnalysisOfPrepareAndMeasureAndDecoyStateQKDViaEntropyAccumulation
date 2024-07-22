function rhoPrime = perturbationChannel(rho, perturbation)
% perturbationChannel Applies a small perturbation to rho if the smallest
% eigenvalue of rho is slightly negative. This makes rho more likely to be
% positive semidefinite (though the perturbation can fail if the minimum
% eigenvalue is too negative). Use perturbationChannelEpsilon to determine
% the optimal perturbation value
%
% Inputs:
% * rho: The density matrix to be perturbed
% * perturbation: the amount to perturb rho by (optional; will be
%   automatically calculated if not provided)
%
% Outputs:
% * rhoPrime: The resulting rho after perturbation
% * epsilon: The magnitude of the perturbation
%
% See also FW2StepSolver, perturbationChannelEpsilon 
arguments
    rho (:,:) double
    perturbation (1,1) double = NaN
end
    % if no perturbation value is specified, calculate a value
    if isnan(perturbation)
        % calculate a value to perturb by
        defaultPerturbation = 1e-14;
        dim = size(rho, 1);
        eigMin = lambda_min(rho);
        perturbation = 0;
        rho = (rho + rho')/2;
        if eigMin <= 0
            if real(trace(rho)) > 0
                perturbation = (eigMin  * dim)/(eigMin*dim - real(trace(rho))) + defaultPerturbation;
                 % check again if perturbation is in the range where the
                 %   theorem can apply
          	if perturbation < 0 || perturbation > 1/(exp(1)*(dim-1))
                ME = MException('perturbation_channel:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                fprintf("**** Error: Perturbation failed. The smallest eigenvalue is less than zero. ****\n")
                throw(ME);
           	end % end of if for validity of realEpsilon
            end
        end
    end
    % perturb by perturbation amount
    dim = size(rho, 1);
    rhoPrime = (1-perturbation) * rho + perturbation * real(trace(rho))*eye(dim)/dim;
    rhoPrime = (rhoPrime + rhoPrime')/2;

end