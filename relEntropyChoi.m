%%%%%%%%%%%%%%%%%%%%% Relative entropy %%%%%%%%%%%%%%%%%%%%%%%%%

function fval = relEntropyChoi(ChoiMat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB)
    %Check if J is PSD
    minEigJ = eigs(ChoiMat,1,'smallestreal');
        
    %Calculate rho
    rho = PartialTrace(kron(eye(dimA),ChoiMat)*(kron(PartialTranspose(stateGenRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
    
    %check validity of rho and perform perturbation if not valid
    rho = perturbationChannel(rho);
    
    % if there is no Kraus operator, then proceed the calculation without the G
    % map.
    if nargin == 2 || isempty(krausOperators)
        
        zRho = 0;
        for jMapElement = 1:numel(keyMap)
            zRho = zRho + keyMap{jMapElement}*rho*keyMap{jMapElement};
        end % calculate the Z(\rho)
        
        %check validity of zRho and perform perturbation if not valid
        zRho = perturbationChannel(zRho);
    
        fval = real(trace(rho*(logm(rho)-logm(zRho)))); % calculate the quantum relative entropy
    else
        % for the case there is a post-selection map.
        
        gRho = krausFunc(rho,krausOperators); % calculate G(\rho).
        
        %check validity of gRho and perform perturbation if not valid
        gRho = perturbationChannel(gRho);
    
        zRho = 0;
        for jMapElement = 1:numel(keyMap)
            zRho = zRho + keyMap{jMapElement}*gRho*keyMap{jMapElement};
        end %calculate the Z(G(\rho))
        
        %check validity of zRho and perform perturbation if not valid
        zRho = perturbationChannel(zRho);
        
        fval = real(trace(gRho*(logm(gRho)-logm(zRho)))); % calculate the quantum relative entropy
    end

end