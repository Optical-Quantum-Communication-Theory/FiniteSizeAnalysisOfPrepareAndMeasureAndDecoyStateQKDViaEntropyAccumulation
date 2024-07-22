%%  FUNCTION NAME: relEntropyChoi
%
% % $f(\rho) := D(\mathcal{G}(\rho)||\mathcal{Z}(\mathcal{G}(\rho)))$
%
% Syntax:  dfval = gradientRelEntChoi(ChoiMatrix,keyMap,krausOperators,stateGenRounds,dimA,dimAprime,dimB)
%
% Input: 
%
%  * ChoiMatrix - Choi Matrix describing the channel between Alice and Bob
%
%  * stateGenRound - states sent by alice in generation rounds
% 
%  * keyMap - Alice's key map PVM (If Alice's key map POVM is not projective, use Naimark's extension)
%
%  * krausOperators - The Kraus operators for the post-selection map of
%  Alice and Bob.
%       dimA:               dimension of Alice's quantum system
%       dimAprime:          dimension of quantum system sent to Bob
%       dimB:               dimension of Bob's quantum system
%
% Output:
%  
%  * dfval - the gradient. 
%%


function dfval = gradientRelEntChoi(Choimat,stateGenRounds,keyMap,krausOperators,dimA,dimAprime,dimB)

    %Calculate rho in gen rounds, i.e. rho = Phi(J)
    rho = PartialTrace(kron(eye(dimA),Choimat)*(kron(PartialTranspose(stateGenRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);

    if nargin == 2 || isempty(krausOperators)
        % if there is no post-selection map
        
        rho = perturbationChannel(rho);
        
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*rho*keyMap{j};
        end
        
        zRho = perturbationChannel(zRho);
        
        %Compute partial transposes of rho and zrho
        GlogGrhoPT = PartialTranspose(logm(rho),[2],[dimA,dimB]);
        GlogZGrhoPT = PartialTranspose(logm(zRho),[2],[dimA,dimB]);
        
        %Phi^dagger o log o rho
        dfval1 = PartialTrace(kron(eye(dimAprime),GlogGrhoPT)*(kron(transpose(stateGenRounds),eye(dimB))),[3],[dimA,dimAprime,dimB]);

        %Phi^dagger o log o Z o rho
        dfval2 = PartialTrace(kron(eye(dimAprime),GlogZGrhoPT)*(kron(transpose(stateGenRounds),eye(dimB))),[3],[dimA,dimAprime,dimB]);
        
        %gradf = Phi^dagger o log o rho - Phi^dagger o log o Z o rho = dfval1 -dfval2
        dfval = dfval1 - dfval2;

    else
        % if there is a post-selection map.
        gRho = krausFunc(rho,krausOperators);
        
        %check validity of gRho and perform perturbation if not valid
        gRho = perturbationChannel(gRho);
    
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*gRho*keyMap{j};
        end
        
        %check validity of zRho and perform perturbation if not valid
        zRho = perturbationChannel(zRho);


        %Compute G^dagger o log o G(rho) 
        GdaggerRho = krausFunc(logm(gRho),krausOperators,'transpose');

        %Compute partial transpose of G^dagger o log o G(rho)
        GlogGrhoPT = PartialTranspose(GdaggerRho,[2],[dimA,dimB]);
        
        %Phi^dagger o G^dagger o log o G(rho)
        dfval1 = PartialTrace(kron(eye(dimAprime),GlogGrhoPT)*(kron(transpose(stateGenRounds),eye(dimB))),[2],[dimAprime,dimA,dimB]);
        

        %Compute G^dagger o log o Z o G(rho)
        GdaggerZRho = krausFunc(logm(zRho),krausOperators,'transpose');

        %Compute partial transpose of G^dagger o log o Z o G(rho)
        GlogZGrhoPT = PartialTranspose(GdaggerZRho,[2],[dimA,dimB]);

        %Phi^dagger o G^dagger o log o Z o G(rho)
        dfval2 = PartialTrace(kron(eye(dimAprime),GlogZGrhoPT)*(kron(transpose(stateGenRounds),eye(dimB))),[2],[dimAprime,dimA,dimB]);
        

        %gradf = Phi^dagger o G^dagger o log o G(rho) - Phi^dagger o G^dagger o log o Z o G(rho) = dfval1 -dfval2
        dfval = dfval1 - dfval2;
    end

end