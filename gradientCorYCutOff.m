function dfvec = gradientCorYCutOff(Yield, YieldCutOff, crossover,testprob,nu,dimR,probsDecoy,decoys,n_photon)

    m = size(crossover,1);
    n = size(crossover,2);
    n_decoy = size(crossover,3);
    
    %Calculate the frequency/ probability distribution generated 
    genfrequency = zeros(m,n,n_decoy);
    
    for indexInt = 1:n_decoy
        pmu = Poisson(decoys(indexInt),0:1:n_photon);
        Pmu = kron(eye(m),pmu);

        genfrequency(:,:,indexInt) = probsDecoy(indexInt)*(Pmu*Yield + YieldCutOff(:,:,indexInt));
    end
    
    %Calculate corresponding variance term
    varp = 1/testprob*sum(genfrequency.*(max(crossover, [],'all') - crossover).^2,'all') - (max(crossover, [],'all') - sum(crossover.*genfrequency,'all'))^2;
    
    %Calculate gradient of variance
    gradvarp = zeros(m,n,n_decoy);
    for indexInt = 1:n_decoy
        for indexrow = 1:m
            for indexcolumn = 1:n
                gradvarp(indexrow,indexcolumn,indexInt) = 1/testprob*(max(crossover, [],'all') - crossover(indexrow,indexcolumn,indexInt))^2 + 2*(max(crossover, [],'all') - sum(crossover.*genfrequency,'all'))*crossover(indexrow,indexcolumn,indexInt);
            end
        end
    end
    
    %calculate gradient for each intensity
    
        
    %Preallocate full gradient for each intensity
    derF = zeros(m,n,n_decoy);

    %sum over each intensity
    for indexInt = 1:n_decoy
        %gradient of V for intensity mu
        dVarqInt = (log2(1+2*dimR) + sqrt(2 + varp)) * gradvarp(:,:,indexInt)/sqrt(2 + varp);

        %Full gradient for each intensity added to derF
        derF(:,:,indexInt) = probsDecoy(indexInt)*(nu/(1-nu)*log(2)/2*dVarqInt + crossover(:,:,indexInt));
    end
    
    %Add graient to cell array
    dfvec = vec(-derF);
end

%%%%%%%%%%%%%%%%%%%%% Possonian distribution %%%%%%%%%%%%%%%%%%%%%%%%%
function prob = Poisson(mu,n)
    prob = exp(-mu).*mu.^n./factorial(n);
end