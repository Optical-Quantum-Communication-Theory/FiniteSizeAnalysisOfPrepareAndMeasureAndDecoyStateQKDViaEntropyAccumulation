function dfmat = gradientCorChoi(ChoiMat,crossover,stateTestRounds,observables,testprob,nu,dimA,dimAprime,dimB,dimR)

    m = size(crossover,1);
    n = size(crossover,2);

    rho = PartialTrace(kron(eye(dimA),ChoiMat)*(kron(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),eye(dimB))),[2],[dimA,dimAprime,dimB]);
    
    genfrequency = zeros(m,n);

    for indexrow = 1:m
        for indexcolumn = 1:n
            genfrequency(indexrow,indexcolumn) = trace(observables{indexrow,indexcolumn}'*rho);
        end
    end

    dVqmat = dVarq(genfrequency,crossover,testprob,dimR);

    derF = zeros(dimAprime*dimB);
    SwappedStates = Swap(PartialTranspose(stateTestRounds,[2],[dimA,dimAprime]),[1,2],[dimA,dimAprime]);

    for indexrow = 1:m
        for indexcolumn = 1:n
            B = PartialTrace(kron(eye(dimAprime),observables{indexrow,indexcolumn})*kron(SwappedStates,eye(dimB)),[2],[dimAprime,dimA,dimB]);

            derF = derF + transpose(B)*(crossover(indexrow,indexcolumn) + nu/(1-nu)*log(2)/2*dVqmat(indexrow,indexcolumn));
        end
    end
    dfmat = -derF;
end


function dfmat = dVarq(genfrequency,crossover,testprob,dimR)
    

    m = size(crossover,1);
    n = size(crossover,2);

    varp = 1/testprob*sum(genfrequency.*(max(crossover, [],'all') - crossover).^2,'all') - (max(crossover, [],'all') - sum(crossover.*genfrequency,'all'))^2;

    gradvarp = zeros(m,n);
    for indexrow = 1:m
        for indexcolumn = 1:n
            gradvarp(indexrow,indexcolumn) = 1/testprob*(max(crossover, [],'all') - crossover(indexrow,indexcolumn))^2 + 2*(max(crossover, [],'all') - sum(crossover.*genfrequency,'all'))*crossover(indexrow,indexcolumn);
        end
    end
    dfmat = (log2(1+2*dimR) + sqrt(2 + varp)) * gradvarp/sqrt(2 + varp);
end