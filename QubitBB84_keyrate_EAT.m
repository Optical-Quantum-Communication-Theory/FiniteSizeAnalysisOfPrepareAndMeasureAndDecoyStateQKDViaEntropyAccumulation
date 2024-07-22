function [optdeltacomp,optkeyrate] = QubitBB84_keyrate_EAT (eta, f_EC, depol, n_signals, epsilon_sec)
    %Inputs:
    %
    %eta : loss parameter could be a vector if needed to plot over various
    %loss parameter
    %
    %f_EC : error correction efficiency
    %
    %depol : depolarization parameter
    %
    %n_signal : vector containing the exponent of number of signals to be
    %considered (eg. n_signal = 11:1:12 means the numebr of singals to be
    %considered are 10^11,10^12,10^13 )
    %
    %epsilon_sec : security parameter

    %Output:
    %
    %optkeyrate : vector containing keyrate values
    %optdeltacomp : vector containing the values of completeness parameter

    dimA = 2;
    dimAprime = 2;
    dimB = 3;
    dimR = 2;

    %Channel const
    ed = 0;

    s = 1;

    %defining vectors containing the final results
    optdeltacomp = zeros(length(eta),length(n_signals));
    optkeyrate = zeros(length(eta),length(n_signals));
    optgamma = zeros(length(eta),length(n_signals));
    optnu = zeros(length(eta),length(n_signals));
    pre_keyrate = 0;
    initial_testProb = 0.1;

    %grid search resolution
    grid = 50;
    gammat = 10.^(-linspace(0.1,9,grid));
    nut = 10.^(-linspace(1,12,grid));
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = n_signals
    maxElement = 1;
for l = 1:length(eta)
    ns = 10^n;
if maxElement <= 0
    break;
end
 k = zeros(length(gammat),length(nut));
 dc = zeros(length(gammat),length(nut));
i=1;
pre_keyrate = 0;

%Perform grid search over gamma and nu
for gam = 10.^(-linspace(0.1,9,grid))
    j=1;
    for nu1 = 10.^(-linspace(1,12,grid))
        if k(i,j) ~= -1
            [pre_keyrate,pre_comp] = Connectors (gam,dimA,dimB,dimAprime,dimR,nu1,eta(l), ...
            depol,epsilon_sec,f_EC,ns);
        end
        if pre_keyrate < 0 || k(i,j) == -1
            if j <= 1
                k(i+1,j) = -1;
            end
            if j > 1 && k(i,j-1) > 0
                break
            end
            if j > 1 && k(i,j-1) < 0
                k(i+1,j) = -1;
            end
        end
        if j > 1
            if k(i,j) == -1 && k(i,j-1) > 0
            end
            if pre_keyrate < 0 && pre_keyrate < k(i,j-1)
                break;
            end
        end
        if k(i,j) ~= -1
            k(i,j) = pre_keyrate;
            dc(i,j) = pre_comp;
        end
        if j > 1
            if k(i,j) < k(i,j-1) 
                break;
            end
        end
        j = j + 1;
    end
    if i > 1
        if all (k(i-1,:) <= 0)
            if all (k(i,:) <= k(i-1,:))
            end
        end
    end
    i = i + 1;
end
[maxElement, linearIndex] = max(k(:));
[row, column] = ind2sub(size(k), linearIndex);
optdeltacomp_temp = dc(row,column);
optdeltacomp(l,s) = optdeltacomp_temp;
optkeyrate(l,s) =  maxElement;

optgamma(l,s) = gammat(row);
optnu(l,s) = nut(column);
disp(['(', num2str(l), ',', num2str(n), ')']);

end
s = s + 1;
    end
end