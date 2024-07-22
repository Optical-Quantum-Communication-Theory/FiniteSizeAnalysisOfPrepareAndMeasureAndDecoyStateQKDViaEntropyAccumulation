function [optdeltacomp_decoy,optkeydecoy] = DecoyBB84_keyrate_EAT (eta, f_EC, misalignment, n_signals, epsilon_sec,decoy_intens,decoy_probs,n_photon)
    %Input:
    %
    %eta : loss parameter could be a vector if needed to plot over various
    %loss parameter
    %
    %f_EC : error correction efficiency
    %
    %misalignment : misalignment parameter
    %
    %n_signal : vector containing the exponent of number of signals to be
    %considered (eg. n_signal = 11:1:12 means the numebr of singals to be
    %considered are 10^11,10^12,10^13 )
    %
    %epsilon_sec : security parameter
    %
    %decoy_intens : vector containing decoy intensities (in this case the
    %intensities are 3)
    %
    %decoy_probs : vector containing probability of sending each decoy pulse
    %
    %n_photon = photon number cut-off
    %
    %Output:
    %
    %optkeydecoy : vector containing keyrate values
    %optdeltacomp_decoy : vector containing the values of completeness parameter

    dimA = 2;
    dimAprime = 2;
    dimB = 3;
    dimR = 2;

    %Channel const
    ed = 0;

    s = 1;

    %Decoy parameters:
    decoys = decoy_intens;
    probsDecoy = decoy_probs;
    pd = 0;

    %delta_completeness parameters
    optdeltacomp_decoy = zeros(length(eta),length(n_signals));
    uppbndstest = cell(length(eta),length(n_signals));
    lowbndstest = cell(length(eta),length(n_signals));
    uppbndgen = zeros(length(eta),length(n_signals));
    lowbndgen = zeros(length(eta),length(n_signals));
    pabortATbnd = zeros(length(eta),length(n_signals));

    %defining vectors containing the final results
    optkeydecoy = zeros(length(eta),length(n_signals));
    optgammadecoy = zeros(length(eta),length(n_signals));
    optnudecoy = zeros(length(eta),length(n_signals));

    pre_keyrate = 0;
    initial_testProb = 0.1;

    %grid search resolution
    grid = 50;
    gammat = 10.^(-linspace(0.1,9,grid));
    nut = 10.^(-linspace(1,11,grid));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for n = n_signals
    disp(n)
    maxElement = 1;
for l = 1:length(eta)
    ns = 10^n;
if maxElement <= 0
    break;
end
 k = zeros(length(gammat),length(nut));
i=1;
pre_keyrate = 0;
for gam = 10.^(-linspace(0.1,9,grid))
    j=1;
    for nu1 = 10.^(-linspace(1,11,grid))%+gam/2,4,100))
        if k(i,j) ~= -1
            try      
                [pre_keyrate,~,~,~,pre_ot,pre_fst_opt] = ConnectorsDecoy (gam,dimA,dimB,dimAprime,dimR,nu1,eta(l), ...
            pd,misalignment,probsDecoy(1),probsDecoy(2),probsDecoy(3), ...
            decoys(1),decoys(2),decoys(3),epsilon_sec, ...
            f_EC,ns,n_photon);
            catch
                pre_keyrate = -2;
                k(i,j) = -2;
                cvx_begin
                cvx_clear
            end
        end
        if pre_keyrate < 0 || k(i,j) == -1
            if j <= 1 && i <= length(gammat)-1
                k(i+1,j) = -1;
            end
            if j > 1 && k(i,j-1) > 0 && k(i,j) ~= -2
                break
            end
            if j > 1 && k(i,j-1) < 0 && k(i,j) ~= -2 && i <= length(gammat)-1
                k(i+1,j) = -1;
            end
               
        end
        if j > 1
            if k(i,j) == -1 && k(i,j-1) > 0
            end
            if pre_keyrate < 0 && pre_keyrate < k(i,j-1) && pre_keyrate ~= -2
                break;
            end
        end
        if k(i,j) ~= -1
            k(i,j) = pre_keyrate;
        end
        if j > 1
            if k(i,j) < k(i,j-1) && k(i,j) ~= -2
                break;
            end
        end
        j = j + 1;
    end
    if i > 1
        if all (k(i-1,:) <= 0)
            if all (k(i,:) <= k(i-1,:))
        %       break;  %(add to decrese the run time)
            end
        end
    end
    i = i + 1;
end
[maxElement, linearIndex] = max(k(:));
[row, column] = ind2sub(size(k), linearIndex);
try
[o2,flag2,crossover,acceptfrequency] = ConnectorsDecoy (gammat(row),dimA,dimB,dimAprime,dimR,nut(column),eta(l), ...
            pd,misalignment,probsDecoy(1),probsDecoy(2),probsDecoy(3), ...
            decoys(1),decoys(2),decoys(3),epsilon_sec, ...
            f_EC,ns,n_photon);
catch

end


[iqw,optdeltacomp_temp , uppbndstest_temp , lowbndstest_temp , uppbndgen_temp , lowbndgen_temp , pabortATbnd_temp] = optdeltacom(ns,gammat(row),acceptfrequency,crossover,1e-3,100);
if iqw == 1
    [iqw,optdeltacomp_temp , uppbndstest_temp , lowbndstest_temp , uppbndgen_temp , lowbndgen_temp , pabortATbnd_temp] = optdeltacom(ns,gammat(row),acceptfrequency,crossover,1e-3,0);
end
%Storing completeness parameters
optdeltacomp_decoy(l,s) = optdeltacomp_temp;
uppbndstest{l,s} = uppbndstest_temp;
lowbndstest{l,s} = lowbndstest_temp;
uppbndgen(l,s) = uppbndgen_temp;
lowbndgen(l,s) = lowbndgen_temp;
pabortATbnd(l,s) = pabortATbnd_temp;

optkeydecoy(l,s) =  maxElement-optdeltacomp_temp;
%Checking the flag to see if the FW returned the solved solution
if ~strcmp(flag2, "Solved")
    disp('Variable is not equal to "Solved". Pausing the run.');
end
%calculating the keyrates along with correction terms

optgammadecoy(l,s) = gammat(row);
optnudecoy(l,s) = nut(column);
end
s = s + 1;
end

end