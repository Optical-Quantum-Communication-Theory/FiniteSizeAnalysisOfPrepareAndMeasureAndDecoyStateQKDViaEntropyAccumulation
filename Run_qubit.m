%Transmitivity
eta_List = 10.^(-linspace(0,4,21));
eta_List(1) = 10^(-2);

%Error correction efficiency
f_EC=1.16;

%Depolarization
depol = 0.01;

%Security parameter
epsilon_sec = 1e-8;

%List of Exponent of number of signals sent, i.e. n = 10^(n_signals)
N_list = [6,8,10,12];

%Choose between unique acceptance and realistic acceptance sets
uniqueAcc = false;

%% run the plot

%Store empty matrix of key rates
KeyRate = zeros(N_list,numel(eta_List));

for indexSignal = 1:numel(N_list)
    for indexEta = 1:numel(eta_List)
        %N signal 
        n_signals = N_list(indexSignal);
        %transmittivity
        eta = eta_List(indexEta);
        [optdeltacomp_qubit, optkeyratequbit] = QubitBB84_keyrate_EAT (eta, f_EC, depol, n_signals, epsilon_sec);
        %store result in matrix. Each row corresponds to one value of
        %n_signal
        if uniqueAcc == true
            KeyRate(indexSignal,indexEta) = optkeyratequbit + optdeltacomp_qubit;
        else
            KeyRate(indexSignal,indexEta) = optkeyratequbit;
        end
    end
end

