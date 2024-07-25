%Transmitivity
eta_List = 10.^(-linspace(0,4,21));
eta_List(1) = 10^(-2);

%Error correction efficiency
f_EC=1.16;

%Misalignment (only used for decoy)
misalignment = 0.01;

%Depolarization
depol = 0.01;

%Security parameter
epsilon_sec = 1e-8;

% Decoy intensities
decoy_intens = [0.9,0.02,0.001];

%Probabilities of sending each decoy intensity given a test round, i.e
%p(decoy|test)
decoy_probs = [1/3,1/3,1/3];

%Photonnumber cutoff
n_photon = 10;

%List of Exponent of number of signals sent, i.e. n = 10^(n_signals)
N_list = [7,8,10,12];

%Choose between unique acceptance and realistic acceptance sets
uniqueAcc = false;

%% run the plot

%Store empty matrix of key rates
KeyRate = zeros(numel(N_list),numel(eta_List));

for indexSignal = 1:numel(N_list)
    for indexEta = 1:numel(eta_List)
        %N signal 
        n_signals = N_list(indexSignal);
        %transmittivity
        eta = eta_List(indexEta);
        [optdeltacomp_decoy, optkeydecoy] = DecoyBB84_keyrate_EAT(eta, f_EC, misalignment, n_signals, epsilon_sec,decoy_intens,decoy_probs,n_photon);
        %store result in matrix. Each row corresponds to one value of
        %n_signal
        if uniqueAcc == true
            KeyRate(indexSignal,indexEta) = optkeydecoy + optdeltacomp_decoy;
        else
            KeyRate(indexSignal,indexEta) = optkeydecoy;
        end
    end
end

