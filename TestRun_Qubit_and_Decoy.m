%Transmitivity
eta = 0.398107170553497;

%Error correction efficiency
f_EC=1.16;

%Misalignment (only used for decoy)
misalignment = 0.01;

%Depolarization
depol = 0.01;

%Exponent of number of signals sent, i.e. n = 10^(n_signals)
n_signals = 9;

%Security parameter
epsilon_sec = 1e-8;

% Decoy intensities
decoy_intens = [0.9,0.02,0.001];

%Probabilities of sending each decoy intensity given a test round, i.e
%p(decoy|test)
decoy_probs = [1/3,1/3,1/3];

%Photonnumber cutoff
n_photon = 10;

%%

%Run qubit
[optdeltacomp_qubit, optkeyratequbit] = QubitBB84_keyrate_EAT (eta, f_EC, depol, n_signals, epsilon_sec);

%Run decoy
% [optdeltacomp_decoy, optkeydecoy] = DecoyBB84_keyrate_EAT(eta, f_EC, misalignment, n_signals, epsilon_sec,decoy_intens,decoy_probs,n_photon);
