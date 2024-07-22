eta = 0.398107170553497;
f_EC=1.16;
misalignment = 0.01;
depol = 0.01;
n_signals = 9;
epsilon_sec = 1e-8;
decoy_intens = [0.9,0.02,0.001];
decoy_probs = [1/3,1/3,1/3];
n_photon = 10;

%Run decoy
% [optdeltacomp_decoy, optkeydecoy] = DecoyBB84_keyrate_EAT(eta, f_EC, misalignment, n_signals, epsilon_sec,decoy_intens,decoy_probs,n_photon);

%Run qubit
[optdeltacomp_qubit, optkeyratequbit] = QubitBB84_keyrate_EAT (eta, f_EC, depol, n_signals, epsilon_sec);