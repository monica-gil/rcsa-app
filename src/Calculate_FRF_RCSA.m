function [FRF_RCSA] = Calculate_FRF_RCSA(LongBarModel,InputIRCSA,EvaluateIRCSA)

%% Information
% FRF computation
% [H11]  = [h11 l11; n11 p11]
% [H12]  = [h12 l12; n12 p12]
% [H21]  = [h21 l21; n21 p21]
% [H22]  = [h22 l22; n22 p22]

%% Variable Information
% LongBarModel: structure containing the free-free dynamic model to evaluate.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.
% InputIRCSA: structure containing the FRF calculation configuration.
%   * FreqMin (Hz): minimum frequency of the FRF calculation.
%   * FreqMax (Hz): maximum frequency of the FRF calculation.
%   * DeltaFreq (Hz): frequency increment of the FRF calculation.
% EvaluateIRCSA: structure containing the interface indices for the FRF calculation.
%   * idxP1 (-): degree-of-freedom index of the first interface point.
%   * idxP2 (-): degree-of-freedom index of the second interface point.
% FRF_RCSA: structure containing h, l, n and p receptance terms and omega (rad/s).

%% Input Variables
MMatrix = LongBarModel.M;
KMatrix = LongBarModel.K;
CMatrix = LongBarModel.C;
FreqMin = InputIRCSA.FreqMin;
FreqMax = InputIRCSA.FreqMax;
DeltaFreq = InputIRCSA.DeltaFreq;
idxP1 = EvaluateIRCSA.idxP1;
idxP2 = EvaluateIRCSA.idxP2;

%% Maain
% Create frequency vector
omega = (FreqMin*2*pi):(DeltaFreq*2*pi):(FreqMax*2*pi);
nomega = length(omega);

% Initialize frequency response function y different degree of freedom
FRF_RCSA.h11 = zeros(nomega,1);
FRF_RCSA.h22 = zeros(nomega,1);
FRF_RCSA.h12 = zeros(nomega,1);
FRF_RCSA.h21 = zeros(nomega,1);

FRF_RCSA.l11 = zeros(nomega,1);
FRF_RCSA.l22 = zeros(nomega,1);
FRF_RCSA.l12 = zeros(nomega,1);
FRF_RCSA.l21 = zeros(nomega,1);

FRF_RCSA.n11 = zeros(nomega,1);
FRF_RCSA.n22 = zeros(nomega,1);
FRF_RCSA.n12 = zeros(nomega,1);
FRF_RCSA.n21 = zeros(nomega,1);

FRF_RCSA.p11 = zeros(nomega,1);
FRF_RCSA.p22 = zeros(nomega,1);
FRF_RCSA.p12 = zeros(nomega,1);
FRF_RCSA.p21 = zeros(nomega,1);

% Obtain the response function vector at each frequency
for idxFreq = 1:nomega
    FRF = inv((-(omega(idxFreq))^2)*MMatrix+j*(omega(idxFreq))*CMatrix+KMatrix);

    % hxx FRFs. Displacement-Displacement degree of freedom
    FRF_RCSA.h11(idxFreq,1) = FRF(idxP1-1,idxP1-1);
    FRF_RCSA.h22(idxFreq,1) = FRF(idxP2-1,idxP2-1);
    FRF_RCSA.h12(idxFreq,1) = FRF(idxP1-1,idxP2-1);
    FRF_RCSA.h21(idxFreq,1) = FRF(idxP2-1,idxP1-1);

    % lxx FRFs. Displacement-Rotation degree of freedom
    FRF_RCSA.l11(idxFreq,1) = FRF(idxP1-1,idxP1);
    FRF_RCSA.l22(idxFreq,1) = FRF(idxP2-1,idxP2);
    FRF_RCSA.l12(idxFreq,1) = FRF(idxP1-1,idxP2);
    FRF_RCSA.l21(idxFreq,1) = FRF(idxP2-1,idxP1);

    % nxx FRFs. Rotation-Displacement degree of freedom
    FRF_RCSA.n11(idxFreq,1) = FRF(idxP1,idxP1-1);
    FRF_RCSA.n22(idxFreq,1) = FRF(idxP2,idxP2-1);
    FRF_RCSA.n12(idxFreq,1) = FRF(idxP1,idxP2-1);
    FRF_RCSA.n21(idxFreq,1) = FRF(idxP2,idxP1-1);

    % pxx FRFs. Rotation-Rotation degree of freedom
    FRF_RCSA.p11(idxFreq,1) = FRF(idxP1,idxP1);
    FRF_RCSA.p22(idxFreq,1) = FRF(idxP2,idxP2);
    FRF_RCSA.p12(idxFreq,1) = FRF(idxP1,idxP2);
    FRF_RCSA.p21(idxFreq,1) = FRF(idxP2,idxP1);
end
FRF_RCSA.omega = omega;
FRF_RCSA.Freq = omega/(2*pi);

end % end of function
