function C = DampingMatrixCalculationRealPart(M,K,Psi)

%% Variable Information
% M (kg): mass matrix of the finite-element model.
% K (N/m): stiffness matrix of the finite-element model.
% Psi (%): modal damping ratio expressed as a percentage.
% C (N*s/m): damping matrix assembled from the real parts of modal frequencies.

[FreqNat,MMr,VN] = ModalParameters_Calculation_FRF(M,K,0);

NFreqNat = length(FreqNat);
CMN = zeros(NFreqNat,NFreqNat);
for idxFreq = 1:NFreqNat
    CMN(idxFreq,idxFreq) = 2*(Psi/100)*real(FreqNat(idxFreq)); % 2*(Psi/100)*real(FreqNat(idxFreq)) %!!
end

% C matrix is obtained operating with the normalized damping matrix (CMN):
C = inv((VN)')*(CMN)*inv(VN);

end % end of function
