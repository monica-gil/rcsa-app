function FRF_Point = Calculate_FRF(ToolModelRestrained,InputToolRestrain,EvaluateToolRestrain)

%% Variable Information
% ToolModelRestrained: structure containing the restrained dynamic model.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.
% InputToolRestrain: structure containing the restrained-FRF frequency configuration.
%   * FreqMin (Hz): minimum frequency of the response calculation.
%   * FreqMax (Hz): maximum frequency of the response calculation.
%   * DeltaFreq (Hz): frequency increment of the response calculation.
% EvaluateToolRestrain: structure containing the restrained-FRF evaluation indices.
%   * idxP1 (-): response degree-of-freedom index for the first FRF.DeltaPreq
%   * idxP2 (-): response degree-of-freedom index for the optional second FRF; 0 disables it.
% FRF_PP1: structure containing the FRF evaluated at idxP1.
%   * Magnitude (m/N): magnitude of the receptance at each frequency.
%   * Real (m/N): real part of the receptance at each frequency.
%   * omega (rad/s): angular-frequency vector used for the calculation.
% FRF_PP2: structure containing the FRF evaluated at idxP2.
%   * Magnitude (m/N): magnitude of the receptance at each frequency; it remains zero when idxP2 is 0.
%   * Real (m/N): real part of the receptance at each frequency; it remains zero when idxP2 is 0.
%   * omega (rad/s): angular-frequency vector used for the calculation.

%% Input Variables
MMatrix = ToolModelRestrained.M;
KMatrix = ToolModelRestrained.K;
CMatrix = ToolModelRestrained.C;
FreqMin = InputToolRestrain.FreqMin;
FreqMax = InputToolRestrain.FreqMax;
DeltaPreq = InputToolRestrain.DeltaFreq;
% idxP1 = EvaluateToolRestrain.idxP1;
% idxP2 = EvaluateToolRestrain.idxP2;
PointID = EvaluateToolRestrain.idxPoint; % idx of the point in M, K, C matrixes

NPoint = length(PointID); % number of points where the FRF is evaluated

%% Main
% Create frequency vector
omega = (FreqMin*2*pi):(DeltaPreq*2*pi):(FreqMax*2*pi);
NOmega = length(omega);

% Initialize frequency response function and different degree of freedom
% FRF_PP1.Magnitude = zeros(NOmega,1);
% FRF_PP2.Magnitude = zeros(NOmega,1);
% FRF_PP1.Real = zeros(NOmega,1);
% FRF_PP2.Real = zeros(NOmega,1);
FRF_Point.Magnitude = NaN(NOmega,NPoint);
FRF_Point.Real = NaN(NOmega,NPoint);

% Obtain the response function vector at each frequency

for idxFreq = 1:NOmega
    FRF = inv((-(omega(idxFreq))^2)*MMatrix+j*(omega(idxFreq))*CMatrix+KMatrix);
    for idxPoint = 1:NPoint
        if PointID(idxPoint) > 0
            FRF_Point.Magnitude(idxFreq,idxPoint) = abs(FRF(PointID(idxPoint)-1,PointID(idxPoint)-1));
            FRF_Point.Real(idxFreq,idxPoint) = real(FRF(PointID(idxPoint)-1,PointID(idxPoint)-1));
        end
        % FRF_PP1.Magnitude(idxFreq,1) = abs(FRF(idxP1-1,idxP1-1));
        % FRF_PP1.Real(idxFreq,1) = real(FRF(idxP1-1,idxP1-1));
        % if idxP2 > 0
        %     FRF_PP2.Magnitude(idxFreq,1) = abs(FRF(idxP2-1,idxP2-1));
        %     FRF_PP2.Real(idxFreq,1) = real(FRF(idxP2-1,idxP2-1));
        % end
    end
end
% FRF_PP1.omega = omega;
% FRF_PP2.omega = omega;
% FRF_PP1.Freq = omega./(2*pi);
% FRF_PP2.Freq = omega./(2*pi);
FRF_Point.Freq = omega./(2*pi);

end % end of function
