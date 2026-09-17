function [ToolCoupledRCSAFRF] = RCSA_Response(AnalyticFRF,ExpFRF,InputIRCSA)

%% Variable Information
% AnalyticFRF: structure containing h, l, n and p free-free receptance terms of the tool.
% ExpFRF: structure containing h33, l33, n33 and p33 receptance terms of the holder.
% InputIRCSA: structure containing the IRCSA frequency configuration.
%   * FreqMax (Hz): maximum frequency of the coupled RCSA response.
%   * DeltaFreq (Hz): frequency increment of the coupled RCSA response.
% ToolCoupledRCSAFRF: two-column array with complex tooltip receptance (m/N) and frequency (Hz).

%% Input Variables
InputFreqMax = InputIRCSA.FreqMax;
DeltaFreq = InputIRCSA.DeltaFreq;

%% Main
ExpFreqMax = (length(ExpFRF.h33)-1)*DeltaFreq;
AnalyticFreqMax = (length(AnalyticFRF.h11)-1)*DeltaFreq;
if InputFreqMax ~= ExpFreqMax
   warning = strjoin({'Experimental frequency range differs from', num2str(InputFreqMax)});
   display(warning)   
end
if InputFreqMax > AnalyticFreqMax
   warning = strjoin({'Increase the maximum frequency of analytical free-free FRFs to', num2str(MaxFrec)});
   display(warning)   
end

% Frequency vector with the precision and maximum frequency defined.
w_int = 0:DeltaFreq:InputFreqMax;
NFreq = length(w_int);

% FRF matrix definition by analitical free-free condition (H11, H22, H12 and
% H21) and experimental characterization of the joint point (H33)
ToolCoupledRCSAFRF.G11 = zeros(NFreq,1);
for idxFreq = 1:NFreq
    H11 = [AnalyticFRF.h11(idxFreq,1) AnalyticFRF.l11(idxFreq,1); AnalyticFRF.n11(idxFreq,1) AnalyticFRF.p11(idxFreq,1)];
    H22 = [AnalyticFRF.h22(idxFreq,1) AnalyticFRF.l22(idxFreq,1); AnalyticFRF.n22(idxFreq,1) AnalyticFRF.p22(idxFreq,1)];
    H12 = [AnalyticFRF.h12(idxFreq,1) AnalyticFRF.l12(idxFreq,1); AnalyticFRF.n12(idxFreq,1) AnalyticFRF.p12(idxFreq,1)];
    H21 = [AnalyticFRF.h21(idxFreq,1) AnalyticFRF.l21(idxFreq,1); AnalyticFRF.n21(idxFreq,1) AnalyticFRF.p21(idxFreq,1)];
    H33 = [ExpFRF.h33(idxFreq,1) ExpFRF.l33(idxFreq,1); ExpFRF.n33(idxFreq,1) ExpFRF.p33(idxFreq,1)];
    % RCSA operation to obtain the output response
    G11 = H11-H12*inv(H22+H33)*H21;
    ToolCoupledRCSAFRF.G11(idxFreq,1) = G11(1,1);
end
ToolCoupledRCSAFRF.Freq = w_int;

end % end of function
