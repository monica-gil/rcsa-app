function Holder_IRCSA_FRF = CalculateIRCSA(ExpFRF,AnalyticFRF,InputIRCSA)

%% Information
% [H33] = [h33 l33; n33 p33]
% [G11] = [h11 -; - -]
% [G12] = [h12 -; - -]

%% Variable Information
% ExpFRF: structure containing experimental receptances.
%   * h11, h12, h33: two-row arrays with frequency (Hz) in row 1 and complex receptance (m/N) in row 2.
% AnalyticFRF: structure containing the free-free analytical receptances of the LongBarData.
%   * h11, h12, h21, h22, n12, n21, n22 and p22: complex receptance vectors.
% InputIRCSA: structure containing the IRCSA configuration.
%   * FreqMin (Hz): requested minimum frequency for IRCSA.
%   * FreqMax (Hz): requested maximum frequency for IRCSA.
%   * DeltaFreq (Hz): frequency increment used to interpolate the experimental FRFs.
%   * signIRCSA (-): sign convention; 0 preserves the experimental FRF sign and 1 reverses it.
% Holder_IRCSA_FRF: structure containing h33, l33, n33, p33 and Freq for the identified holder.

%% Input Variables
InputFreqMin = InputIRCSA.FreqMin;
InputFreqMax = InputIRCSA.FreqMax;
DeltaFreq = InputIRCSA.DeltaFreq;
signo = InputIRCSA.signIRCSA;

%% Main
% Read the minimum frequency of the experimental FRFs
MaxFreq_H11 = ExpFRF.h11(1,end);
MaxFreq_H12 = ExpFRF.h12(1,end);
MaxFreq_H33 = ExpFRF.h33(1,end);

MinFreq_H11 = ExpFRF.h11(1,1);
MinFreq_H12 = ExpFRF.h12(1,1);
MinFreq_H33 = ExpFRF.h33(1,1);

ExpFreqMax = min([MaxFreq_H11 MaxFreq_H12 MaxFreq_H33]);
ExpFreqMin = max([MinFreq_H11 MinFreq_H12 MinFreq_H33]);

% Frequency array definition for FRF operations
wint = InputFreqMin:DeltaFreq:InputFreqMax; % 0:DeltaFreq:InputFreqMax

if ExpFreqMax < InputFreqMax
   error = strjoin({'Reduce the maximum frequency for RCSA analysis to value', num2str(ExpFreqMax)});
   display(error)   
end

if ExpFreqMin < InputFreqMin
   error = strjoin({'Increase the minimum frequency for RCSA analysis to value', num2str(ExpFreqMin)});
   display(error)   
end

% IRCSAFreqMax = wint(end);

% Experimental FRF interpolation with the DeltaFreq defined in the funtion.
h11_exp = ((-1)^signo).*(interp1(ExpFRF.h11(1,:),real(ExpFRF.h11(2,:)),wint,'pchip')+j*interp1(ExpFRF.h11(1,:),imag(ExpFRF.h11(2,:)),wint,'pchip'));
h12_exp = ((-1)^signo).*(interp1(ExpFRF.h12(1,:),real(ExpFRF.h12(2,:)),wint,'pchip')+j*interp1(ExpFRF.h12(1,:),imag(ExpFRF.h12(2,:)),wint,'pchip'));
h33_exp = ((-1)^signo).*(interp1(ExpFRF.h33(1,:),real(ExpFRF.h33(2,:)),wint,'pchip')+j*interp1(ExpFRF.h33(1,:),imag(ExpFRF.h33(2,:)),wint,'pchip'));

% Simbolic operation to obtain the coefficients for H33 characterization
nn33 = zeros(length(wint),1);
pp33 = zeros(length(wint),1);
for idxFreq = 1:length(wint)
    u = h11_exp(1,idxFreq);
    v = h12_exp(1,idxFreq);
    a = AnalyticFRF.h11(idxFreq,1);
    b = AnalyticFRF.h12(idxFreq,1);
    c = AnalyticFRF.h21(idxFreq,1);
    d = AnalyticFRF.h22(idxFreq,1);
    e = AnalyticFRF.n12(idxFreq,1);
    f = AnalyticFRF.n21(idxFreq,1);
    g = AnalyticFRF.n22(idxFreq,1);
    k = h33_exp(1,idxFreq)+d;
    m = AnalyticFRF.p22(idxFreq,1);

    beta = (-k*f*v+f*b*k+k*g*u-k*g*a-d*b*f+g*c*b)/(c*b+d*u-d*a-v*c);
    delta = (k*g^2*a^2+k*f^2*b^2-u*f^2*d^2-2*u*k*f*v*g-f^2*b^2*d+f^2*d^2*a+k*f^2*v^2+g^2*k*u^2-f^2*c*b*d+f^2*v*d*b+...
            f*b^2*g*c-g^2*a*c*b+f^2*d*v*c+c^2*f*g*b-c^2*f*g*v-2*k*f^2*v*b-2*u*k*g^2*a+u*g^2*c*b-f*v*g*c*b-c*f*g*d*a+...
            g*a*d*b*f-2*k*f*b*g*a+2*k*f*v*g*a+2*u*g*f*k*b-u*g*d*b*f+u*c*f*g*d)/(c*b+d*u-d*a-v*c)^2;

    nn33(idxFreq,1) = beta-g;
    pp33(idxFreq,1) = delta-m;
end

Holder_IRCSA_FRF.h33 = h33_exp.';
Holder_IRCSA_FRF.n33 = nn33;
Holder_IRCSA_FRF.l33 = Holder_IRCSA_FRF.n33;
Holder_IRCSA_FRF.p33 = pp33;
Holder_IRCSA_FRF.Freq = wint; 

end % end of function
    
