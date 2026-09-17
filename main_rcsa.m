%% Receptance Coupling Substructure Analysis (RCSA)

% 1.0 initial version (Josu Peña)
% 2.0 19/05/2026 (Monica Gil)
% 3.0 09/09/2026 chatGPT (Monica Gil)

%% Information
% RCSA calculation for a beam with 3 sections
% Fixed connection between the tool and the tool holder

%% Variable Information
% Steel
%   * E: Young modulus (Pa)
%   * Rho: density (kg/m^3)
%   * Nu: Poisson coefficient (-)
% LongBarData: long bar used for IRCSA
%   * Length (mm): length of each section (array)
%   * DiamExt (mm): outer diameter of each section (array)
%   * NElem (-): number of elements in each section (array)
%   * DiamInt (mm): inner diameter of each section (unique value for all sections)
%   * Psi: relative damping
%   * E: Young modulus (Pa)
%   * Rho: density (kg/m^3)
%   * Nu: Poisson coefficient (-)
% Filename: name of the MAT file containing the experimental FRFs
% ToolData: grinding tool or boring bar where RCSA is applied
%   * Length (mm): length of each section (array)
%   * DiamExt (mm): outer diameter of each section (array)
%   * NElem (-): number of elements in each section (array)
%   * DiamInt (mm): inner diameter of each section (unique value for all sections)
%   * Psi: relative damping
%   * E: Young modulus (Pa)
%   * Rho: density (kg/m^3)
%   * Nu: Poisson coefficient (-)
% HeadData: head of the tool modelled as COMN2 element (mass element)
%   * Mass (kg)
%   * Length (mm)
%   * CoGOffset (mm): distance from the head end to its centre of gravity.
%   * Diam (mm)
%   * NElem (-)
%   * E: Young modulus (Pa)
%   * Rho: density (kg/m^3)
%   * Nu: Poisson coefficient (-)
% LongBarModel: dynamic model of the free-free LongBarData.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.
% ToolModelFree: dynamic model of the free-free ToolData and HeadData.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.
% ToolModelRestrained: dynamic model of the ToolData and HeadData with the restraint applied.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.
% InputIRCSA:
%   * FreqMin (Hz)
%   * FreqMax (Hz)
%   * DeltaFreq (Hz)
%   * signIRCSA (-): sign convention; 0 preserves the experimental FRF sign and 1 reverses it.
% EvaluateIRCSA: evaluation points for the LongBarData free-free FRFs.
%   * NodeP1: first evaluation node.
%   * DoFP1: degree of freedom of the first evaluation node.
%   * idxP1: index of the first evaluation degree of freedom.
%   * NodeP2: second evaluation node.
%   * DoFP2: degree of freedom of the second evaluation node.
%   * idxP2: index of the second evaluation degree of freedom.
% InputRCSA:
%   * NodeP1: node P1 at which the FRF is calculated
%   * DoFP1: dof of the node P1 where the FRF is calculated
%   * idxP1: index of P1 where the FRF is calculated
%   * NodeP2: node P2 at which the FRF is calculated
%   * DoFP2: dof of the node P2 where the FRF is calculated
%   * idxP2: index of P2 where the FRF is calculated
% InputToolRestrain (FRF adquisition):
%   * FreqMin (Hz)
%   * FreqMax (Hz)
%   * DeltaFreq (Hz)
% EvaluateToolRestrain: evaluation points for the restrained-tool FRFs.
%   * Node: (array) node of each point at which the FRF is calculated
%   * DoF: (array) dof of each point where the FRF is calculated
%   * idxPoint: (array) index of each point where the FRF is calculated
% Result:
%   * FreqIni (Hz): initial frequency used to search for the RCSA peak
%   * BandWidth (Hz): frequency bandwidth used to identify the RCSA peak
% NumDoFNode: number of dof of each node

% ToolData
% HeadData
% ToolHead -> ToolData + HeadData

%% Clear
clear variables
% close all

%% Path
addpath('src')
addpath('data')

%% Input Data | General
% Materials
Steel.E = 2.068*10^11; % Young's modulus of the steel (Pa)
Steel.Rho = 7820; % Density of the steel (kg/m^3)
Steel.Nu = 0.29; % (-)

%% Input Parameters | IRCSA
% Folder = '';
Filename = 'FRFs_IRCS_DIRY.mat';

% LongBarData section: [section 1,  section 2, ...]
LongBarData.Length = 405*1e-3; % [268 16.5 16.5]*1e-03; % (m) % Length of the bar employed for IRCSA
LongBarData.DiamExt = 62.8*1e-3; % [58 31.75 12]*10^-3; % (m) % Diameter of the bar employed for IRCSA
LongBarData.DiamInt = 18*1e-3; % 5*10^-3; % (m)
LongBarData.NElem = 20;
% RCSA information
LongBarData.Psi = 0.1;
% Material
LongBarData.E = 2.068*10^11; %(Pa)
LongBarData.Rho = 7820; %(kg/m^3)
LongBarData.Nu = 0.29;

% Frequency range for IRCSA
InputIRCSA.FreqMin = 0; % (Hz)
InputIRCSA.FreqMax = 1600; % 2000; % (Hz)
InputIRCSA.DeltaFreq = 0.1; % (Hz) precision
InputIRCSA.signIRCSA = 0; % Sign 0->Positive ; 1->Negative

% Evaluation points for the LongBarData free-free FRFs
EvaluateIRCSA.NodeP1 = 1; % 3rd dof of the union node of section 1 and 2
EvaluateIRCSA.DoFP1 = 3; % 3rd dof of NodeP1 
EvaluateIRCSA.NodeP2 = sum(LongBarData.NElem) + 1; % 3rd dof of the last node
EvaluateIRCSA.DoFP2 = 3; % 3rd dof of NodeP2 

% Frequency range for the FRF
InputToolRestrain.FreqMin = 0; % (Hz)
InputToolRestrain.FreqMax = 400; % (Hz)
InputToolRestrain.DeltaFreq = 0.1; % (Hz) precision

Result.FreqIni = 20; % (Hz)
Result.BandWidth = 20; % (Hz)

%% Input Data | ToolData
% ToolData sections: [section1 section2, ....]
ToolData.Length = 405*1e-03; % [105.5 12.5 20]; [93 10 20]*1e-03; % (m)
ToolData.DiamExt = 62.8*1e-03; %  [63 41.5 20]; [63 45 30]*1e-03; % (m)
ToolData.DiamInt = 18e-03; % (m) % internal hole
ToolData.NElem = 32;

% Bar damping
ToolData.Psi = 1; % 0.1; % Relative damping percentage of the integral boring bar
% Material
ToolData.Rho = Steel.Rho;
ToolData.E = Steel.E;
ToolData.Nu = Steel.Nu;

% Grinding wheel mass
HeadData.Mass = 7.118; % 0.05; 0.11(kg)
HeadData.Length = 2*21.8*1e-3; % 30; 20(m)
HeadData.CoGOffset = 21.8*1e-3; % HeadData.Length/2; % 10*1e-3; %% (m)
HeadData.Diam = ToolData.DiamExt; % (m)
HeadData.NElem = 2; % Fix value. Don't modify
% HeadData mass and properties
HeadData.Rho = 7820*10^-3; % (kg/m^3)
HeadData.E = 2.2*10^11; % (Pa)
HeadData.Nu = 0.29; % (-)

% Evaluation points for the ToolData free-free FRFs 
InputRCSA.NodeP1 = 3; % 2; % 2nd node
InputRCSA.DoFP1 = 3; % 3rd dof
InputRCSA.NodeP2 = sum(ToolData.NElem) + HeadData.NElem + 1; % sum(ToolData.NElem) + HeadData.NElem + 1
InputRCSA.DoFP2 = 3; % 3rd dof

% Evaluation points for the restrained-tool FRFs
EvaluateToolRestrain.Node = InputRCSA.NodeP1; % EvaluateIRCSA.NodeP1; [3 0]
EvaluateToolRestrain.DoF = InputRCSA.DoFP1; % [P1, P2, ...] [3 3];

NumDoFNode = 3; % (fixed) number of degree of freedom in each node

%% Input Results
% figure axes limits
FigOpts.Xlim = [30 2000];
FigOpts.Ylim = [0 8e-7];

%% Calculate Other Data
% IRCSA | Long bar
LongBarData.NElemTotal = sum(LongBarData.NElem); 
LongBarData.Node = LongBarData.NElemTotal + 1;
LongBarData.NDof = NumDoFNode*LongBarData.Node;

% Evaluation points for the LongBarData free-free FRFs
EvaluateIRCSA.idxP1 = EvaluateIRCSA.NodeP1*EvaluateIRCSA.DoFP1; 
EvaluateIRCSA.idxP2 = EvaluateIRCSA.NodeP2*EvaluateIRCSA.DoFP2; 

% RCSA | ToolData 
ToolData.NElemTotal = sum(ToolData.NElem) + HeadData.NElem; % total number of elements (tool + head)
ToolData.NNode = ToolData.NElemTotal + 1;
ToolData.NDOF = ToolData.NNode*NumDoFNode; % number of dgf

% Evaluation points for the ToolData restrained FRFs
InputRCSA.idxP1 = InputRCSA.NodeP1*InputRCSA.DoFP1; % node 2; 3rd dof
InputRCSA.idxP2 = InputRCSA.NodeP2*InputRCSA.DoFP2;

% Evaluation points for the restrained-tool FRFs
EvaluateToolRestrain.idxPoint = EvaluateToolRestrain.Node.*EvaluateToolRestrain.DoF; 

%% EXPERIMENTAL (SHORT AND LONG BARS)
% Load experimental FRFs 
% addpath(genpath('src\IRCSA_Experimental'));
load(Filename)
ExpFRF.h11 = Dummy.H11;
ExpFRF.h12 = Dummy.H12;
ExpFRF.h33 = Dummy.H33;
clear H11 H12 H33;

%% FREE-FREE BEAM MODEL for IRCSA
% mass, stiffnes and damping matrixes
LongBarModel = RCSABarDynamicProperties(LongBarData);

% free-free analytic FRF IRCSA
LongBarAnalyticFRF = Calculate_FRF_RCSA(LongBarModel,InputIRCSA,EvaluateIRCSA);

%% IRCSA
% h33,n33,l33 and p33 FRF adquisition.
HolderIRCSAFRF = CalculateIRCSA(ExpFRF,LongBarAnalyticFRF,InputIRCSA);

%% FREE-FREE BEAM MODEL FOR THE TOOL
% mass, stiffnes and damping matrixes
ToolModelFree = ToolDynamicProperties(ToolData,HeadData);

% free-free analytic FRF 
ToolAnalyticFRF = Calculate_FRF_RCSA(ToolModelFree,InputIRCSA,InputRCSA);

%% RCSA
% RCSA calculation by H33 experimental and H11, H12, H22 analytical set of FRFs
[ToolCoupledRCSAFRF] = RCSA_Response(ToolAnalyticFRF,HolderIRCSAFRF,InputIRCSA);

%% GRINDING TOOL WITH RESTRAINS
% mass, stiffnes and damping matrixes
ToolModelRestrained = ToolRestrainDynamicProperties(ToolModelFree,ToolData,NumDoFNode);

% FRF in the tooltip and also the initial zone of the bar (experimental).
[ToolRestrainedFRF] = Calculate_FRF(ToolModelRestrained,InputToolRestrain,EvaluateToolRestrain);
[Kst] = StaticStiffness_Calculation(ToolModelRestrained.K);

%% RESULTS. PLOT WITH THE FREQUENCY RESPONSE FUNTION OF THE BORING BAR.
% display experimental FRF for IRCSA and calculated FRF of the tool

DisplayFigureFRF(ExpFRF,ToolCoupledRCSAFRF,ToolRestrainedFRF,FigOpts)

%%
% % RESULTS
% % Restrained model
% [AmplRestain,idxPeaksRestain] = findpeaks(ToolRestrainedFRF.Magnitude(:,1));
% Result.FreqRestain = (ToolRestrainedFRF.Freq(idxPeaksRestain));
% Result.AmplRestain = AmplRestain;
% 
% % RCSA model
% FreqIni = Result.FreqIni;
% DeltaFreqIRCSA = InputIRCSA.DeltaFreq;
% idxFreqIni = FreqIni/DeltaFreqIRCSA;
% [valueMax, idxMax] = max(abs(ToolCoupledRCSAFRF.G11(idxFreqIni:end)));
% idxMax = idxMax + idxFreqIni; % update idx value for complete frequency range
% 
% BandWidth = Result.BandWidth;
% deltaIdxBandWidth = BandWidth/DeltaFreqIRCSA;
% idxFreqIni2 = idxMax - deltaIdxBandWidth;
% idxFreqEnd2 = idxMax + deltaIdxBandWidth;
% [AmplRCSA,idxMaxPeaks] = findpeaks(abs(ToolCoupledRCSAFRF.G11(idxFreqIni2:idxFreqEnd2)),'MinPeakProminence', 0.15*valueMax);
% idxMaxPeaks = idxFreqIni2 + idxMaxPeaks - 1; % update idx value for complete frequency range
% Result.AmplRCSA = max(AmplRCSA);
% Result.FreqRCSA = abs(ToolCoupledRCSAFRF.Freq(idxMaxPeaks(1)));

