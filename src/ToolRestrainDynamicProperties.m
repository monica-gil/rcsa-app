function ToolModelRestrained = ToolRestrainDynamicProperties(ToolModelFree,ToolData,NumDoFNode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input Data
FreeM = ToolModelFree.M;
FreeK = ToolModelFree.K;
Psi = ToolData.Psi;

%% Main
[NDoF,~] = size(FreeM);
NNode = NDoF/NumDoFNode;
NElemTotal = NNode - 1;
NRestrain = NNode*NumDoFNode;

RestrainM = FreeM;
RestrainK = FreeK;

% C Matrix creation as proporcional damping
Mintegral = RestrainM(1:(NElemTotal*NumDoFNode),1:(NElemTotal*NumDoFNode));
Kintegral = RestrainK(1:(NElemTotal*NumDoFNode),1:(NElemTotal*NumDoFNode));

% Natural frequency and reflected modal mass calculation of the integral boring bar for TMD optimal tuning
% RestrainC = DampingMatrixCalculation(Mintegral,Kintegral,Psi);
RestrainC = DampingMatrixCalculationRealPart(Mintegral,Kintegral,Psi); %!!

RestrainM(NRestrain-2:1:NRestrain,:) = [];
RestrainM(:,NRestrain-2:1:NRestrain) = [];

RestrainK(NRestrain-2:1:NRestrain,:) = [];
RestrainK(:,NRestrain-2:1:NRestrain) = [];

%% Output Struct
ToolModelRestrained.M = RestrainM;
ToolModelRestrained.K = RestrainK;
ToolModelRestrained.C = RestrainC;

end % end of function
