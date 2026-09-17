function ToolModel = ToolDynamicProperties(ToolData,HeadData)
%% Information
% Calculate M, K, C matrices of the tool + head system
% The tool has a set of sections
% The head is modelled with a mass element

%% Variable Information
% ToolData: structure containing the tool finite-element model data.
%   * DiamExt (m): outer diameter of each tool section.
%   * DiamInt (m): inner diameter common to all tool sections.
%   * Length (m): length of each tool section.
%   * NElem (-): number of elements in each tool section.
%   * NElemTotal (-): total number of tool and head elements.
%   * NDOF (-): total number of degrees of freedom.
%   * Psi (%): modal damping ratio.
%   * E (Pa): Young's modulus.
%   * Rho (kg/m^3): material density.
%   * Nu (-): Poisson's ratio.
% HeadData: structure containing the tool-head finite-element model data.
%   * Mass (kg): concentrated mass of the tool head.
%   * Diam (m): outer diameter of the tool head.
%   * Length (m): total length of the tool head.
%   * CoGOffset (m): distance from the head end to its centre of gravity.
%   * NElem (-): number of head elements.
%   * E (Pa): Young's modulus.
%   * Rho (kg/m^3): material density.
%   * Nu (-): Poisson's ratio.
% ToolModel: structure containing the free-free dynamic model of the tool and head.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.

%% Input Data

%% Main
% Beam properties vector
[ToolHeadInertia,ToolHeadLength,ToolHeadArea,ToolHeadKappa] = BeamProperties_ToolFRF(ToolData,HeadData); % area, inertia, kappa and length for each element

% Initialize mass and stiffness matrix of the toolhead and cover plate FEM model
ToolModel.M = zeros(ToolData.NDOF,ToolData.NDOF);
ToolModel.K = zeros(ToolData.NDOF,ToolData.NDOF);

for idxElem = 1:ToolData.NElemTotal
    if idxElem <= HeadData.NElem
        E = HeadData.E; % tool head elements
        Rho = HeadData.Rho; % tool head elements
        Nu = HeadData.Nu; % tool head elements
    else
        E = ToolData.E; % tool elements
        Rho = ToolData.Rho; % tool elements
        Nu = ToolData.Nu; % tool elements
    end
    ElemArea = ToolHeadArea(idxElem);
    ElemInertia = ToolHeadInertia(idxElem);
    ElemLength = ToolHeadLength(idxElem);
    ElemKappa = ToolHeadKappa(idxElem);
    [MElement, KElement] = MassStiffnessMatrix_Timoshenko(E,Rho,Nu,ElemArea,ElemInertia,ElemLength,ElemKappa,0);

    idxDoFElem = [3*idxElem-2,3*idxElem-1,3*idxElem,3*idxElem+1,3*idxElem+2,3*idxElem+3];
    % idxDoFElem = 3*idxElem-2:3*idxElem+3;
    ToolModel.M(idxDoFElem,idxDoFElem) = ToolModel.M(idxDoFElem,idxDoFElem) + MElement;
    ToolModel.K(idxDoFElem,idxDoFElem) = ToolModel.K(idxDoFElem,idxDoFElem) + KElement;
end

% Add mass and inertia of the toolhead in its gravity center
ToolModel.M((HeadData.NElem*3-1),(HeadData.NElem*3-1)) = ToolModel.M((HeadData.NElem*3-1),(HeadData.NElem*3-1)) + HeadData.Mass;

% % C Matrix creation as proporcional damping
% Mintegral = M(1:(ToolData.NElemTotal*3),1:(ToolData.NElemTotal*3));
% Kintegral = K(1:(ToolData.NElemTotal*3),1:(ToolData.NElemTotal*3));
% 
% % Natural frequency and reflected modal mass calculation of the integral boring bar for TMD optimal tuning
% [ToolHead.FreqNat,ToolHead.MMr,ToolHead.VN] = ModalParameters_Calculation_FRF(Mintegral,Kintegral,0);
% NFreqNat = length(ToolHead.FreqNat);
% ToolHead.CMN = zeros(NFreqNat,NFreqNat);
% for idxFreq = 1:NFreqNat
%     ToolHead.CMN(idxFreq,idxFreq) = 2*(ToolData.Psi/100)*ToolHead.FreqNat(idxFreq);
% end
% % C Matrix is obtained operating with the normalized damping matrix (CMN):
% C = inv((ToolHead.VN)')*(ToolHead.CMN)*inv(ToolHead.VN);

% FRFs FREE-FREE FOR RCSA CALCULATION (free-free)
% Damping matrix of free-free model
ToolModel.C = DampingMatrixCalculationRealPart(ToolModel.M,ToolModel.K,ToolData.Psi); % DampingMatrixCalculation(M,K,ToolData.Psi)

end % end of function
