function LongBarModel = RCSABarDynamicProperties(LongBarData)
%% Information
% Calculate M, K, C matrices of IRCSA (long bar)

%% Variable Information
% LongBarData: structure containing the long-bar finite-element model data.
%   * DiamExt (m): outer diameter of each section.
%   * DiamInt (m): inner diameter common to all sections.
%   * Length (m): length of each section.
%   * NElem (-): number of elements in each section.
%   * NElemTotal (-): total number of finite elements.
%   * NDof (-): total number of degrees of freedom.
%   * Psi (%): modal damping ratio.
%   * E (Pa): Young's modulus.
%   * Rho (kg/m^3): material density.
%   * Nu (-): Poisson's ratio.
% LongBarModel: structure containing the free-free dynamic model of the long bar.
%   * M (kg): assembled mass matrix.
%   * K (N/m): assembled stiffness matrix.
%   * C (N*s/m): assembled damping matrix.

%% Input Variables
DiamExt = LongBarData.DiamExt;
DiamInt = LongBarData.DiamInt;
Len = LongBarData.Length;
NDof = LongBarData.NDof;
NElem = LongBarData.NElem;
NElemTotal = LongBarData.NElemTotal;
Psi = LongBarData.Psi;
E = LongBarData.E;
Rho = LongBarData.Rho;
Nu = LongBarData.Nu;

%% Main
[BarInertia,BarLength,BarArea,BarKappa] = BeamProperties_RCSA(DiamExt,DiamInt,Len,NElem,Nu);

LongBarModel.M = zeros(NDof,NDof);
LongBarModel.K = zeros(NDof,NDof);
for idxElem = 1:NElemTotal
    ElemArea = BarArea(idxElem);
    ElemInertia = BarInertia(idxElem);
    ElemLength = BarLength(idxElem);
    ElemKappa = BarKappa(idxElem);
    [MElement, KElement] = MassStiffnessMatrix_Timoshenko(E,Rho,Nu,ElemArea,ElemInertia,ElemLength,ElemKappa,0);
    
    idxDoFElem = [3*idxElem-2,3*idxElem-1,3*idxElem,3*idxElem+1,3*idxElem+2,3*idxElem+3];
    LongBarModel.M(idxDoFElem,idxDoFElem) = LongBarModel.M(idxDoFElem,idxDoFElem) + MElement;
    LongBarModel.K(idxDoFElem,idxDoFElem) = LongBarModel.K(idxDoFElem,idxDoFElem) + KElement;
end

% Damping matrix construction as proporcional damping
LongBarModel.C = DampingMatrixCalculationRealPart(LongBarModel.M,LongBarModel.K,Psi); % DampingMatrixCalculationRealPart(M,K,ToolData.Psi)

end % end of function
