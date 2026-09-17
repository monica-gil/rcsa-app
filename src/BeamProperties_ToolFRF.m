function [ToolHeadInertia,ToolHeadLength,ToolHeadArea,ToolHeadKappa] = BeamProperties_ToolFRF(ToolData,HeadData)
%% Information
% CALCULO DE LAS PROPIEDADES FISICAS DEL MODELO BEAMS. BARRA ENTERIZA SIN EL PASIVO
% Parámetros obtenidos:
%   -Vector longitudes
%   -Vector areas
%   -Vector inercias

%% Variable Information
% ToolData: structure containing the tool geometry and material data.
%   * DiamExt (m): outer diameter of each tool section.
%   * DiamInt (m): inner diameter common to all tool sections.
%   * Length (m): length of each tool section.
%   * NElem (-): number of finite elements in each tool section.
%   * Nu (-): Poisson's ratio of the tool material.
% HeadData: structure containing the tool-head geometry and material data.
%   * Diam (m): outer diameter of the tool head.
%   * Length (m): total length of the tool head.
%   * CoGOffset (m): distance from the head end to its centre of gravity.
%   * NElem (-): number of finite elements in the tool head.
%   * Nu (-): Poisson's ratio of the tool-head material.
% ToolHeadInertia (m^4): second moment of area for each tool-head or tool element.
% ToolHeadLength (m): length of each tool-head or tool element.
% ToolHeadArea (m^2): cross-sectional area of each tool-head or tool element.
% ToolHeadKappa (-): shear correction factor of each tool-head or tool element.

%% Input Data
% Grinding tool
DiamExt = ToolData.DiamExt;
DiamInt = ToolData.DiamInt;
ToolLength = ToolData.Length;
NElem = ToolData.NElem;
ToolNu = ToolData.Nu;

% Grinding tool head
DiamHead = HeadData.Diam;
LengthHead = HeadData.Length;
CoGOffset = HeadData.CoGOffset;
NElemHead = HeadData.NElem;
HeadNu = HeadData.Nu;

%% Other Input Data
NSection = length(DiamExt);
% NElemTotal = sum(NElem) + NElemHead;

%% Main
% Grinding tool
ToolSection.Inertia = NaN(1,NSection);
ToolSection.Area = NaN(1,NSection);
ToolSection.Kappa = NaN(1,NSection);
for idxSection = 1:NSection
    ToolSection.Inertia(idxSection) = pi*(DiamExt(idxSection)^4-DiamInt^4)/64; % Inertia of the boring bar
    ToolSection.Area(idxSection) = pi*(DiamExt(idxSection)^2-DiamInt^2)/4; % Area of the boring bar
    ToolSection.Kappa(idxSection) = (6*(1+ToolNu)*(1+(DiamInt/DiamExt(idxSection))^2)^2)/((7+6*ToolNu)*(1+(DiamInt/DiamExt(idxSection))^2)^2+(20+12*ToolNu)*(DiamInt/DiamExt(idxSection))^2);
end

% Grinding head
HeadData.Inertia = pi*DiamHead^4/64; % Inertia of the beams employed in the toolhead
HeadData.Area = pi*(DiamHead^2)/4; % Area of the beams employed in the toolhead
HeadData.Kappa = (6*(1+HeadNu))/(7+6*HeadNu); % Shear coefficient of the beams employed in the toolhead

% Fill the vectors of the tool head
ToolAndHead.Length = repmat((LengthHead-CoGOffset)/(NElemHead-1),NElemHead,1);
ToolAndHead.Length(end,1) = CoGOffset;
ToolAndHead.Area = repmat(HeadData.Area,NElemHead,1);
ToolAndHead.Inertia = repmat(HeadData.Inertia,NElemHead,1);
ToolAndHead.Kappa = repmat(HeadData.Kappa,NElemHead,1);

for idxSection = NSection:-1:1
    ToolAndHead.Length = [ToolAndHead.Length; repmat(ToolLength(idxSection)/NElem(idxSection),NElem(idxSection),1)]; % !!
    ToolAndHead.Area = [ToolAndHead.Area; repmat(ToolSection.Area(idxSection),NElem(idxSection),1)];
    ToolAndHead.Inertia = [ToolAndHead.Inertia; repmat(ToolSection.Inertia(idxSection),NElem(idxSection),1)];
    ToolAndHead.Kappa = [ToolAndHead.Kappa; repmat(ToolSection.Kappa(idxSection),NElem(idxSection),1)];
end

ToolHeadLength = ToolAndHead.Length;
ToolHeadArea = ToolAndHead.Area;
ToolHeadInertia = ToolAndHead.Inertia;
ToolHeadKappa = ToolAndHead.Kappa;

end % end of function
