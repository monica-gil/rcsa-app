function [BarInertia,BarLength,BarArea,BarKappa] = BeamProperties_RCSA(DiamExt,DiamInt,Len,NElem,Nu)

%% Information
% CALCULO DE LAS PROPIEDADES FISICAS DEL MODELO BEAMS. BARRA ENTERIZA SIN EL PASIVO
% Parámetros obtenidos:
%   -Vector longitudes
%   -Vector areas
%   -Vector Inercias

%% Variable Information
% DiamExt (m): outer diameter of each LongBarData section.
% DiamInt (m): inner diameter common to all LongBarData sections.
% Len (m): length of each LongBarData section.
% NElem (-): number of finite elements in each LongBarData section.
% Nu (-): Poisson's ratio of the LongBarData material.
% BarInertia (m^4): second moment of area for each finite element.
% BarLength (m): length of each finite element.
% BarArea (m^2): cross-sectional area of each finite element.
% BarKappa (-): shear correction factor of each finite element.

%% Input Data

%% Other Input Data
NSection = length(DiamExt); % number of sections (diam, length)

%% Main
% Properties of each section
BarSection.Inertia = NaN(1,NSection);
BarSection.Area = NaN(1,NSection);
BarSection.Kappa = NaN(1,NSection);
for idxSection = 1:NSection
    BarSection.Inertia(idxSection) = pi*(DiamExt(idxSection)^4-DiamInt^4)/64; % Inertia of the boring bar
    BarSection.Area(idxSection) = pi*(DiamExt(idxSection)^2-DiamInt^2)/4; % Area of the boring bar
    BarSection.Kappa(idxSection) = (6*(1+Nu)*(1+(DiamInt/DiamExt(idxSection))^2)^2)/((7+6*Nu)*(1+(DiamInt/DiamExt(idxSection))^2)^2+(20+12*Nu)*(DiamInt/DiamExt(idxSection))^2);
end

% Properties of each element
BarLength =  [];
BarArea = [];
BarInertia = [];
BarKappa = [];
for idxSection = NSection:-1:1
    BarLength = [BarLength; repmat(Len(idxSection)/NElem(idxSection),NElem(idxSection),1)];
    BarArea = [BarArea; repmat(BarSection.Area(idxSection),NElem(idxSection),1)];
    BarInertia = [BarInertia; repmat(BarSection.Inertia(idxSection),NElem(idxSection),1)];
    BarKappa = [BarKappa; repmat(BarSection.Kappa(idxSection),NElem(idxSection),1)];
end
end % end of function
