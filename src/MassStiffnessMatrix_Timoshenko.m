function [Mgeneral, Kgeneral] = MassStiffnessMatrix_Timoshenko(E,ro,nu,A,I,L,kappa,theta)
%% Information
% CALCULO DE LA MATRIZ DE MASAS Y RIGIDEZ EN EL SISTEMA DE COORDENADAS
% GENERAL A PARTIR DE DATOS DE ENTRADA DE PROPIEDADES DE LOS BEAM
% Entradas:
%   - Propiedades de material (Modulo elástico y densidad)
%   - Propiedades del elemento beam (Area, inercia y longitud)
%   - Angulo de giro para paso de coordenadas locales a generales

%% Variable Information
% E (Pa): Young's modulus of the element material.
% ro (kg/m^3): density of the element material.
% nu (-): Poisson's ratio of the element material.
% A (m^2): cross-sectional area of the beam element.
% I (m^4): second moment of area of the beam element.
% L (m): length of the beam element.
% kappa (-): shear correction factor of the beam element.
% theta (rad): rotation from local to global coordinates.
% Mgeneral (kg): 6-by-6 element mass matrix in global coordinates.
% Kgeneral (N/m): 6-by-6 element stiffness matrix in global coordinates.

%% Input Data

%% Main
phi = (24*(1/kappa)*(1+nu)*I)/(A*L^2);

% Matrices locales
%     Mlocal = (ro*A*L/840)*[280             0                          0               140             0                            0;
%                           0     312+588*phi+280*phi^2     (44+77*phi+35*phi^2)*L     0      108+252*phi+175*phi^2     (-26-63*phi-35*phi^2)*L;
%                           0    (44+77*phi+35*phi^2)*L     (8+14*phi+7*phi^2)*L^2     0      (26+63*phi+35*phi^2)*L    (-6-14*phi-7*phi^2)*L^2;
%                          140             0                          0               280               0                           0;
%                           0     108+252*phi+175*phi^2     (26+63*phi+35*phi^2)*L     0     (312+588*phi+280*phi^2)     (-44-77*phi-35*phi^2)*L;
%                           0    (-26-63*phi-35*phi^2)*L    (-6-14*phi-7*phi^2)*L^2    0     (-44-77*phi-35*phi^2)*L     (8+14*phi+7*phi^2)*L^2];

Mlocal_1 = (ro*A*L)/(210*(1+phi)^2)*[70*(1+phi)^2                     0                           0                  35*(1+phi)^2                     0                           0;
    0                   (70*phi^2+147*phi+78)     (35*phi^2+77*phi+44)*L/4          0                  (35*phi^2+63*phi+27)        -(35*phi^2+63*phi+26)*L/4;
    0                (35*phi^2+77*phi+44)*L/4     (7*phi^2+14*phi+8)*L^2/4          0               (35*phi^2+63*phi+26)*L/4       -(7*phi^2+14*phi+6)*L^2/4;
    35*(1+phi)^2                        0                           0                  70*(1+phi)^2                     0                           0;
    0                  (35*phi^2+63*phi+27)        (35*phi^2+63*phi+26)*L/4         0                   (70*phi^2+147*phi+78)      -(35*phi^2+77*phi+44)*L/4;
    0               -(35*phi^2+63*phi+26)*L/4      -(7*phi^2+14*phi+6)*L^2/4       0                -(35*phi^2+77*phi+44)*L/4      (7*phi^2+14*phi+8)*L^2/4];



Mlocal_2 = (ro*I/(30*(1+phi)^2*L))*[     0                             0                           0                     0                          0                           0;
    0                            36                     -(15*phi-3)*L               0                         -36                   -(15*phi-3)*L;
    0                      -(15*phi-3)*L           (10*phi^2+5*phi+4)*L^2           0                      (15*phi-3)*L       (5*phi^2-5*phi-1)*L^2;
    0                             0                           0                     0                          0                           0;
    0                            -36                     (15*phi-3)*L               0                          36                   (15*phi-3)*L;
    0                      -(15*phi-3)*L           (5*phi^2-5*phi-1)*L^2            0                      (15*phi-3)*L       (10*phi^2+5*phi+4)*L^2];


Mlocal = Mlocal_1 + Mlocal_2;

Klocal = (E*I/L^3)*[A*L^2/I         0                0           -A*L^2/I      0                    0;
    0       12/(1+phi)        6*L/(1+phi)        0       -12/(1+phi)         6*L/(1+phi);
    0       6*L/(1+phi)   (4+phi)*L^2/(1+phi)    0       -6*L/(1+phi)     (2-phi)*L^2/(1+phi);
    -A*L^2/I         0                0            A*L^2/I       0                   0;
    0       -12/(1+phi)      -6*L/(1+phi)        0       12/(1+phi)         -6*L/(1+phi);
    0       6*L/(1+phi)   (2-phi)*L^2/(1+phi)    0      -6*L/(1+phi)      (4+phi)*L^2/(1+phi)];


% Cambio de coordenada de locales a generales
R = [cos(theta)  sin(theta)   0       0           0        0;
    -sin(theta)  cos(theta)   0       0           0        0;
    0           0         1       0           0        0;
    0           0         0   cos(theta)  sin(theta)   0;
    0           0         0   -sin(theta) cos(theta)   0;
    0           0         0       0           0        1];

% Matriz Masa y Rigidez tras el cambio de coordenadas
Mgeneral = R'*Mlocal*R;
Kgeneral = R'*Klocal*R;

end % end of function
