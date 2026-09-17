% FUNCION PARA CALCULAR LOS PARÁMETROS MODALES DEL MODELO DINÁMICO
% Entradas:
%   - Matriz de masas
%   - Matriz de rigidez
% Salidas:
%   -Frecuencia natural
%   -Masa reflejada (En el cdg del agujero) y en la punta de hta.
%   -Vector normalizado (VN)
%____________________________________________________________________________

function [W,MMr_Tooltip,VN] = ModalParameters_Calculation_FRF(M,K,Nmodes)

%% Variable Information
% M (kg): mass matrix of the finite-element model.
% K (N/m): stiffness matrix of the finite-element model.
% Nmodes (-): modal-selection flag; 1 returns only the first mode, otherwise all modes.
% W (rad/s): natural angular frequencies ordered from lowest to highest.
% MMr_Tooltip (kg): reflected modal mass calculated from the first mode.
% VN (-): mass-normalized modal-vector matrix; it is 0 when Nmodes equals 1.

%% Main
%     [V,D] = eigs(K,M,15,'sm'); % Problema autovalores por mÃ©todos numÃ©ricos y Ãºnicamente considerando un nÃºmero determinado de modos
%     [V D] = eig(inv(M)*(K));% Solve eigenvalues and eigenvectors.
[V D] = eig(K,M);% Solve eigenvalues and eigenvectors.

% Put in order from lowest to highest the frequencies and the modes.
d = diag(D);
[d,ind] = sort(d);
D = diag(d);
V = V(:,ind);

[m,n] = size(M); % Obtain the number of degree of freedom
if Nmodes == 1
    W = D(1,1)^0.5;
    mm1 = (V(:,1))'*M*(V(:,1));
    MMr_Tooltip = mm1/((V(2,1))^2);
    VN = 0;
else
    % Obtain the natural frequency.
    for r = 1:m   % desde 1 al numero de frecuencias
        W(r) = D(r,r)^0.5;   % Frecuencias naturales en radianes/s
    end
    % Unity modal mass normalization
    MM = (V')*M*(V);
    for i = 1:m    % desde 1 al numero de modos
        VN(:,i) = V(:,i)/MM(i,i)^0.5;
    end
    % Reflected modal mass calcualtion. In the centre of the hole
    MMr_Tooltip = 1/((VN(2,1))^2);
end
end % end of function
