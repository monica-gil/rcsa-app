%FUNCION PARA CALCULAR LA RIGIDEZ ESTÁTICA EN LA PUNTA DE LA HERRAMIENTA
% Entradas:
%   - Matriz de rigidez (K)
% Salidas:
%   -Rigidez estática en la punta(Kst)
%____________________________________________________________________________

function [Kst]=StaticStiffness_Calculation(K)
    %% Variable Information
    % K (N/m): stiffness matrix of the restrained finite-element model.
    % Kst (N/m): static stiffness evaluated at the second degree of freedom.

    [a,b]=size(K);
    F=zeros(1,a);
    F(1,2)=1000;
    X=F*inv(K);
    Kst=F(1,2)/X(1,2);
    
end
