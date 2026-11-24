% majorante
% Trabalho de grupo (3 alunos), LEBiom, IST, 2026/2027
% Estima o majorante de |f'| ou de |f''| num troco [a,b] (algoritmo 1), com
% diferencas finitas centrais calculadas em np pontos do troco.
%
% Entradas
%   f     funcao, tem de aceitar vetores
%   a, b  extremos do troco
%   k     ordem da derivada, 1 ou 2
%   np    numero de pontos do troco onde se calcula a diferenca finita (5 por omissao)
% Saidas
%   D     maior valor absoluto da diferenca finita nos np pontos
%   nev   numero de vezes que f foi calculada

function [D,nev] = majorante(f,a,b,k,np)

if nargin<5, np=5; end
u=eps/2;                       % unidade de arredondamento do MATLAB, cerca de 1.1e-16
if np==1
    z=(a+b)/2;                 % so o ponto medio (linspace(a,b,1) daria b)
else
    z=linspace(a,b,np);        % np pontos igualmente espacados, com os dois extremos
end
if k==1
    d=(3*u)^(1/3);             % passo otimo (15) com M3=1, cerca de 6.9e-6
    Dz=(f(z+d)-f(z-d))/(2*d);              % diferenca central (12)
    nev=2*length(z);
else
    d=(48*u)^(1/4);            % passo otimo para a segunda derivada com M4=1, cerca de 2.7e-4
    Dz=(f(z+d)-2*f(z)+f(z-d))/d^2;         % diferenca central (13)
    nev=3*length(z);
end
D=max(abs(Dz));                % os np valores sao calculados de uma vez (vetorial)
end
