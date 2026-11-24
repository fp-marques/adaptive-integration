% integ_adapt
% Trabalho de grupo (3 alunos), LEBiom, IST, 2026/2027
% Calcula o integral de f com integracao adaptativa nao iterativa (algoritmo 2).
% Em cada um dos 3 trocos iniciais escolhe o numero de subintervalos a partir
% do majorante da derivada e aplica a regra do ponto a esquerda ou do trapezio.
%
% Entradas
%   f        funcao integranda, tem de aceitar vetores
%   xt       extremos dos 3 trocos iniciais, [a x1 x2 b]
%   tol      tolerancia pedida para o integral
%   regra    'esq' para o ponto a esquerda ou 'trap' para o trapezio
%   desenha  true para desenhar f e a discretizacao (opcional)
%   np       pontos de cada troco usados no majorante (opcional, 5 por omissao)
% Saidas
%   Ih       valor aproximado do integral
%   neval    numero de vezes que f foi calculada, incluindo as diferencas finitas
%   x        nos da discretizacao final
%   m        numero de subintervalos de cada troco

function [Ih,neval,x,m] = integ_adapt(f,xt,tol,regra,desenha,np)

if nargin<5, desenha=false; end
if nargin<6, np=5; end
a=xt(1); b=xt(end);
N=length(xt)-1;          % numero de trocos (N=3)
h=zeros(1,N); m=zeros(1,N); I=zeros(1,N);
neval=0;
for i=1:N
    h(i)=xt(i+1)-xt(i);
    ei=tol*h(i)/(b-a);                           % tolerancia do troco (9)
    if strcmp(regra,'esq')
        [D,nev]=majorante(f,xt(i),xt(i+1),1,np); % majorante de |f'|
        E=0.5*D*h(i)^2;                          % erro do troco sem subdividir
    else
        [D,nev]=majorante(f,xt(i),xt(i+1),2,np); % majorante de |f''|
        E=D*h(i)^3/12;
    end
    % como nas aulas: so se subdivide se o erro do troco passar a sua tolerancia
    if E<=ei
        m(i)=1;
    elseif strcmp(regra,'esq')
        m(i)=ceil(E/ei);                         % (11)
    else
        m(i)=ceil(sqrt(E/ei));                   % (10)
    end
    % regra composta no troco i, com passo constante em todo o troco
    xi=linspace(xt(i),xt(i+1),m(i)+1);
    hs=h(i)/m(i);
    if strcmp(regra,'esq')
        fi=f(xi(1:end-1));                       % o ultimo no nao entra (5)
        I(i)=hs*sum(fi);
    else
        fi=f(xi);
        I(i)=hs*(sum(fi)-(fi(1)+fi(end))/2);     % (6)
    end
    neval=neval+nev+length(fi);
    % nos da discretizacao, para devolver e desenhar (a fronteira entre trocos so entra uma vez)
    if i==1
        x=xi;
    else
        x=[x xi(2:end)];
    end
end
Ih=sum(I);               % Ih = I1 + I2 + I3

if desenha
    % eixo da esquerda: f e as fronteiras dos 3 trocos iniciais
    yyaxis left
    % 1000 pontos chegam para a curva parecer continua, porque unir os pontos por retas
    % afasta a curva de f no maximo cerca de 6e-6 nas tres funcoes, muito menos do que a
    % altura do grafico. Estes calculos de f so servem para o desenho e nao entram em neval
    xf=linspace(a,b,1000);
    plot(xf,f(xf),'b-','LineWidth',1.5); hold on
    ylabel('f(x)')
    yl=ylim;
    for i=2:N
        plot([xt(i) xt(i)],yl,'k--','LineWidth',1)
    end
    % eixo da direita: passo de cada troco, que mostra onde a malha e mais fina.
    % Os nos nao se marcam porque, mesmo no trapezio, as marcas ficariam sobrepostas
    yyaxis right
    hs=h./m;
    stairs(xt,[hs hs(end)],'-','LineWidth',1.5)
    set(gca,'YScale','log')   % escala logaritmica porque o passo vai de cerca de 1e-6 a 1e-2 entre as duas regras
    % limites em decadas inteiras, para os degraus nao ficarem colados as margens
    y1=10^floor(log10(min(hs))); y2=10^ceil(log10(max(hs)));
    if y1==y2, y2=10*y2; end
    ylim([y1 y2])
    ylabel('passo h_i/m_i')
    xlim([a b])
    xlabel('x'); grid on
    title(sprintf('m = [%d %d %d], %d calculos de f',m,neval))
end
end
