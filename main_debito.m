% main_debito
% Trabalho de grupo (3 alunos), LEBiom, IST, 2026/2027
% Calcula a area sob a curva de termodiluicao com as regras do ponto a esquerda
% e do trapezio, o debito cardiaco pela formula de Stewart Hamilton e o erro
% relativo de cada regra. Faz tambem a figura da curva com as duas aproximacoes.

clear all
clc
% dados da termodiluicao (enunciado)
t=[0 2 4 6 8 10 12 14 16];                   % tempo (s)
dT=[0.0 0.1 0.4 0.9 1.2 0.8 0.3 0.1 0.0];    % variacao de temperatura (C)
Aref=7.25;                                   % area de referencia (C.s)
h=t(2)-t(1);                                 % espacamento constante, h=2 s

% regras compostas (5) e (6) com os 9 pontos da tabela (N=8)
Aesq=h*sum(dT(1:end-1));                     % o ultimo ponto nao entra
Atrap=(sum(dT)-(dT(1)+dT(end))/2)*h;         % como no trap.m das aulas
A=[Aesq Atrap];

Q=100./A;                                    % Stewart-Hamilton (16)
Er=abs(Aref-A)/Aref*100;                     % erro relativo percentual (17)

nomeregra={'Ponto a esquerda','Trapezio'};
fprintf('%-17s %10s %12s %10s\n','regra','area','Q (L/min)','Er (%)')
for r=1:2
    fprintf('%-17s %10.4f %12.4f %10.2f\n',nomeregra{r},A(r),Q(r),Er(r))
end
fprintf('%-17s %10.4f %12.4f\n','referencia',Aref,100/Aref)

% curva e as duas aproximacoes
fig=figure('DefaultAxesFontSize',11);
stairs(t,dT,'-','Color',[235 104 52]/255,'LineWidth',1.5); hold on   % ponto a esquerda: degraus
plot(t,dT,'-','Color',[42 120 214]/255,'LineWidth',1.5)              % trapezio: segmentos de reta
plot(t,dT,'ko','MarkerFaceColor','k')
xlabel('t (s)'); ylabel('\DeltaT (°C)'); grid on
legend('Ponto a esquerda','Trapezio','Dados','Location','northeast')
try
    theme(fig,'light')
catch
end
set(fig,'Color','w')   % fundo branco para o relatorio
exportgraphics(fig,'fig4_debito.png','Resolution',300,'BackgroundColor','white')
