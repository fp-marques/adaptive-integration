% main_integ
% Trabalho de grupo (3 alunos), LEBiom, IST, 2026/2027
% Aplica a integracao adaptativa as tres funcoes de teste com as duas regras.
% Mostra a tabela para e = 1e-5 e a comparacao com o majorante do erro teorico
% (Tabela 2), e faz as figuras da discretizacao, do numero de calculos de f e do
% erro em funcao da tolerancia.

clear all
clc
close all
% funcoes de teste
fun={@(x) sin(x), @(x) exp(-x.^2), @(x) log(log(x))};
nome={'sin(x)','exp(-x^2)','ln(ln x)'};
% trocos iniciais [a x1 x2 b]: mais curtos onde as derivadas sao maiores
% ln(ln x) so existe para x>1 e as derivadas explodem perto de 1, por isso comeca em 1.5
XT={[0 pi/3 2*pi/3 pi], [0 0.5 1.2 2], [1.5 2 3 5]};
regras={'esq','trap'};
nomeregra={'Ponto a esquerda','Trapezio'};

% valores de referencia: sin tem primitiva conhecida, nas outras usa-se o integral do MATLAB com tolerancia 1e-13
Iref=zeros(1,3);
Iref(1)=-cos(pi)+cos(0);
for k=2:3
    Iref(k)=integral(fun{k},XT{k}(1),XT{k}(end),'RelTol',1e-13,'AbsTol',1e-13);
end

% tolerancias com meia decada entre cada uma. Abaixo de 1e-6 o ponto a esquerda
% precisa de dezenas de milhoes de nos, que ocupam centenas de MB (anexo C)
tols=10.^(-2:-0.5:-6);
nt=length(tols);
Ih=zeros(3,2,nt);         % valor do integral
neval=zeros(3,2,nt);      % numero de calculos de f
erro=zeros(3,2,nt);       % erro efetivo |I-Ih|
for k=1:3
    for r=1:2
        for it=1:nt
            [Ih(k,r,it),neval(k,r,it)]=integ_adapt(fun{k},XT{k},tols(it),regras{r});
            erro(k,r,it)=abs(Iref(k)-Ih(k,r,it));
        end
    end
end

% tabela para e = 1e-5
tol5=1e-5;
fprintf('e = %.0e\n',tol5)
fprintf('%-10s %-17s %18s %18s %10s %8s   %s\n','funcao','regra','Ih','I referencia','|I-Ih|','cumpre','m e calculos de f')
for k=1:3
    for r=1:2
        [I,ne,~,m]=integ_adapt(fun{k},XT{k},tol5,regras{r});
        E=abs(Iref(k)-I);
        if E<tol5, ok='sim'; else, ok='nao'; end
        fprintf('%-10s %-17s %18.12f %18.12f %10.2e %8s   m=[%d %d %d], %d\n', ...
            nome{k},nomeregra{r},I,Iref(k),E,ok,m,ne)
    end
end

% Tabela 2: comparacao do erro efetivo com o majorante do erro teorico Eteor, que e a
% soma dos erros estimados dos 3 trocos depois da subdivisao, com os D_i do majorante.m
fprintf('\nErro efetivo e majorante do erro teorico\n')
fprintf('%-10s %-17s %8s %10s %14s\n','funcao','regra','e','Eteor/e','|I-Ih|/Eteor')
for k=1:3
    xt=XT{k};
    for r=1:2
        for it=1:nt
            [Iht,~,~,m]=integ_adapt(fun{k},xt,tols(it),regras{r});
            Eteor=0;
            for i=1:3
                h=xt(i+1)-xt(i);
                if r==1
                    D=majorante(fun{k},xt(i),xt(i+1),1);
                    Eteor=Eteor+0.5*D*h^2/m(i);      % erro do troco dividido em m(i) subintervalos
                else
                    D=majorante(fun{k},xt(i),xt(i+1),2);
                    Eteor=Eteor+D*h^3/12/m(i)^2;     % erro do troco dividido em m(i) subintervalos
                end
            end
            fprintf('%-10s %-17s %8.0e %10.3f %14.3g\n',nome{k},nomeregra{r},tols(it), ...
                Eteor/tols(it),abs(Iref(k)-Iht)/Eteor)
        end
    end
end

% figura 1: f e a discretizacao para e = 1e-5 (a pedida no enunciado)
fig1=figure('Position',[50 50 1200 900],'DefaultAxesFontSize',10);
for k=1:3
    for r=1:2
        subplot(3,2,2*(k-1)+r)
        integ_adapt(fun{k},XT{k},tol5,regras{r},true);
        letra=char('a'+2*(k-1)+r-1);   % letras (a) a (f), cada linha e uma funcao e cada coluna uma regra
        title({['(' letra ') ' nome{k} ' - ' nomeregra{r}], get(get(gca,'Title'),'String')})
    end
end

% figura 2: numero de calculos de f em funcao da tolerancia
% Os eixos estao em escala logaritmica porque a tolerancia, os calculos de f e o erro variam
% ao longo de varias ordens de grandeza. Assim cabem todos no mesmo grafico e as relacoes
% previstas por (10) e (11) aparecem como retas, de declive -1 no ponto a esquerda e -0.5 no
% trapezio, que se leem diretamente no grafico
cor={[235 104 52]/255, [42 120 214]/255};
fig2=figure('Position',[50 50 1100 350],'DefaultAxesFontSize',10);
for k=1:3
    subplot(1,3,k)
    for r=1:2
        loglog(tols,squeeze(neval(k,r,:)),'o-','Color',cor{r},'LineWidth',1.5); hold on
    end
    set(gca,'XDir','reverse')     % da tolerancia menos exigente para a mais exigente
    xlabel('tolerancia e'); ylabel('calculos de f'); title(['(' char('a'+k-1) ') ' nome{k}]); grid on
    legend(nomeregra,'Location','northwest')
end

% figura 3: erro efetivo em funcao da tolerancia (abaixo da reta tracejada cumpre a tolerancia),
% em escala logaritmica pelo mesmo motivo da figura 2
fig3=figure('Position',[50 50 1100 350],'DefaultAxesFontSize',10);
for k=1:3
    subplot(1,3,k)
    for r=1:2
        loglog(tols,squeeze(erro(k,r,:)),'o-','Color',cor{r},'LineWidth',1.5); hold on
    end
    loglog(tols,tols,'k--')
    set(gca,'XDir','reverse')
    xlabel('tolerancia e'); ylabel('|I - I_h|'); title(['(' char('a'+k-1) ') ' nome{k}]); grid on
    legend([nomeregra {'|I-I_h| = e'}],'Location','southwest')
end

% fundo branco para o relatorio: o tema claro do MATLAB R2025 deixa o fundo da figura cinzento,
% por isso depois do tema poe-se a cor da figura a branco
nomesfig={'fig1_discretizacao.png','fig2_calculos.png','fig3_erro.png'};
figs=[fig1 fig2 fig3];
for i=1:3
    try
        theme(figs(i),'light')
    catch
    end
    set(figs(i),'Color','w')
    exportgraphics(figs(i),nomesfig{i},'Resolution',300,'BackgroundColor','white')
end
