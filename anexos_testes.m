% anexos_testes
% Trabalho de grupo (3 alunos), LEBiom, IST, 2026/2027
% Testes dos anexos do relatorio. O anexo A ve o efeito do numero de pontos do
% majorante, o anexo B o erro das diferencas finitas em funcao do passo e o
% anexo C o tempo e a memoria para tolerancias menores.

clear all
clc
close all
fun={@(x) sin(x), @(x) exp(-x.^2), @(x) log(log(x))};
nome={'sin(x)','exp(-x^2)','ln(ln x)'};
XT={[0 pi/3 2*pi/3 pi], [0 0.5 1.2 2], [1.5 2 3 5]};
regras={'esq','trap'};
nomeregra={'Ponto a esquerda','Trapezio'};
Iref=zeros(1,3);
Iref(1)=2;
for k=2:3
    Iref(k)=integral(fun{k},XT{k}(1),XT{k}(end),'RelTol',1e-13,'AbsTol',1e-13);
end
cor={[235 104 52]/255, [42 120 214]/255};
marc={'o','s','^'};
corf={[235 104 52]/255, [42 120 214]/255, [27 175 122]/255};   % uma cor por funcao (anexo B)

%% Anexo A: numero de pontos usados no majorante, com e = 1e-5
tol=1e-5;
nps=[1 2 3 5 9 17];      % com 1 usa so o ponto medio, de 2 em diante cada conjunto contem o anterior
npref=1001;              % majorante praticamente exato, usado como referencia
razao=zeros(3,2,length(nps));
fprintf('Anexo A (e = %.0e)\n',tol)
fprintf('%-10s %-17s %5s %28s %10s %8s %10s\n','funcao','regra','np','m','M/Mref','|I-Ih|','cumpre')
for k=1:3
    for r=1:2
        [~,~,~,mref]=integ_adapt(fun{k},XT{k},tol,regras{r},false,npref);
        for j=1:length(nps)
            [I,~,~,m]=integ_adapt(fun{k},XT{k},tol,regras{r},false,nps(j));
            razao(k,r,j)=sum(m)/sum(mref);
            E=abs(Iref(k)-I);
            if E<tol, ok='sim'; else, ok='nao'; end
            fprintf('%-10s %-17s %5d %28s %10.4f %8.2e %10s\n',nome{k},nomeregra{r},nps(j), ...
                sprintf('[%d %d %d]',m),razao(k,r,j),E,ok)
        end
    end
end
% eixo de np em escala logaritmica porque np quase duplica de um valor para o seguinte,
% o que deixa os pontos igualmente espacados no grafico
figA=figure('Position',[50 50 700 450],'DefaultAxesFontSize',10);
for k=1:3
    for r=1:2
        semilogx(nps,squeeze(razao(k,r,:)),['-' marc{k}],'Color',cor{r},'LineWidth',1.5); hold on
    end
end
set(gca,'XTick',nps)
xlabel('pontos do majorante, n_p'); ylabel('M / M_{ref}'); grid on
legend({'sin, esq','sin, trap','exp, esq','exp, trap','ln ln, esq','ln ln, trap'},'Location','southeast')

%% Anexo B: erro das diferencas finitas em funcao do passo
% derivadas exatas escritas a mao (so para verificar, nao entram no metodo)
f1={@(x) cos(x), @(x) -2*x.*exp(-x.^2), @(x) 1./(x.*log(x))};
f2={@(x) -sin(x), @(x) (4*x.^2-2).*exp(-x.^2), @(x) -(log(x)+1)./(x.^2.*log(x).^2)};
f3={@(x) -cos(x), @(x) (-8*x.^3+12*x).*exp(-x.^2), ...
    @(x) (2*log(x).^2+3*log(x)+2)./(x.^3.*log(x).^3)};
f4={@(x) sin(x), @(x) (16*x.^4-48*x.^2+12).*exp(-x.^2), ...
    @(x) -(6*log(x).^3+11*log(x).^2+12*log(x)+6)./(x.^4.*log(x).^4)};
x0=[1 0.5 1.5];          % ponto de teste de cada funcao
u=eps/2;
d1=(3*u)^(1/3); d2=(48*u)^(1/4);   % os passos usados no trabalho (M3 = M4 = 1)
delta=logspace(-12,-1,100);
% escala logaritmica nos dois eixos porque delta vai de 1e-12 a 1e-1 e o erro varia ao longo
% de varias ordens de grandeza. Assim os dois ramos do V aparecem como retas
figB=figure('Position',[50 50 1100 420],'DefaultAxesFontSize',10);
fprintf('\nAnexo B\n')
fprintf('%-10s %6s %8s %8s %10s %10s %10s %10s\n','funcao','x','M3','M4','f1 nosso','f1 real','f2 nosso','f2 real')
for k=1:3
    f=fun{k}; x=x0(k);
    xg=linspace(XT{k}(1),XT{k}(end),200001);
    M3=max(abs(f3{k}(xg))); M4=max(abs(f4{k}(xg)));
    r1=(3*u/M3)^(1/3); r2=(48*u/M4)^(1/4);       % passos otimos com os M reais
    e1=@(d) abs((f(x+d)-f(x-d))./(2*d)-f1{k}(x));
    e2=@(d) abs((f(x+d)-2*f(x)+f(x-d))./d.^2-f2{k}(x));
    subplot(1,2,1)
    loglog(delta,e1(delta),'-','Color',corf{k},'LineWidth',1.2); hold on
    subplot(1,2,2)
    loglog(delta,e2(delta),'-','Color',corf{k},'LineWidth',1.2); hold on
    fprintf('%-10s %6.2f %8.3g %8.3g %10.1e %10.1e %10.1e %10.1e\n',nome{k},x,M3,M4,e1(d1),e1(r1),e2(d2),e2(r2))
end
subplot(1,2,1)
yl=ylim; plot([d1 d1],yl,'k--')
xlabel('\delta'); ylabel('erro de (12)'); title('(a) primeira derivada'); grid on
legend([nome {'\delta usado'}],'Location','north')
subplot(1,2,2)
yl=ylim; plot([d2 d2],yl,'k--')
xlabel('\delta'); ylabel('erro de (13)'); title('(b) segunda derivada'); grid on
legend([nome {'\delta usado'}],'Location','north')

%% Anexo C: tempo e memoria para tolerancias menores
fprintf('\nAnexo C\n')
fprintf('%-10s %-17s %8s %12s %10s %12s\n','funcao','regra','e','soma m','tempo (s)','x (MB)')
for k=1:3
    for r=1:2
        if r==1, ts=[1e-5 1e-6 1e-7]; else, ts=[1e-5 1e-6 1e-7 1e-8]; end
        for t=ts
            tic
            [~,~,x,m]=integ_adapt(fun{k},XT{k},t,regras{r});
            tempo=toc;
            info=whos('x');                          % memoria ocupada pelos nos
            fprintf('%-10s %-17s %8.0e %12d %10.2f %12.3g\n',nome{k},nomeregra{r},t,sum(m),tempo,info.bytes/1e6)
            clear x
        end
    end
end

% fundo branco para o relatorio
nomesfig={'figA_pontos_majorante.png','figB_passo_df.png'};
figs=[figA figB];
for i=1:2
    try
        theme(figs(i),'light')
    catch
    end
    set(figs(i),'Color','w')
    exportgraphics(figs(i),nomesfig{i},'Resolution',300,'BackgroundColor','white')
end
