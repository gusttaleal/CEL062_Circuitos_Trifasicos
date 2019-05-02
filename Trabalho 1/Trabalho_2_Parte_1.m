%##########################################################################
%#               UNIVERSIDADE FEDERAL DE JUIZ DE FORA                     #
%#              GUSTAVO LEAL SILVA E SOUZA - 201469055B                   #
%##########################################################################
%% Função Quadrada (Função Par > bn = 0)
clear all; close all; clc;

syms t              % Variavel de tempo

% Correçao na Escala de Tempo | Inicio
    nppc = 256;                 % Número de pontos por ciclo
    f1 = 10;                    % Frequência fundamental
    Nc = 3;                     % Número de ciclos
    Ts = 1/(f1*nppc);           % Periodo de amostragem
    tempo = (0:Nc*nppc-1)*Ts;   % Linha de tempo 
% Correçao na Escala de Tempo | Fim

% Parâmetros da função
V = 2;              % Função de t, Tensão
I = 0.1;            % Função de t, Corrente

T = 0.1;            % Período da função
w = 2*pi/T;         % Frequência Angular

N = 15;              % Número referente aos coeficientes 'n' da série
n = 1:N;

a = 0.025;          % Limite de Integração
b = 0.075;          % Limite de Integração
c = 0.1;            % Limite de Integração

%   TENSÃO   
% Coeficientes a0, an e bn da série de Fourier
    V_a0 = (1/T)*(int(V,t,0,a)+int(V,t,b,c));
    V_an = (2/T)*(int(V*cos(n*w*t),t,0,a)+int(V*cos(n*w*t),t,b,c));
    V_bn = (2/T)*(int(V*sin(n*w*t),t,0,a)+int(V*sin(n*w*t),t,b,c));
      
% Série de Fourier
    SF_V = V_a0 + sum(V_an.*cos(n*w*t)) + sum(V_bn.*sin(n*w*t));

% Forma Alternativa
    V_An = sqrt((V_an.^2)+(V_bn.^2));
    V_teta = Teta(V_an,V_bn,N);

% Valor Eficaz de onda
    V_rms = sqrt(V_a0^2 + (sum((V_An.^2))/2));
    
% THD para Tensão
    THDv = (sqrt( sum((V_An.^2))-(V_An(1)^2)) / V_An(1) )*100;
    
%   CORRENTE
% Coeficientes a0, an e bn da série de Fourier
    I_a0 = (1/T)*(int(I,t,0,a)+int(I,t,b,c));
    I_an = (2/T)*(int(I*cos(n*w*t),t,0,a)+int(I*cos(n*w*t),t,b,c));
    I_bn = (2/T)*(int(I*sin(n*w*t),t,0,a)+int(I*sin(n*w*t),t,b,c));

% Série de Fourier
    SF_I = I_a0 + sum(I_an.*cos(n*w*t)) + sum(I_bn.*sin(n*w*t));
    
% Forma Alternativa
    I_An = sqrt((I_an.^2)+(I_bn.^2));
    I_teta = Teta(I_an,I_bn,N);

% Valor Eficaz de onda
    I_rms = sqrt(I_a0^2 + (sum((I_An.^2))/2));
    
% THD para Tensão
    THDi = (sqrt( sum((I_An.^2))-(I_An(1)^2)) / I_An(1) )*100;

% POTÊNCIA
    P = V_a0*I_a0 + (sum(V_An.*I_An.*cos(deg2rad(V_teta-I_teta)))/2);
% Fator de Potência
    FP = P/(V_rms*I_rms);
% Representação
    disp_plot(V_a0, V_an, V_bn, SF_V,V_An, V_teta, V_rms, THDv,...
I_a0, I_an, I_bn, SF_I, I_An, I_teta, I_rms, THDi, P, FP, N, n, t, tempo);
%% Gráfico 1
clear all; close all; clc;

syms x

ez1 = ezplot(2*rectangularPulse(x*20), [-.05 .05]); hold on; grid on;
ez2 = ezplot(2*rectangularPulse((x+.1)*20), [-.15 -.05]);
ez3 = ezplot(2*rectangularPulse((x-.1)*20), [.05 .15]);

stem(0,3,'k')
stem(1,0,'k')

set(ez1,'color','b');
set(ez2,'color','b');
set(ez3,'color','b');

title('Função V(t) - Onda Quadrada')
ylabel('Tensão [ V ]')
xlabel('Tempo [ s ]')
axis([-.15 .15 -.2 2.2])
%% Função Triangulo(Função Par > bn = 0)

clear all; close all; clc;

syms t          % Variavel de tempo

% Correçao na Escala de Tempo | Inicio
    nppc = 256;                 % Número de pontos por ciclo
    f1 = 10;                    % Frequência fundamental
    Nc = 3;                     % Número de ciclos
    Ts = 1/(f1*nppc);           % Periodo de amostragem
    tempo = (0:Nc*nppc-1)*Ts;   % Linha de tempo 
% Correçao na Escala de Tempo | Fim

% Parâmetros da função
    V1 = 40*t;      % Função de t, Tensão
    V2 = -40*t+4;

    I1 = 2*t;       % Função de t, Corrente
    I2 = -2*t+0.2;

    T = 0.1;        % Período da função
    w = 2*pi/T;     % Frequência Angular

    N = 15;          % Número referente aos coeficientes 'n' da série
    n = 1:N;     
    
% Limites de Integração
    a = 0.05;      
    b = 0.1;
        
%   TENSÃO   
% Coeficientes a0, an e bn da série de Fourier
    V_a0 = (1/T)*(int(V1,t,0,a)+int(V2,t,a,b));
    V_an = (2/T)*(int(V1*cos(n*w*t),t,0,a)+int(V2*cos(n*w*t),t,a,b));
    V_bn = (2/T)*(int(V1*sin(n*w*t),t,0,a)+int(V2*sin(n*w*t),t,a,b));
    
% Série de Fourier
    SF_V = V_a0 + sum(V_an.*cos(n*w*t)) + sum(V_bn.*sin(n*w*t));

% Forma Alternativa
    V_An = sqrt((V_an.^2)+(V_bn.^2));
    V_teta = Teta(V_an,V_bn,N);

% Valor Eficaz de onda
    V_rms = sqrt(V_a0^2 + (sum((V_An.^2))/2));
    
% THD para Tensão
    THDv = (sqrt( sum((V_An.^2))-(V_An(1)^2)) / V_An(1) )*100;
    
%   CORRENTE
% Coeficientes a0, an e bn da série de Fourier
    I_a0 = (1/T)*(int(I1,t,0,a)+int(I2,t,a,b));
    I_an = (2/T)*(int(I1*cos(n*w*t),t,0,a)+int(I2*cos(n*w*t),t,a,b));
    I_bn = (2/T)*(int(I1*sin(n*w*t),t,0,a)+int(I2*sin(n*w*t),t,a,b));

% Série de Fourier
    SF_I = I_a0 + sum(I_an.*cos(n*w*t)) + sum(I_bn.*sin(n*w*t));
    
% Forma Alternativa
    I_An = sqrt((I_an.^2)+(I_bn.^2));
    I_teta = Teta(I_an,I_bn,N);

% Valor Eficaz de onda
    I_rms = sqrt(I_a0^2 + (sum((I_An.^2))/2));
    
% THD para Tensão
    THDi = (sqrt( sum((I_An.^2))-(I_An(1)^2)) / I_An(1) )*100;

% POTÊNCIA
    P = V_a0*I_a0 + (sum(V_An.*I_An.*cos(deg2rad(V_teta-I_teta)))/2);
% Fator de Potência
    FP = P/(V_rms*I_rms); 
    
% Representação
    disp_plot(V_a0, V_an, V_bn, SF_V,V_An, V_teta, V_rms, THDv,...
I_a0, I_an, I_bn, SF_I, I_An, I_teta, I_rms, THDi, P, FP, N, n, t, tempo);
%% Gráfico 2
clear all; close all; clc;

a = [2 0 2 0 2 0 2];
b = [-.15 -.1 -.05 0 .05 .1 .15];

plot(b,a, 'b'); hold on; grid on;

stem(0,3,'k')
stem(1,0,'k')

title('Função V(t) - Onda Triangular')
ylabel('Tensão [ V ]')
xlabel('Tempo [ s ]')
axis([-.15 .15 -.2 2.2])
