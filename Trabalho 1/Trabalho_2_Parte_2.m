%##########################################################################
%#               UNIVERSIDADE FEDERAL DE JUIZ DE FORA                     #
%#              GUSTAVO LEAL SILVA E SOUZA - 201469055B                   #
%##########################################################################
clear all; close all; clc;

nppc = 256;                         % Número de pontos por ciclo
f1 = 60;                            % Frequência fundamental
Nc = 1;                             % Número de ciclos
Ts = 1/(f1*nppc);                   % Periodo de amostragem
Fs = 1/Ts;                          % Frequência de amostragem
tempo = (0:Nc*nppc-1)*Ts;           % Linha de tempo

delt_f = Fs/(nppc*Nc);
ini = 256+1;
fim = ini+nppc;

sim barras15;                       % Start para simulação do SIMULINK

t_5 = tensao1.signals.values;       % Medidor de Tensão linha 5
t_10 = tensao1.signals.values;      % Medidor de Tensão linha 10

c_5 = corrente1.signals.values;     % Medidor de Corrente linha 5
c_10 = corrente1.signals.values;    % Medidor de Corrente linha 10

T_5 = 2*fft(t_5(ini:fim-1))/length(t_5(ini:fim-1));
T_10 = 2*fft(t_10(ini:fim-1))/length(t_10(ini:fim-1));

C_5 = 2*fft(c_5(ini:fim-1))/length(c_5(ini:fim-1));
C_10 = 2*fft(c_10(ini:fim-1))/length(c_10(ini:fim-1));

% Plot
subplot(2,2,1);
plot(tempo,(fft(T_5)),'b')
title('Transformada de Fourier V(t) Linhas 5');
ylabel('Tensão [ V ]');
xlabel('Tempo [ s ]');
grid on;

subplot(2,2,2);
plot(tempo,(fft(C_5)),'b');
title('Transformada de Fourier I(t) Linhas 5');
ylabel('Corrente [ A ]');
xlabel('Tempo [ s ]');
grid on;

subplot(2,2,3);
plot(tempo,(fft(T_10)),'b');
title('Transformada de Fourier V(t) Linhas 10');
ylabel('Tensão [ V ]');
xlabel('Tempo [ s ]');
grid on;

subplot(2,2,4);
plot(tempo,(fft(C_10)),'b');
title('Transformada de Fourier I(t) Linhas 10');
ylabel('Corrente [ A ]');
xlabel('Tempo [ s ]');
grid on;