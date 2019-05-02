%##########################################################################
%#               UNIVERSIDADE FEDERAL DE JUIZ DE FORA                     #
%#              GUSTAVO LEAL SILVA E SOUZA - 201469055B                   #
%##########################################################################
function disp_plot( V_a0, V_an, V_bn, SF_V,V_An, V_teta, V_rms, THDv,...
  I_a0, I_an, I_bn, SF_I, I_An, I_teta, I_rms, THDi, P, FP, N, n, t, tempo)
    %   DISPLAY
    disp('--- Tens�o ---'); 
    disp(' ');
    fprintf('V_a0 = %.3f\n', double(V_a0));
    disp(' ');
    fprintf('V_a%d = %.3f\n', cat(1,n,double(V_an)));
    disp(' ');
    fprintf('V_b%d = %.3f\n', cat(1,n,double(V_bn)));
    disp(' ');
    disp(SF_V);
    fprintf('V_A%d = %.3f\n',cat(1,n,double(V_An)));
    disp(' ');
    fprintf('V_T%d = %.3f�\n',cat(1,n,double(V_teta)));
    disp(' ');
    fprintf('V_rms = %.3f\n',double(V_rms));
    disp(' ');
    fprintf('V_THD = %.3f\n',double(THDv));
    
    disp('--- Corrente ---'); 
    disp(' ');
    fprintf('I_a0 = %.3f\n', double(I_a0));
    disp(' ');
    fprintf('I_a%d = %.3f\n', cat(1,n,double(I_an)));
    disp(' ');
    fprintf('I_b%d = %.3f\n', cat(1,n,double(I_bn)));
    disp(' ');
    disp(SF_I);
    fprintf('I_A%d = %.3f\n',cat(1,n,double(I_An)));
    disp(' ');
    fprintf('I_T%d = %.3f�\n',cat(1,n,double(I_teta)));
    disp(' ');
    fprintf('I_rms = %.3f\n',double(I_rms));
    disp(' ');
    fprintf('I_THD = %.3f\n',double(THDi));   
    
    disp('--- Pot�ncia ---');
    disp(' ');
    fprintf('P = %.3f\n',double(P));
    fprintf('FP = %.3f\n',double(FP));
    
%   PLOT    
    subplot(2,2,1)   
    stem(n,V_An,'*b');
    title('Fun��o V(t) - Onda Quadrada | M�dulo')
    ylabel('Tens�o [ V ]')
    xlabel('Tempo [ s ]')
    grid on;
    
    subplot(2,2,2)   
    stem(n,V_teta,'*b');
    title('Fun��o V(t) - Onda Quadrada | Fase')
    ylabel('Tens�o [ teta � ]')
    xlabel('Tempo [ s ]')
    axis([0 N+1 -90 90]);
    grid on;
    
    subplot(2,2,3)    
    stem(n,I_An,'*b');
    title('Fun��o I(t) - Onda Quadrada | M�dulo')
    ylabel('Corrente [ A ]')
    xlabel('Tempo [ s ]')
    grid on;
    
    subplot(2,2,4)    
    stem(n,I_teta,'*b');
    axis([0 N+1 -90 90]);
    title('Fun��o I(t) - Onda Quadrada | Fase')
    ylabel('Corrente [ teta � ]')
    xlabel('Tempo [ s ]')
    grid on;
    
    figure
    plot(tempo,double(subs(SF_V,t,tempo)),'b')
    grid on;
    title('Fun��o V(t) - Onda Quadrada | S�rie de Fourier')
    ylabel('Tens�o [ V ]')
    xlabel('Tempo [ s ]')

end

