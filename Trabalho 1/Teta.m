%##########################################################################
%#               UNIVERSIDADE FEDERAL DE JUIZ DE FORA                     #
%#              GUSTAVO LEAL SILVA E SOUZA - 201469055B                   #
%##########################################################################
function Teta = Teta(an,bn,N)

Teta = ones(1,N);

for n = 1:N
    % Se +an + jbn ou +an - jbn
    if(an(n)>0 && bn(n)>0 || an(n)>0 && bn(n)<0)
        Teta(n) = rad2deg(-atan(bn(n)/an(n)));
    end
    % Se -an + jbn
    if(an(n)<0 && bn(n)>0)
        Teta(n) = rad2deg(-atan(bn(n)/an(n)))+180;
    end
    % Se -an - jbn
    if(an(n)<0 && bn(n)<0)
        Teta(n) = rad2deg(-atan(bn(n)/an(n)))-180;
    end
   
    if (an(n)~=0 && bn(n)==0)
        Teta(n) = 0;
    end
    
    if (an(n)==0 && bn(n)~=0)
        Teta(n) = 90;
    end
    
    if (an(n)==0 && bn(n)==0)
        Teta(n) = 0;
    end
end

end

