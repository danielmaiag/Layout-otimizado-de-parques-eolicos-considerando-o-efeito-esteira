function [FOB,ef,P_gerada] = F_objetivo(M,n,fk,linhas,colunas,alt_rotor,diam_rotor,coef_empuxo,f_axial,rugosidade,vel_vento_livre,arraste,raio_esteira)
    %% Cálculo de velocidades e potências nos aerogeradores de acordo com a direção do vento

    %% Vento a 0º
    if fk(1) ~= 0

        % 1) ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO

        % Definição da quantidade de turbinas por coluna
        b = sum(M);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        pos_turb0=zeros(size(M));
        for k=1:colunas;
          for w=1:linhas;
            if M(w,k)==1
                pos_turb0(w,k)=w;
            else
                pos_turb0(w,k)=0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = linhas;
        a=1;
        while a~=0;
            a=0;
            for k=1:colunas;
              for w=1:linhas-1;
                if pos_turb0(w,k) == 0 && pos_turb0(w+1,k) ~= 0
                        pos_turb0(w,k) = pos_turb0(w+1,k);
                        pos_turb0(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end

        % 2) CÁLCULO DA VELOCIDADE MÉDIA DO VENTO EM CADA TURBINA

        velocidade0=zeros(linhas,colunas);
        for k=1:colunas % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b(k)==1;
                velocidade0(1,k)=vel_vento_livre;
            else if b(k)>1;
                    velocidade0(1,k)=vel_vento_livre;
                    for w=(2:b(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma0=0;
                        vel_parciais0=zeros(b(k)-1,b(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais0(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb0(w,k)-pos_turb0(u,k)))*5*diam_rotor))/(raio_esteira)))^2));
                            soma0=soma0+((1-(vel_parciais0(u,w-1)/12))^2);
                        end
                        velocidade0(w,k)=12*(1-sqrt(soma0));
                    end
                end
            end
        end

        % 3) CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P0 = zeros(linhas,colunas);
        for k=1:colunas
            for w=1:linhas
                if velocidade0(w,k) <= 2.3
                    P0(w,k)=0;
                else if velocidade0(w,k) >= 2.3 && velocidade0(w,k) <= 12.8
                        P0(w,k)=0.3*velocidade0(w,k)^3;
                    else if velocidade0(w,k) >= 12.8 && velocidade0(w,k) <= 18
                            P0(w,k)=630;
                        else if velocidade0(w,k)>=18
                                P0(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total0=sum(sum(P0));
    else P_total0=0;
    end

    %% VENTO A 45º

    if fk(2) ~= 0

        % A matriz M será rotacionada de forma que a lógica utilizada para o
        % cálculo das velocidades possa ser aproveitada.
        % A matriz M45 é, portanto, a matriz M reescrita de forma que os
        % geradores dentro de uma esteira sejam representados na mesma coluna.

        M45=zeros(linhas,(2*linhas-1));
        for k=(1:(2*linhas-1))
              if k <= linhas
                for w=1:k
                    M45(w,k) = M(w,(k+1-w));
                end
              end
              if k > linhas
                d=k-(linhas-1);
                for w=1:(2*linhas-k) 
                    M45(w,k) = M(d,(k+1-d));
                    d=d+1;
                end
              end
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------

        % Definição da quantidade de turbinas por coluna
        b45 = sum(M45);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        [L,C] = size(M45);
        pos_turb45=zeros(size(M45));
        for k=1:C;
          for w=1:L;
            if M45(w,k)==1
                pos_turb45(w,k)=w;
            else
                pos_turb45(w,k)=0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = L;
        a=1;
        while a~=0;
            a=0;
            for k=1:C;
              for w=1:L-1;
                if pos_turb45(w,k) == 0 && pos_turb45(w+1,k) ~= 0
                        pos_turb45(w,k) = pos_turb45(w+1,k);
                        pos_turb45(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end
        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade45=zeros(L,C);
        for k=1:C % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b45(k)==1;
                velocidade45(1,k)=vel_vento_livre;
            else if b45(k)>1;
                    velocidade45(1,k)=vel_vento_livre;
                    for w=(2:b45(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma45=0;
                        vel_parciais45=zeros(b45(k)-1,b45(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais45(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb45(w,k)-pos_turb45(u,k)))*5*diam_rotor*sqrt(2)))/(raio_esteira)))^2));
                            soma45=soma45+((1-(vel_parciais45(u,w-1)/12))^2);
                        end
                        velocidade45(w,k)=12*(1-sqrt(soma45));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA
        % ------------------------------------------------------------------------------------------
        % Cálculo da potência gerada pela configuração de geradores considerada.
        P45 = zeros(L,C);
        for k=1:C
            for w=1:L
                if velocidade45(w,k) <= 2.3
                    P45(w,k)=0;
                else if velocidade45(w,k) >= 2.3 && velocidade45(w,k) <= 12.8
                        P45(w,k)=0.3*velocidade45(w,k)^3;
                    else if velocidade45(w,k) >= 12.8 && velocidade45(w,k) <= 18
                            P45(w,k)=630;
                        else if velocidade45(w,k)>=18
                                P45(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total45=sum(sum(P45));
    else P_total45=0;
    end

    %% VENTO A 90º

    if fk(3) ~= 0

            % A matriz M será rotacionada de forma que a lógica utilizada para o
            % cálculo das velocidades possa ser aproveitada.
            % A matriz M90 é, portanto, a matriz M reescrita de forma que os
            % geradores dentro de uma esteira sejam representados na mesma coluna.

            M90=zeros(linhas,colunas);
        for k = 1:colunas
            d=linhas;
            for w = 1:linhas
                M90(w,k)=M(k,d);
                d=d-1;
            end
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------
        % Definição da quantidade de turbinas por coluna
        b90 = sum(M90);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        pos_turb90=zeros(size(M90));
        for k=1:colunas;
          for w=1:linhas;
            if M90(w,k)==1
                pos_turb90(w,k)=w;
            else
                pos_turb90(w,k)=0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = linhas;
        a=1;
        while a~=0;
            a=0;
            for k=1:colunas;
              for w=1:linhas-1;
                if pos_turb90(w,k) == 0 && pos_turb90(w+1,k) ~= 0
                        pos_turb90(w,k) = pos_turb90(w+1,k);
                        pos_turb90(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end

        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade90=zeros(linhas,colunas);
        for k=1:colunas % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b90(k)==1;
                velocidade90(1,k)=vel_vento_livre;
            else if b90(k)>1;
                    velocidade90(1,k)=vel_vento_livre;
                    for w=(2:b90(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma90=0;
                        vel_parciais90=zeros(b90(k)-1,b90(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais90(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb90(w,k)-pos_turb90(u,k)))*5*diam_rotor))/(raio_esteira)))^2));
                            soma90=soma90+((1-(vel_parciais90(u,w-1)/12))^2);
                        end
                        velocidade90(w,k)=12*(1-sqrt(soma90));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P90 = zeros(linhas,colunas);
        for k=1:colunas
            for w=1:linhas
                if velocidade90(w,k) <= 2.3
                    P90(w,k)=0;
                else if velocidade90(w,k) >= 2.3 && velocidade90(w,k) <= 12.8
                        P90(w,k)=0.3*velocidade90(w,k)^3;
                    else if velocidade90(w,k) >= 12.8 && velocidade90(w,k) <= 18
                            P90(w,k)=630;
                        else if velocidade90(w,k)>=18
                                P90(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total90=sum(sum(P90));
    else P_total90=0;
    end

    %% VENTO A 135º

    if fk(4) ~= 0

        % A matriz M será rotacionada de forma que a lógica utilizada para o
        % cálculo das velocidades possa ser aproveitada.
        % A matriz M135 é, portanto, a matriz M reescrita de forma que os
        % geradores dentro de uma esteira sejam representados na mesma coluna.

        M135=zeros(linhas,(2*linhas-1));
        for k=(1:(2*linhas-1))
              if k < linhas
                d=linhas;
                p=k;
                for w=1:k
                    M135(w,k) = M(d,p);
                    d=d-1;
                    p=p-1;
                end
              else if k == linhas
                      d=k;
                      for w=1:k
                              M135(w,k)=M(d,d);
                              d=d-1;
                      end
              else if k > linhas
                d=linhas;
                p=(2*linhas-k);
                for w=1:(2*linhas-k) 
                    M135(w,k) = M(p,d);
                    d=d-1;
                    p=p-1;
                end
                  end
                  end
              end
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------
        % Definição da quantidade de turbinas por coluna
        b135 = sum(M135);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        [L,C] = size(M135);
        pos_turb135=zeros(size(M135));
        for k=1:C;
          for w=1:L;
            if M135(w,k)==1
                pos_turb135(w,k)=w;
            else
                pos_turb135(w,k)=0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = L;
        a=1;
        while a~=0;
            a=0;
            for k=1:C;
              for w=1:L-1;
                if pos_turb135(w,k) == 0 && pos_turb135(w+1,k) ~= 0
                        pos_turb135(w,k) = pos_turb135(w+1,k);
                        pos_turb135(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end
        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade135=zeros(L,C);
        for k=1:C % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b135(k)==1;
                velocidade135(1,k)=vel_vento_livre;
            else if b135(k)>1;
                    velocidade135(1,k)=vel_vento_livre;
                    for w=(2:b135(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma135=0;
                        vel_parciais135=zeros(b135(k)-1,b135(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais135(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb135(w,k)-pos_turb135(u,k)))*5*diam_rotor*sqrt(2)))/(raio_esteira)))^2));
                            soma135=soma135+((1-(vel_parciais135(u,w-1)/12))^2);
                        end
                        velocidade135(w,k)=12*(1-sqrt(soma135));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P135 = zeros(L,C);
        for k=1:C
            for w=1:L
                if velocidade135(w,k) <= 2.3
                    P135(w,k)=0;
                else if velocidade135(w,k) >= 2.3 && velocidade135(w,k) <= 12.8
                        P135(w,k)=0.3*velocidade135(w,k)^3;
                    else if velocidade135(w,k) >= 12.8 && velocidade135(w,k) <= 18
                            P135(w,k)=630;
                        else if velocidade135(w,k)>=18
                                P135(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total135=sum(sum(P135));
    else P_total135=0;
    end

    %% VENTO A 180º 

    if fk(5) ~= 0

        % A matriz M será rotacionada de forma que a lógica utilizada para o
        % cálculo das velocidades possa ser aproveitada.
        % A matriz M180 é, portanto, a matriz M reescrita de forma que os
        % geradores dentro de uma esteira sejam representados na mesma coluna.

        M180=zeros(linhas,colunas);
        for k = 1:colunas
            d=linhas;
            for w = 1:linhas
                M180(w,k)=M(d,k);
                d=d-1;
            end
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------
        % Definição da quantidade de turbinas por coluna
        b180 = sum(M180);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        pos_turb180=zeros(size(M180));
        for k=1:colunas;
          for w=1:linhas;
            if M180(w,k)==1
                pos_turb180(w,k)=w;
            else
                pos_turb180(w,k)=0;
            end
          end
        end
        % Reorganização da matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = linhas;
        a=1;
        while a~=0;
            a=0;
            for k=1:colunas;
              for w=1:linhas-1;
                if pos_turb180(w,k) == 0 && pos_turb180(w+1,k) ~= 0
                        pos_turb180(w,k) = pos_turb180(w+1,k);
                        pos_turb180(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end

        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade180=zeros(linhas,colunas);
        for k=1:colunas % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b180(k)==1;
                velocidade180(1,k)=vel_vento_livre;
            else if b180(k)>1;
                    velocidade180(1,k)=vel_vento_livre;
                    for w=(2:b180(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma180=0;
                        vel_parciais180=zeros(b180(k)-1,b180(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais180(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb180(w,k)-pos_turb180(u,k)))*5*diam_rotor))/(raio_esteira)))^2));
                            soma180=soma180+((1-(vel_parciais180(u,w-1)/12))^2);
                        end
                        velocidade180(w,k)=12*(1-sqrt(soma180));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P180 = zeros(linhas,colunas);
        for k=1:colunas
            for w=1:linhas
                if velocidade180(w,k) <= 2.3
                    P180(w,k)=0;
                else if velocidade180(w,k) >= 2.3 && velocidade180(w,k) <= 12.8
                        P180(w,k)=0.3*velocidade180(w,k)^3;
                    else if velocidade180(w,k) >= 12.8 && velocidade180(w,k) <= 18
                            P180(w,k)=630;
                        else if velocidade180(w,k)>=18
                                P180(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total180=sum(sum(P180));
    else P_total180=0;
    end

    %% VENTO A 225º 

    if fk(6) ~= 0 

        % A matriz M será rotacionada de forma que a lógica utilizada para o
        % cálculo das velocidades possa ser aproveitada.
        % A matriz M225 é, portanto, a matriz M reescrita de forma que os
        % geradores dentro de uma esteira sejam representados na mesma coluna.

        M225=zeros(linhas,(2*linhas-1));
        for k=(1:(2*linhas-1))
              if k <= linhas
                for w=1:k
                    M225(w,k) = M((k+1-w),w);
                end
              end
              if k > linhas
                d=k-(linhas-1);
                for w=1:(2*linhas-k) 
                    M225(w,k) = M((k+1-d),d);
                    d=d+1;
                end
              end
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------
        % Definição da quantidade de turbinas por coluna
        b225 = sum(M225);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        [L,C] = size(M225);
        pos_turb225=zeros(size(M225));
        for k=1:C;
          for w=1:L;
            if M225(w,k)==1
                pos_turb225(w,k)=w;
            else
                pos_turb225(w,k)=0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = L;
        a=1;
        while a~=0;
            a=0;
            for k=1:C;
              for w=1:L-1;
                if pos_turb225(w,k) == 0 && pos_turb225(w+1,k) ~= 0
                        pos_turb225(w,k) = pos_turb225(w+1,k);
                        pos_turb225(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end

        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade225=zeros(L,C);
        for k=1:C % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b225(k)==1;
                velocidade225(1,k)=vel_vento_livre;
            else if b225(k)>1;
                    velocidade225(1,k)=vel_vento_livre;
                    for w=(2:b225(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma225=0;
                        vel_parciais225=zeros(b225(k)-1,b225(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais225(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb225(w,k)-pos_turb225(u,k)))*5*diam_rotor*sqrt(2)))/(raio_esteira)))^2));
                            soma225=soma225+((1-(vel_parciais225(u,w-1)/12))^2);
                        end
                        velocidade225(w,k)=12*(1-sqrt(soma225));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P225 = zeros(L,C);
        for k=1:C
            for w=1:L
                if velocidade225(w,k) <= 2.3
                    P225(w,k)=0;
                else if velocidade225(w,k) >= 2.3 && velocidade225(w,k) <= 12.8
                        P225(w,k)=0.3*velocidade225(w,k)^3;
                    else if velocidade225(w,k) >= 12.8 && velocidade225(w,k) <= 18
                            P225(w,k)=630;
                        else if velocidade225(w,k)>=18
                                P225(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total225=sum(sum(P225));
    else P_total225=0;
    end

    %% VENTO A 270º

    if fk(7) ~= 0;

        % A matriz M será rotacionada de forma que a lógica utilizada para o
        % cálculo das velocidades possa ser aproveitada.
        % A matriz M270 é, portanto, a matriz M reescrita de forma que os
        % geradores dentro de uma esteira sejam representados na mesma coluna.

        M270=zeros(linhas,colunas);
        d=linhas;
        for k = 1:colunas
            for w = 1:linhas
                M270(w,k)=M(d,w);
            end
            d=d-1;
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------
        % Definição da quantidade de turbinas por coluna
        b270 = sum(M270);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        pos_turb270=zeros(size(M270));
        for k=1:colunas;
          for w=1:linhas;
            if M270(w,k)==1
                pos_turb270(w,k)=w;
            else
                pos_turb270(w,k)=0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = linhas;
        a=1;
        while a~=0;
            a=0;
            for k=1:colunas;
              for w=1:linhas-1;
                if pos_turb270(w,k) == 0 && pos_turb270(w+1,k) ~= 0
                        pos_turb270(w,k) = pos_turb270(w+1,k);
                        pos_turb270(w+1,k)=0;
                        a=1;  
                end
              end
            end
        end

        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade270=zeros(linhas,colunas);
        for k=1:colunas % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b270(k)==1;
                velocidade270(1,k)=vel_vento_livre;
            else if b270(k)>1;
                    velocidade270(1,k)=vel_vento_livre;
                    for w=(2:b270(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma270=0;
                        vel_parciais270=zeros(b270(k)-1,b270(k)-1);
                        for u=1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais270(u,w-1)=vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb270(w,k)-pos_turb270(u,k)))*5*diam_rotor))/(raio_esteira)))^2));
                            soma270=soma270+((1-(vel_parciais270(u,w-1)/12))^2);
                        end
                        velocidade270(w,k)=12*(1-sqrt(soma270));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P270 = zeros(linhas,colunas);
        for k=1:colunas
            for w=1:linhas
                if velocidade270(w,k) <= 2.3
                    P270(w,k)=0;
                else if velocidade270(w,k) >= 2.3 && velocidade270(w,k) <= 12.8
                        P270(w,k)=0.3*velocidade270(w,k)^3;
                    else if velocidade270(w,k) >= 12.8 && velocidade270(w,k) <= 18
                            P270(w,k)=630;
                        else if velocidade270(w,k)>=18
                                P270(w,k)=0;
                            end
                        end
                    end
                end
            end
        end
        P_total270=sum(sum(P270));
        else P_total270=0;
    end

    %% VENTO A 315º

    if fk(8) ~= 0   

        % A matriz M será rotacionada de forma que a lógica utilizada para o
        % cálculo das velocidades possa ser aproveitada.
        % A matriz M315 é, portanto, a matriz M reescrita de forma que os
        % geradores dentro de uma esteira sejam representados na mesma coluna.

        d = linhas;
        M315 = zeros(linhas,(2*linhas-1));
        for k = (1:(2*linhas-1))
              memoria = d;
              if k < linhas
                for w = 1:k
                    M315(w,k) = M(d,w);
                    d = d+1;
                end
                d = memoria-1;
              else if k == linhas
                      for w = 1:k
                              M315(w,k) = M(w,w);
                      end
              else if k > linhas
                p = (k-(linhas-1));
                for w = 1:(2*linhas-k) 
                    M315(w,k) = M(w,p);
                    p = p+1;
                end
                  end
                  end
              end
        end

        % -----------------------------------------------------------------------------------------------
        % ANÁLISE DA DISTÂNCIA ENTRE AS TURBINAS NO TERRENO
        % -----------------------------------------------------------------------------------------------
        % Definição da quantidade de turbinas por coluna
        b315 = sum(M315);
        % Criação de matriz com identificadores das linhas onde estão posicionadas
        % as turbinas
        [L,C] = size(M315);
        pos_turb315 = zeros(size(M315));
        for k = 1:C;
          for w = 1:L;
            if M315(w,k) == 1
                pos_turb315(w,k) = w;
            else
                pos_turb315(w,k) = 0;
            end
          end
        end
        % Reorganização a matriz de posições, passando os algarismos diferentes de 0
        % para cima. A variável 'a' é o critério de parada.
        w = L;
        a = 1;
        while a ~= 0;
            a = 0;
            for k = 1:C;
              for w = 1:L-1;
                if pos_turb315(w,k) == 0 && pos_turb315(w+1,k) ~= 0
                        pos_turb315(w,k) = pos_turb315(w+1,k);
                        pos_turb315(w+1,k) = 0;
                        a = 1;  
                end
              end
            end
        end

        % ------------------------------------------------------------------------------------------------
        % MODELAGEM MATEMÁTICA DO PROBLEMA
        % ---------------------------------------------------------------------------------------------
        % Cálculo das velocidades médias do vento em cada turbina

        velocidade315 = zeros(L,C);
        for k = 1:C % Este laço passa por todas as colunas da matriz que representa o terreno.
            if b315(k) == 1;
                velocidade315(1,k) = vel_vento_livre;
            else if b315(k) > 1;
                    velocidade315(1,k) = vel_vento_livre;
                    for w = (2:b315(k)) %Este laço passa por cada uma das turbinas em uma mesma coluna.
                        soma315 = 0;
                        vel_parciais315 = zeros(b315(k)-1,b315(k)-1);
                        for u = 1:(w-1)
                            % Criação de vetor auxiliar que irá computar a velocidade
                            % parcial em cada turbina devido às demais turbinas
                            vel_parciais315(u,w-1) = vel_vento_livre*(1-((2*f_axial)/(1+((arraste*(abs((pos_turb315(w,k)-pos_turb315(u,k)))*5*diam_rotor*sqrt(2)))/(raio_esteira)))^2));
                            soma315 = soma315+((1-(vel_parciais315(u,w-1)/12))^2);
                        end
                        velocidade315(w,k) = 12*(1-sqrt(soma315));
                    end
                end
            end
        end
        % --------------------------------------------------------------------------------------------
        % CÁLCULOS RELACIONADOS À POTÊNCIA

        % Cálculo da potência gerada pela configuração de geradores considerada.
        P315 = zeros(L,C);
        for k = 1:C
            for w = 1:L
                if velocidade315(w,k) <= 2.3
                    P315(w,k) = 0;
                else if velocidade315(w,k) >= 2.3 && velocidade315(w,k) <= 12.8
                        P315(w,k) = 0.3*velocidade315(w,k)^3;
                    else if velocidade315(w,k) >= 12.8 && velocidade315(w,k) <= 18
                            P315(w,k) = 630;
                        else if velocidade315(w,k) >= 18
                                P315(w,k) = 0;
                            end
                        end
                    end
                end
            end
        end
        P_total315 = sum(sum(P315));
    else P_total315 = 0;
    end
    %% Cálculo da potência total e da eficiência da usina

    % Criação do vetor P, que contém as potências máximas para cada direção do
    % vento considerada.
    P = [P_total0,P_total45,P_total90,P_total135,P_total180,P_total225,P_total270,P_total315];

    % Cálculo da potência total, ponderada pela probabilidade de incidência dos
    % ventos em cada direção.
    P_gerada = sum(P.*fk);

    % Cálculo da potência gerada por um único gerador livre do efeito esteira
    if vel_vento_livre <= 2.3
        PL = 0;
    else if vel_vento_livre >= 2.3 && vel_vento_livre <= 12.8
            PL = 0.3*vel_vento_livre^3;
        else if vel_vento_livre >= 12.8 && vel_vento_livre <= 18
                PL = 630;
            else if vel_vento_livre >= 18
                    PL = 0;
                end
            end
        end
    end

    % Potência total na usina considerada se não houvesse o efeito esteira
    P_livre = n*PL;

    % ---------------------------------------------------------------------------------------------
    % CÁLCULO DA EFICIÊNCIA DO PARQUE EÓLICO PROJETADO

    ef = P_gerada/P_livre;

    %% Definição da função objetivo

    custo = (n*((2/3)+((1/3)*exp(-0.00174*(n^2)))));
    FOB = custo/P_gerada;
end