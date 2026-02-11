clc
clear all
close all
% ----------------------------------------------------------------------------------------------
% ENTRADA DE DADOS NO PROGRAMA

% Parâmetros do modelo

linhas = 10; % Número de linhas da matriz que representa o terreno
colunas = 10; % Número de colunas da matriz que representa o terreno
nmax = 100; % Número máximo de turbinas
alt_rotor = 60; % Altura do rotor [m]
diam_rotor = 40; % Diâmetro do rotor [m]
coef_empuxo = 0.88; % Coeficiente de empuxo
f_axial = 0.3268; % Fator de indução axial
rugosidade = 0.30; % Rugosidade do solo
vel_vento_livre = 12; % Velocidade do vento livre [m/s]
arraste = (0.5/log(alt_rotor/rugosidade)); 
raio_esteira = 20*(sqrt((1-f_axial)/(1-(2*f_axial)))); % Raio da esteira

% Probabilidade de ocorrência do vento em direções com intervalo de 45º.

%dir:   0º    45º   90º  135º  180º  225º  270º  315º
  fk = [0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125];
% fk = [0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000];
% fk = [1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000];

% O vetor fk deve ser preenchido obrigatoriamente com 8 posições.
% A soma de todas as posições deve totalizar 1.

% Parâmetros do WOA
N_agentes = 1500; % Número de agentes de busca - 1500
lim_iteracoes = 120; % Limite de iterações - 120

% ----------------------------------------------------------------------------------------------
% VALIDAÇÃO DE DADOS

if length(fk)<8 || sum(fk) ~= 1
    errordlg('A densidade de probabilidade informada é inválida.');
    error('A densidade de probabilidade informada é inválida.');
end
if nmax>(linhas*colunas)
    errordlg('O número de aerogeradores informados é incompatível com o terreno escolhido.');
    error('O número de aerogeradores informados é incompatível com o terreno escolhido.');
end

% ----------------------------------------------------------------------------------------------
% CRIAÇÃO DA POPULAÇÃO INICIAL

for i = 1:N_agentes
    Matriz_eolica = zeros(linhas,colunas);
    % Sorteio das posições onde haverão turbinas
    coef_linhas=zeros(nmax,1);
    coef_colunas=zeros(nmax,1);
    for k=1:nmax;
      coef_linhas(k,1) = randi(linhas,1); % vetor de índice de linhas
      coef_colunas(k,1) = randi(colunas,1); % vetor de índice de colunas
    end
    coef = [coef_linhas,coef_colunas]; % matriz de índice das células que receberão turbinas
    % Implantação das turbinas na matriz que representa o terreno
    for k=1:length(coef(:,1));
           Matriz_eolica(coef(k,1),coef(k,2))=1;
    end
    Ag_busca{i}=Matriz_eolica;
    % Cálculo da quantidade real de turbinas implantadas no terreno
    N_turb(i)=sum(sum(Ag_busca{i}));
end
% ----------------------------------------------------------------------------------------------
% OBTENÇÃO DE PARÂMETROS PARA CADA UMA DAS MATRIZES GERADAS

FOB = zeros(1,N_agentes);
ef = zeros(1,N_agentes);
P_gerada = zeros(1,N_agentes);
for k = 1:N_agentes
    M = Ag_busca{k};
    n = N_turb(k);
    [FOB(k),ef(k),P_gerada(k)] = F_objetivo(M,n,fk,linhas,colunas,alt_rotor,diam_rotor,coef_empuxo,f_axial,rugosidade,vel_vento_livre,arraste,raio_esteira);
end
% ----------------------------------------------------------------------------------------------
% OBTENÇÃO DA MELHOR SOLUÇÃO ATUAL

FOB_min = min(FOB);
for k = 1:N_agentes
    if FOB(k) == FOB_min;
        pos_melhor_agente = k;
    end
end
melhor_AB = Ag_busca{pos_melhor_agente};

% ----------------------------------------------------------------------------------------------
% PROCESSAMENTO DO BWOA
it_atual = 2;
cte_a = zeros(1,lim_iteracoes);
coef_A = zeros(1,lim_iteracoes);
coef_C = zeros(1,lim_iteracoes);
while it_atual <= lim_iteracoes
    % No laço for abaixo, cada um dos agentes de busca será atualizado de acordo com o WOA
    i = 1;
    while i <= N_agentes
        % Definição dos parâmetros iniciais para o agente de busca i.
        cte_a(it_atual) = 2-(it_atual*((2)/lim_iteracoes)); % O valor de a decresce linearmente de 2 até 0
        cte_r1 = rand;
        cte_r2 = rand;
        cte_p = rand;
        cte_l = (2*rand)-1; % O valor de l varia entre -1 e 1
        cte_b = 1; % Valor arbitrário
        coef_A(it_atual) = (2*cte_a(it_atual)*cte_r1)-cte_a(it_atual);
        coef_C(it_atual) = 2*cte_r2;
        % Escolha da equação de atualização da posição do agente de busca i
        % de acordo com os parâmetros iniciais
        if cte_p < 0.5
            if abs(coef_A(it_atual)) < 1 % Mecanismo de encolhimento do cerco
                D = abs(coef_C(it_atual) * melhor_AB - Ag_busca{i});
            else if abs(coef_A(it_atual)) >= 1 % Fase de exploração
                    AB = randi([1,N_agentes]);
                    AB_aleatorio = Ag_busca{AB};
                    D = abs(coef_C(it_atual) * AB_aleatorio - Ag_busca{i});
                end
            end
        else if cte_p >= 0.5 % Atualização de posições pelo modelo espiral
                D = abs(melhor_AB - Ag_busca{i});
            end
        end
        % Aplica a função de transferência e atualiza o agente de busca
        Cstep = (1./(1+exp(-10*(coef_A(it_atual)*(D)-0.5)))); %FT
        M = Ag_busca{i};
        for l = 1:linhas
            for c = 1:colunas
                if rand < Cstep(l,c)
                    M(l,c) = ~M(l,c);
                else
                    M(l,c) = M(l,c);
                end
            end
        end
        Ag_busca{i} = M;
        % Verifica se o agente de busca gerado está dentro do limite de
        % aerogeradores
          if sum(sum(Ag_busca{i})) <= nmax
              i = i+1;
          end
    end
    % Calcula a função objetivo da nova geração
    for k = 1:N_agentes
        M = Ag_busca{k};
        n = sum(sum(M));
        [FOB(k),ef(k),P_gerada(k)] = F_objetivo(M,n,fk,linhas,colunas,alt_rotor,diam_rotor,coef_empuxo,f_axial,rugosidade,vel_vento_livre,arraste,raio_esteira);
    end
    % Encontra o melhor agente de busca da geração atual
    FOB_min_G_atual = min(FOB);
    for k = 1:N_agentes
        if FOB(k) == FOB_min_G_atual;
            pos_melhor_agente_G_atual = k; % Índice do melhor agente de busca da geração
        end
    end
    % Atualiza o melhor AB se houver um melhor
    if FOB_min_G_atual <= FOB_min(it_atual-1)
        FOB_min(it_atual) = FOB_min_G_atual;
        melhor_AB = Ag_busca{pos_melhor_agente_G_atual};
    else
        FOB_min(it_atual) = FOB_min(it_atual-1);
    end
    % Atualiza a iteração atual
    it_atual = it_atual+1;
    % Salva o melhor agente de busca da iteração
    melhor_AB_cell{it_atual} = melhor_AB;
end
% Escolhe o melhor agente de busca gerado
FOB_otimo = min(FOB_min);
for k = 1:lim_iteracoes
    if FOB_min(k) == FOB_otimo
        pos_otimo = k;
    end
end
config_otima = melhor_AB_cell{pos_otimo};
% plot da convergência do WOA
plot(1:lim_iteracoes,FOB_min)
title('Convergência do WOA')
xlabel('Iteração')
ylabel('F. objetivo')

% -----------------------------------------------------------------------------------------------
% Apresentação dos resultados
for k = 0:9
   for m = 0:9
       if config_otima(m+1,k+1) == 1
            subplot(1,2,1)
            plot(k+0.5,9-m+0.5,'bo','MarkerFaceColor', 'b');
            title('Configuração ótima')
            hold on
            grid on
       end
   end
end

% plot da convergência do WOA
subplot(1,2,2)
plot(1:lim_iteracoes,FOB_min,'b','LineWidth', 1.2)
title('Convergência do WOA');
xlim([1,lim_iteracoes])
xlabel('Iteração')
ylabel('F. objetivo')
fig = gcf; % Obtém o identificador da figura atual
set(gcf, 'Position', get(0,'Screensize'))

% ---------------------------------------------------------------------------------------------------
% VISUALIZAÇÃO DOS PARÂMETROS INTERNOS DO ALGORITMO

% figure(2)
% subplot(1,3,1)
% plot((1:lim_iteracoes),cte_a,'b')
% title('Constante a')
% subplot(1,3,2)
% plot((1:lim_iteracoes),coef_A,'b')
% title('Coeficiente A')
% subplot(1,3,3)
% plot((1:lim_iteracoes),coef_C,'b')
% title('Coeficiente C')