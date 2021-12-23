% Equipe
% JOÃO GUILHERME SALES DE OLIVEIRA (20190034570)
% GABRYEL MARTINS RAPOSO DE ALENCAR (20190029338 )

% limpando o console e as variáveis anteriores
clear
clc

% utilizaremos variáveis simbólicas para cálculos de parâmetros genéricos
% com o tempo
syms t

% lendo o arquivo de entrada
entradaDados = readmatrix('entrada_circ_Q3.txt')

%lendo a frequência de oscilação da fonte
w = max(entradaDados(:,13))

% vamos obter o número máximo de nós pegando o máixmo valor da coluna de
% nós correspondente. O número de ramos é dado pela quantidade de linhas da
% matriz de entrada
matrizNos = entradaDados(:,2:3)
numNos = max(max(matrizNos))
numRamos = max(max(entradaDados(:,1)))

% matriz de incidência completa
matrizIncidenciaCompleta = zeros(numNos, numRamos)

for r = 1:numRamos
    no_saida = matrizNos(r,1)
    no_entrada = matrizNos(r,2)
    %percorrermos a matriz dos nós e pegamos os nós de entrada/saída
    %referente aos ramos para atribuir os valores da nossa matriz de
    %incidência
    matrizIncidenciaCompleta(no_saida,r) = 1
    matrizIncidenciaCompleta(no_entrada,r) = -1
end

A = matrizIncidenciaCompleta(1:numNos-1,:) 

A_T = transpose(A)

% Levando o circuito para laplace

%resistor
R = double(entradaDados(:, 4))
%indutor
L = double(entradaDados(:,5))
Xl = (i*w) * L
%capacitância
C = double(entradaDados(:,6))
Xc = 1./(C*(w*i))
Xc(Xc(:) == Inf) = 0 %tratando a exceção de divisão por zero e adicionando...
% zero caso o elemento passado seja nulo.

%impedância equivalente dos ramos
Zeq = R + Xl + Xc
Yeq = 1./Zeq(:)
Yeq(Yeq(:) == Inf) = 1

%%% Vamos criar a matriz de admitância de ramo (é um matriz diagonal)
Yb = diag(Yeq)
% assim conseguimos calcular Fn

% Uma vez calculada a matriz de admitância de Ramo, precisa-se separar os vetores de 
% fontes independentes de tensão e de corrente de modo a termos informações
% para o calcular de Is.
% OBS.: iremos considerar as condições iniciais

vs = zeros(numRamos,1);
Js = zeros(numRamos,1);

for r = 1:numRamos
    vs(r) = (entradaDados(r,7)*cos(entradaDados(r,8)) + entradaDados(r,7)*i*sin(entradaDados(r,8)))
    Js(r) =  (entradaDados(r,10)*cos(entradaDados(r,11)) + entradaDados(r,10)*i*sin(entradaDados(r,11)))
end

% calculando a matriz de admitancia de nó
Yn = A * Yb * A_T
Is = (A*Yb)*vs - (A*Js)

% Fasores das Tensões de Nó
E = (inv(Yn)) * Is

%Fasores da Tensão de ramo
V = A_T * E

%Fasores de Corrente nos Ramos:
J = Js + Yb*V - Yb*vs

%nesse ponto, precisamos calcular a função que iremos plotar, obtendo o
%módulo e fase de cada fasor obtido

moduloE = abs(E)
faseE = angle(E)

moduloV = abs(V)
faseV = angle(V)

moduloJ = abs(J)
faseJ = angle(J)

% %gráficos
tmax = 20
passo = 1e-2
t = 0:passo:tmax 
% %e_t = subs(e,t)
% %v_t = subs(v,t)
% %j_t = subs(j,t)

%tensão no capacitor de (1/5)F
v2_t = moduloV(1) * cos(w*t + faseV(1))

%gráfico da tensão no capacitor de 1/5 F 
figure
plot(t,v2_t)
axis([0 tmax -60 60])
set(gca,'FontSize',16)
xlabel('tempo (s)','Interpreter','LaTex','FontSize',18)
ylabel('Tensão (V)','Interpreter','LaTex','FontSize',18)
grid()
title ('Tensão no capacitor de 1/5 Faraday')