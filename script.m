% Equipe:
% João Guilherme Sales de Oliveira (20190034570)
% Gabriel Reposo de ALencar...

clear

% utilizaremos variáveis simbólicas para cálculos de parâmetros genéricos e
% com laplace
syms s
syms t

% lendo o arquivo de entrada
entradaDados = readmatrix('data_in.txt')

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
Xl = s * L
%capacitância
C = double(entradaDados(:,7))
Xc = 1./(s*C)
Xc(Xc(:) == Inf) = 0 %tratando a exceção de divisão por zero e adicionando...
% zero caso o elemento passado seja nulo.

%impedância equivalente dos ramos
Zeq = R + Xl + Xc
Yeq = 1./Zeq(:)

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
    vs(r) = entradaDados(r,9) + entradaDados(r,5)*entradaDados(r,6)  %%% validar  derminamos para cada ramo por isso o (r)
    

    Js(r) = entradaDados(r,10) + entradaDados(r,6)*entradaDados(r,7)  %%% validar
end

vs = laplace(sym(vs))
Js = laplace(sym(Js))


% calculando a matriz de admitancia de nó
Yn = A * Yb * A_T
Is = A*Yb*vs - A*Js


% Cálculo das Tensões de Nó
E = inv(Yn) * Is

%Cálculo da Tensão de ramo
V = A_T * E

%Corrente nos Ramo:
J = Js + Yb*V - Yb*vs


% Nesse ponto, temos as respostas do circuito em Laplace, portanto,
% precisamos convertê-las para o Tempo


%corrente nos ramos no tempo:
j = ilaplace(J)

%tensão nos ramos no tempo:
v = ilaplace(V)

%tensão nos nós no tempo:
e = ilaplace(E)




%gráficos
tmax = 40e-3
passo = 1e-4 
t = 0:passo:tmax 
%e_t = subs(e,t)
%v_t = subs(v,t)
%j_t = subs(j,t)

%tensão no resistor de 10 ohms
v0 = 10*j(3)
v0_t = subs(v0, t)

%gráfico da tensão no resistor de 10 ohm, desejado
figure
plot(t,v0_t)
axis([0 tmax -60 60])
set(gca,'FontSize',16)
xlabel('tempo (s)','Interpreter','LaTex','FontSize',18)
ylabel('Tensão (V)','Interpreter','LaTex','FontSize',18)
grid()
title ('Tensão no resistor de 10 ohms')















    









   







    



