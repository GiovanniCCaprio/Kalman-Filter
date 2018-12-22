close all, clear all, clc

% Inicialização
%load('reta.mat')
%load('curva_suave.mat')
load('manobra')
posicao = x_m';
media = 0; %média da distribuição gaussiana
variancia = 200;%variância da aceleração-estimado
sigma = sqrt(variancia);%desvio padrão
a_x = normrnd(media,sigma);%aceleracao em x com média nula
a_y = normrnd(media,sigma);%aceleracao em y com média nula
dT = 4; %Intervalo entre amostras
v_x = abs(posicao(1,1)-posicao(2,1))/dT;%velocidade em x
v_y = abs(posicao(1,2)-posicao(2,2))/dT;%velocidade em y
t = length(x_m(1,:)); %números de amostras de observação
cov_v=(1200).^2; %covariação
sigma_pos=sqrt(cov_v);%desvio padrão
R=[cov_v 0 0 0;0 variancia 0 0;0 0 cov_v 0;0 0 0 variancia];
%R=[variance_pos 0 0 0;0 variance 0 0;0 0 variance_pos 0;0 0 0 variance];
Q = [600^2 0 0 0;0 50^2 0 0;0 0 600^2 0;0 0 0 50^2]; %Erro da velocidade)
% Planta
F = [1 dT 0 0; 0 1 0 0; 0 0 1 dT; 0 0 0 1];
G = [(1/2)*(dT)^2 0;dT 0; 0 (1/2)*(dT)^2; 0 dT];
H =[1 0 0 0; 0 1 0 0;0 0 1 0; 0 0 0 1];
figure(1)
plot(0,0,'or'); %Posicionamento do radar em coordenadas
legend('Radar')
xlabel('Coordenada X(m)'), ylabel('Coordenada Y(m)')
grid on
hold on
figure(1)
plot(x_c(:,1),x_c(:,2),'r') %traçado real(x_c)
legend('Radar','Trajeto real do objeto')
figure(1)
plot(x_m(1,:),x_m(2,:),'-.') %traçado observado(x_m)
legend('Radar','Trajeto real do objeto','Trajeto observado do objeto')
% Variáveis Filtro de Kalman
vel=[];
Erro=0; %Erro que soma a cada iteração
errocov=[];
P = F*Q*F'; %Erro Inicial da Covariancia
D=diag(P);
P=diag(D); %zera os termos cruzados da matriz P
x = [-30000; v_x;-30000; v_y]; %estado inicial x=[x; x_dot; y; y_dot]
% Filtro de Kalman
for i = 1:t
xkp = (F*x) + [0;a_x;0;a_y];%z de kalman a cada iteração
% B*[acel_x;acel_y];%[0;acel_x;0;acel_y];
%atualização das acelerações(modeladas como ruídos gaussianos)
%para cada ciclo
a_x = normrnd(media,sigma);
a_y = normrnd(media,sigma);
Mn = P*H'/(H*P*H'+R); %ganho
%medições + respectivos ruídos do sensor
Y = F*[posicao(i,1);xkp(2,1);posicao(i,2);xkp(4,1)];%+[normrnd(mean,sigma_pos);0;normrnd(mean,sigma_pos);0];
%calculo do estado atual
x = (xkp) + Mn*(Y - H*x);
%atualizando a covariância do processo
P = (eye(4) - Mn*H)*P; %atualização covariancia do erro
errocov=[errocov; H*P*H'];%covariância do erro de estimação do estado
P = F*P*F' + Q;
D=diag(P);
P=diag(D);%zera os termos cruzados da matriz P
%Somatório para o EMQ
Erro = Erro + sqrt((x(1,1)- x_c(i,1)).^2 + (x(3,1)- x_c(i,2)).^2 );
%estimativa da velocidade
vel=[vel; (sqrt((xkp(2,1)).^2 + (xkp(4,1)).^2))];
figure(1)
plot(x(1,1), x(3,1),'^k');
legend('Radar','Trajeto real do objeto','Trajeto observado do objeto', 'Estimativa da posição')
end
%Erro Médio Quadrático
Erro_MQ = Erro/t
% Velocidade Média Estimada
V_media = sum(vel)/(max(size(vel)))
% Variancia da velocidade Média
V_var = var(vel)
% Plotando a covariância do erro da posição
figure(3)
title('Covariância do erro da posição')
ylabel('Covariância do erro')
grid on
hold on
errocov_pos=[];
for i=1:4:(4*t)
errocov_pos = [errocov_pos;errocov(i,1)];
end
figure(3)
plot(errocov_pos)