close all, clear all, clc

% Inicializa��o
%load('reta.mat')
%load('curva_suave.mat')
load('manobra')
posicao = x_m';
media = 0; %m�dia da distribui��o gaussiana
variancia = 200;%vari�ncia da acelera��o-estimado
sigma = sqrt(variancia);%desvio padr�o
a_x = normrnd(media,sigma);%aceleracao em x com m�dia nula
a_y = normrnd(media,sigma);%aceleracao em y com m�dia nula
dT = 4; %Intervalo entre amostras
v_x = abs(posicao(1,1)-posicao(2,1))/dT;%velocidade em x
v_y = abs(posicao(1,2)-posicao(2,2))/dT;%velocidade em y
t = length(x_m(1,:)); %n�meros de amostras de observa��o
cov_v=(1200).^2; %covaria��o
sigma_pos=sqrt(cov_v);%desvio padr�o
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
plot(x_c(:,1),x_c(:,2),'r') %tra�ado real(x_c)
legend('Radar','Trajeto real do objeto')
figure(1)
plot(x_m(1,:),x_m(2,:),'-.') %tra�ado observado(x_m)
legend('Radar','Trajeto real do objeto','Trajeto observado do objeto')
% Vari�veis Filtro de Kalman
vel=[];
Erro=0; %Erro que soma a cada itera��o
errocov=[];
P = F*Q*F'; %Erro Inicial da Covariancia
D=diag(P);
P=diag(D); %zera os termos cruzados da matriz P
x = [-30000; v_x;-30000; v_y]; %estado inicial x=[x; x_dot; y; y_dot]
% Filtro de Kalman
for i = 1:t
xkp = (F*x) + [0;a_x;0;a_y];%z de kalman a cada itera��o
% B*[acel_x;acel_y];%[0;acel_x;0;acel_y];
%atualiza��o das acelera��es(modeladas como ru�dos gaussianos)
%para cada ciclo
a_x = normrnd(media,sigma);
a_y = normrnd(media,sigma);
Mn = P*H'/(H*P*H'+R); %ganho
%medi��es + respectivos ru�dos do sensor
Y = F*[posicao(i,1);xkp(2,1);posicao(i,2);xkp(4,1)];%+[normrnd(mean,sigma_pos);0;normrnd(mean,sigma_pos);0];
%calculo do estado atual
x = (xkp) + Mn*(Y - H*x);
%atualizando a covari�ncia do processo
P = (eye(4) - Mn*H)*P; %atualiza��o covariancia do erro
errocov=[errocov; H*P*H'];%covari�ncia do erro de estima��o do estado
P = F*P*F' + Q;
D=diag(P);
P=diag(D);%zera os termos cruzados da matriz P
%Somat�rio para o EMQ
Erro = Erro + sqrt((x(1,1)- x_c(i,1)).^2 + (x(3,1)- x_c(i,2)).^2 );
%estimativa da velocidade
vel=[vel; (sqrt((xkp(2,1)).^2 + (xkp(4,1)).^2))];
figure(1)
plot(x(1,1), x(3,1),'^k');
legend('Radar','Trajeto real do objeto','Trajeto observado do objeto', 'Estimativa da posi��o')
end
%Erro M�dio Quadr�tico
Erro_MQ = Erro/t
% Velocidade M�dia Estimada
V_media = sum(vel)/(max(size(vel)))
% Variancia da velocidade M�dia
V_var = var(vel)
% Plotando a covari�ncia do erro da posi��o
figure(3)
title('Covari�ncia do erro da posi��o')
ylabel('Covari�ncia do erro')
grid on
hold on
errocov_pos=[];
for i=1:4:(4*t)
errocov_pos = [errocov_pos;errocov(i,1)];
end
figure(3)
plot(errocov_pos)