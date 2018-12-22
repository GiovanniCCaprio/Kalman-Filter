close all, clc


%Dados
curva = 2;
if curva == 1
load('reta.mat')
elseif curva == 2
load('curva_suave.mat')
else
load('manobra.mat')
end

posicao = x_m';

% Inicializa��o
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
    R=cov_v*eye(2);
    %R=[variance_pos 0 0 0;0 variance 0 0;0 0 variance_pos 0;0 0 0 variance];
    Q = 10^2*eye(4); %Erro da velocidade)

% Planta
    F = [1 dT 0 0; 0 1 0 0; 0 0 1 dT; 0 0 0 1];
    G = [dT; (1/2)*(dT)^2 ; dT ;(1/2)*(dT)^2];
    H = [1 0 0 0; 0 0 1 0];

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
    %P = F*Q*F'; %Erro Inicial da Covariancia
    P = [4*10^5 10^4 0 0;10^4 4*10^5 0 0;0 0 4*10^5 10^4;0 0 10^4 4*10^5];   
    x = [x_m(1,1); v_x;x_m(2,1); v_y]; %estado inicial x=[x; x_dot; y; y_dot]
    xx = zeros(2,t);
% Filtro de Kalman    
for i = 1:t
    
   Mn = P*H'/(H*P*H'+R); %ganho 
   xx = [x_m(1,i);x_m(2,i)] - H*x;
  
   %calculo do estado atual
   x = x + Mn*(xx); 
   %atualizando a covari�ncia do processo
   P = (eye(4) - Mn*H)*P; %atualiza��o covariancia do erro
   errocov=[errocov; H*P*H'];%covari�ncia do erro de estima��o do estado    
   P = F*P*F' + Q;
    
   
   %Somat�rio para o EMQ
   Erro = Erro + sqrt((x(1,1)- x_c(i,1)).^2 + (x(3,1)- x_c(i,2)).^2 );
   %estimativa da velocidade
   vel=[vel; (sqrt((x(2,1)).^2 + (x(4,1)).^2))];
   
   
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

