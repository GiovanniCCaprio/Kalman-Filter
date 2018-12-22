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

% Inicialização
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
plot(x_c(:,1),x_c(:,2),'r') %traçado real(x_c)
legend('Radar','Trajeto real do objeto')

figure(1)
plot(x_m(1,:),x_m(2,:),'-.') %traçado observado(x_m)
legend('Radar','Trajeto real do objeto','Trajeto observado do objeto')

% Variáveis Filtro de Kalman
    vel=[];
    Erro=0; %Erro que soma a cada iteração
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
   %atualizando a covariância do processo
   P = (eye(4) - Mn*H)*P; %atualização covariancia do erro
   errocov=[errocov; H*P*H'];%covariância do erro de estimação do estado    
   P = F*P*F' + Q;
    
   
   %Somatório para o EMQ
   Erro = Erro + sqrt((x(1,1)- x_c(i,1)).^2 + (x(3,1)- x_c(i,2)).^2 );
   %estimativa da velocidade
   vel=[vel; (sqrt((x(2,1)).^2 + (x(4,1)).^2))];
   
   
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

