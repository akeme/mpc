%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MPC FOR MOBILE ROBOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

Ts = 0.1;             % sampling time
Tsim = 85
h = 0.1;           % integrating time


x0 = [0; 0; 0];  % x,y,theta
X = x0;

% Runge kutta

% Y(n+1) = Y(n) + (1/6)(k1 +2 k2 +2 k3 +k4);
% k1 = h f(X(n),Y(n))
% k2 = h f(X(n)+h/2, Y(n)+k1/2)
% k3 = h f(X(n)+h/2, Y(n) +k2/2)
% k4 = h f(X(n)+h, Y(n)+k3)

% The next loop difines the reference

for k = 1:round(Tsim/h)
    v(k) = 0.3;
    if k<= 419;
    w(k) = 0.15;
    else
        w(k) = -0.15;
    end
 
 k1 = h*[v(k)*cos(x0(3)); v(k)*sin(x0(3));w(k)];
 k2 = h*[v(k)*cos(x0(3)+k1(3)/2); v(k)*sin(x0(3)+k1(3)/2);w(k)];
 k3 = h*[v(k)*cos(x0(3)+k2(3)/2); v(k)*sin(x0(3)+k2(3)/2);w(k)];
 k4 = h*[v(k)*cos(x0(3)+k3(3)); v(k)*sin(x0(3)+k3(3));w(k)];
 
 X = [X, X(:,k) + (1/6)*(k1+2*k2+2*k3+k4)];
 x0 = X(:,k+1);
% 
% X1 = [X1, X1(:,k)+ h*([v(k)*cos(xa(3));v(k)*sin(xa(3));w(k)])];
%  xa = X1(:,k+1);
%  
end

Xr = X;
Ur = [v;w];
 plot(X(1,:),X(2,:),'r--')
 hold on
% plot(X1(1,:),X1(2,:),':')
% legend('runge_kutta','euler')
grid

vr = v;
wr = w;

tic
N_mpc = 5;
Q_mpc = [1 0 0;
     0 1 0;
     0 0 .1];% no segundo exemplo foi usado 0.1
 R_mpc = 0.1*eye(2);

  vmin_mpc = -0.4;
 vmax_mpc = 0.4;
 wmin_mpc = -0.4;
 wmax_mpc = 0.4;
% %relaxando as restrições / (retirando-as)
% vmin = -10;
%  vmax = 10;
%  wmin = -10;
%  wmax = 10;
 
 lb_mpc = [];
 ub_mpc = [];
 for k = 1:N_mpc
     Qb_mpc(3*k-2:3*k,3*k-2:3*k) = Q_mpc;
     Rb_mpc(2*k-1:2*k,2*k-1:2*k)=R_mpc;
     lb_mpc = [lb_mpc;[vmin_mpc;wmin_mpc]];
     ub_mpc = [ub_mpc;[vmax_mpc;wmax_mpc]]
 end



x0 = [0; -1; pi/2];  % x,y,theta
x0_mpc = x0;
X_mpc = x0_mpc;
v_mpc = [];
w_mpc = [];

for k = 1:round(Tsim/h)-N_mpc;
    A_mpc = [];
       xb_mpc = [];
   aux = eye(3); 
   for j = 0:N_mpc-1
       Aj_mpc{j+1} = [1 0 -Ts*vr(k+j)*sin(Xr(3,k+j));
               0 1  Ts*vr(k+j)*cos(Xr(3,k+j));
               0 0  1];
       Bj_mpc{j+1} = [cos(Xr(3,k+j))*Ts  0;
             sin(Xr(3,k+j))*Ts  0;
             0                Ts];
         aux = Aj_mpc{j+1}*aux;
          A_mpc = [A_mpc;aux];
         xb_mpc = [xb_mpc;Ur(:,k+j)];
   end
   
   B_mpc = [];

   for j = 1:N_mpc
       aux2  = eye(3);
       Baux = [];
       for i = 1:j
           if i<j
           Baux = [aux2*Bj_mpc{j-i+1},Baux];
           aux2 = aux2*Aj_mpc{j-i+1};
           else
              Baux = [aux2*Bj_mpc{j-i+1},Baux,zeros(3,(N_mpc-j)*2)]; 
           end
       end
       B_mpc = [B_mpc;Baux];
       
   end
   
   H_mpc = 2*(B_mpc'*Qb_mpc*B_mpc+Rb_mpc);
   f_mpc = 2*B_mpc'*Qb_mpc*A_mpc*(X_mpc(:,k)-Xr(:,k));
   options = optimset('Display','off','Largescale','off');
   %Uopt = quadprog(H,f,[],[],[],[],[lb-xb],[ub-xb],[],options);
   Uopt_mpc = -inv(H_mpc)*f_mpc;
   U_mpc(:,k) = Ur(:,k)+Uopt_mpc(1:2,1);
  v_r(k) = Ur(1,k);
    w_r(k) = Ur(2,k);
   v_mpc(k) = U_mpc(1,k);
   w_mpc(k) = U_mpc(2,k);
   
 k1 = h*[v_mpc(k)*cos(x0_mpc(3)); v_mpc(k)*sin(x0_mpc(3));w_mpc(k)];
 k2 = h*[v_mpc(k)*cos(x0_mpc(3)+k1(3)/2); v_mpc(k)*sin(x0_mpc(3)+k1(3)/2);w_mpc(k)];
 k3 = h*[v_mpc(k)*cos(x0_mpc(3)+k2(3)/2); v_mpc(k)*sin(x0_mpc(3)+k2(3)/2);w_mpc(k)];
 k4 = h*[v_mpc(k)*cos(x0_mpc(3)+k3(3)); v_mpc(k)*sin(x0_mpc(3)+k3(3));w_mpc(k)];
 
 X_mpc = [X_mpc, X_mpc(:,k) + (1/6)*(k1+2*k2+2*k3+k4)];
 x0_mpc = X_mpc(:,k+1);
   k
   
   Erro_mpc(:,k) = Xr(:,k+1)- X_mpc(:,k+1); 
   E1_qmpc(k) = Erro_mpc(1,k)*Erro_mpc(1,k);
   E2_qmpc(k) = Erro_mpc(2,k)*Erro_mpc(2,k);
   E3_qmpc(k) = Erro_mpc(3,k)*Erro_mpc(3,k);
   
    E1_tmpc(k) = k*abs(Erro_mpc(1,k));
    E2_tmpc(k) = k*abs(Erro_mpc(2,k));
    E3_tmpc(k) = k*abs(Erro_mpc(3,k));
    
    E1_tqmpc(k) = k*E1_qmpc(k);
    E2_tqmpc(k) = k*E2_qmpc(k);
    E3_tqmpc(k) = k*E3_qmpc(k);
end

tempo_mpc = toc

 plot(X_mpc(1,:),X_mpc(2,:),'b')
 legend('referência','trajetória do robô')
title('MPC Linearizado')
xlabel('x[m]')
ylabel('y[m]')

 figure 
 ne = length(Erro_mpc);
 tempo = Ts*(0:ne-1);
 
 plot(tempo,v_mpc)
 hold on
plot(tempo,v_r,'b--')
plot (tempo,w_mpc,'r')
plot(tempo,w_r,'b--')
legend('v','w','vr','wr')
title('Velocidades MPC linearizado')
xlabel('t[s]')
ylabel('v[m/s] w[rad/s]')

 
 figure

 E1_mpc = Erro_mpc(1,:);
 E2_mpc = Erro_mpc(2,:);
 E3_mpc = Erro_mpc(3,:);

plot(tempo,E1_mpc,tempo,E2_mpc,tempo,E3_mpc);

ERRO_TOTAL_mpc = E1_mpc*E1_mpc'*Q_mpc(1,1)+E2_mpc*E2_mpc'*Q_mpc(2,2)+E3_mpc*E3_mpc'*Q_mpc(3,3);
legend('x','y','theta')
title('Erro MPC linearizado')
xlabel('t[s]')
ylabel('x[m] y[m] theta x[rad]')

grid

 
commandwindow

%indices de desempenho

% Integral do erro absoluto
IAE1_mpc = trapz(abs(E1_mpc))
IAE2_mpc = trapz(abs(E2_mpc))
IAE3_mpc = trapz(abs(E3_mpc))

%  IE1 = sum(E1)/length(E1); %essas contas estão certas?
%  IE2 = sum(E2)/length(E2);
%  IE3 = sum(E3)/length(E3);

% Integral do erro ao quadrado

ISE1_mpc = trapz(E1_qmpc)
ISE2_mpc = trapz(E2_qmpc)
ISE3_mpc = trapz(E3_qmpc)


ITAE1_mpc = trapz(E1_tmpc)
ITAE2_mpc = trapz(E2_tmpc)
ITAE3_mpc = trapz(E3_tmpc)


ITSE1_mpc = trapz(E1_tqmpc)
ITSE2_mpc = trapz(E2_tqmpc)
ITSE3_mpc = trapz(E3_tqmpc)

% desvio1 = std(E1);
% desvio2 = std(E2);
% desvio3 = std(E3);

% erro_max1 = max(abs(E1));
% erro_max2 = max(abs(E2));
% erro_max3 = max(abs(E3));


