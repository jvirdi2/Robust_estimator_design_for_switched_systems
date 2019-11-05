% Robust design for switched systems with unknown switching time
% Implementation by Jasvir Virdi
% LMIs are solved using YALMIP (SDPT3 solver)

% The code is divided into sections namely: Plant dynamics, computing
% observer gains, simulation of plant dynamics with unknown inputs, fast
% estimator simulation, slow estimator simulation and finally the hybrid
% estimator which is highlighted in the paper. 

% The code has to be run in the following order:1) Plant dynamics, 2)
% Computing observer gains, 3)Simulation of plant dynamics with unknown
% 4) Either fast estimator or slow estimator or hybrid estimator whichever
% you wish to run. 

% One can observe that fast estimator has a peaking phenomenon issue and slow
% estimator doesn't work well for faster switching. Hybrid Estimator gives
% the best of both the worlds where there is no peaking phenomenon and
% smaller switching times can be used.

%% Plant dynamics
A1=[[-1 2 2];[0 -2 1];[-1 0 -3]];
B1=[0;0;1];
E1=[0;1;0];
C1=[[1 0 0];[0 1 0]];

A2=[[-2 1 0];[-3 -1 1];[1 -2 -1]];
B2=[1;0;0];
E2=[-1;0;0];
C2=[[1 1 0];[1 0 1]];
phi=2*eye(3);


U1_paper=-E1*pinv(C1*E1);
V1_paper=eye(2)-(C1*E1)*pinv(C1*E1);

U2_paper=-E2*pinv(C2*E2);
V2_paper=eye(2)-(C2*E2)*pinv(C2*E2);

%% Computing observer gains
P=sdpvar(3,3);
Ktilda_1= sdpvar(3,2,'full','real');
Ytilda_1= sdpvar(3,2,'full','real');
Ktilda_2= sdpvar(3,2,'full','real');
Ytilda_2= sdpvar(3,2,'full','real');

% (alpha=500)

F = [P >= 0, ...
     P*A1+ A1'*P+ P*U1_paper*C1*A1+ A1'*C1'*U1_paper'*P+ Ytilda_1*V1_paper*C1*A1+ A1'*C1'*V1_paper'*Ytilda_1'- Ktilda_1*C1- C1'*Ktilda_1'+2*500*P <= 0, ...
     phi*P*phi-5*P<=0,...
     P*A2+ A2'*P+ P*U2_paper*C2*A2+ A2'*C2'*U2_paper'*P+ Ytilda_2*V2_paper*C2*A2+ A2'*C2'*V2_paper'*Ytilda_2'- Ktilda_2*C2- C2'*Ktilda_2'+ 2*500*P<= 0, ...
     P*A1+ A1'*P+ P*U1_paper*C1*A1+ A1'*C1'*U1_paper'*P+ Ytilda_1*V1_paper*C1*A1+ A1'*C1'*V1_paper'*Ytilda_1'- Ktilda_1*C1- C1'*Ktilda_1'+2*0*P <= 0, ...
     P*A2+ A2'*P+ P*U2_paper*C2*A2+ A2'*C2'*U2_paper'*P+ Ytilda_2*V2_paper*C2*A2+ A2'*C2'*V2_paper'*Ytilda_2'- Ktilda_2*C2- C2'*Ktilda_2'+ 2*0*P<= 0, ...
     ]; 

optimize(F);
Pval=value(P);
Ktilda1_val=value(Ktilda_1);
Ytilda1_val=value(Ytilda_1);
K1=inv(Pval)*Ktilda1_val;
Y1=inv(Pval)*Ytilda1_val;
Y1(:,2)=zeros(3,1);

Ktilda2_val=value(Ktilda_2);
Ytilda2_val=value(Ytilda_2);
K2=inv(Pval)*Ktilda2_val;
Y2=inv(Pval)*Ytilda2_val;


J1_m=U1_paper+Y1*V1_paper;
M1_m=eye(3)+J1_m*C1;
H1_m=M1_m*A1-K1*C1;
L1_m=K1*(eye(2,2)+C1*J1_m)-M1_m*A1*J1_m;
G1_m=M1_m*B1;

J2_m=U2_paper+Y2*V2_paper;
M2_m=eye(3)+J2_m*C2;
H2_m=M2_m*A2-K2*C2;
L2_m=K2*(eye(2,2)+C2*J2_m)-M2_m*A2*J2_m;
G2_m=M2_m*B2;

% alpha=1
F = [P >= 0, ...
     P*A1+ A1'*P+ P*U1_paper*C1*A1+ A1'*C1'*U1_paper'*P+ Ytilda_1*V1_paper*C1*A1+ A1'*C1'*V1_paper'*Ytilda_1'- Ktilda_1*C1- C1'*Ktilda_1'+2*1*P <= 0, ...
     phi*P*phi-5*P<=0,...
     P*A2+ A2'*P+ P*U2_paper*C2*A2+ A2'*C2'*U2_paper'*P+ Ytilda_2*V2_paper*C2*A2+ A2'*C2'*V2_paper'*Ytilda_2'- Ktilda_2*C2- C2'*Ktilda_2'+ 2*1*P<= 0, ...
     P*A1+ A1'*P+ P*U1_paper*C1*A1+ A1'*C1'*U1_paper'*P+ Ytilda_1*V1_paper*C1*A1+ A1'*C1'*V1_paper'*Ytilda_1'- Ktilda_1*C1- C1'*Ktilda_1'+2*0*P <= 0, ...
     P*A2+ A2'*P+ P*U2_paper*C2*A2+ A2'*C2'*U2_paper'*P+ Ytilda_2*V2_paper*C2*A2+ A2'*C2'*V2_paper'*Ytilda_2'- Ktilda_2*C2- C2'*Ktilda_2'+ 2*0*P<= 0, ...
     ]; 
 
optimize(F);
Pval_c=value(P);
Ktilda1_val_c=value(Ktilda_1);
Ytilda1_val_c=value(Ytilda_1);
K1_c=inv(Pval_c)*Ktilda1_val_c;
Y1_c=inv(Pval_c)*Ytilda1_val_c;
Y1_c(:,2)=zeros(3,1);

Ktilda2_val_c=value(Ktilda_2);
Ytilda2_val_c=value(Ytilda_2);
K2_c=inv(Pval_c)*Ktilda2_val_c;
Y2_c=inv(Pval_c)*Ytilda2_val_c;

J1_c=U1_paper+Y1_c*V1_paper;
M1_c=eye(3)+J1_c*C1;
H1_c=M1_c*A1-K1_c*C1;
L1_c=K1_c*(eye(2,2)+C1*J1_c)-M1_c*A1*J1_c;
G1_c=M1_c*B1;

J2_c=U2_paper+Y2_c*V2_paper;
M2_c=eye(3)+J2_c*C2;
H2_c=M2_c*A2-K2_c*C2;
L2_c=K2_c*(eye(2,2)+C2*J2_c)-M2_c*A2*J2_c;
G2_c=M2_c*B2;

%% Simulation of plant dynamics with unknown inputs 
x_initial=[1;1;1];
tf=50;
dt=0.0002;
actual_state=[x_initial,zeros(3,tf/dt)];
actual_output=[[1;1],zeros(2,tf/dt)];
time_for_switch=2.5; 
% tf=50 time_for_switch=2.5, hybrid works, fast (peaking but works), slow doesn't work well
% need design parameter tuning
% tf=300 time_for_switch=50, hybrid, fast (peaking but works) and slow all work
% tf=50  time_for_switch=0.5  divergence (hybrid), fast estimator only (peaking but works)
time_for_u_step=0.5; 
control_vec=zeros(1,tf/dt+1);
mode_vec=zeros(1,tf/dt+1);
counter1=1;
counter2=1;
for index=1:((tf/dt))
    t=(index-1)*dt;
    v=sin(t); % Unknown input
    if mod(floor(t/time_for_switch),2)==0
        mode=1;
    else
        mode=2;
    end
    mode_vec(index)=mode;
    if mod(floor(t/time_for_u_step),2)==0
        u=0.5;
    else
        u=0;
    end
    control_vec(index)=u;
    
    if mode==1
        if t<=time_for_switch
            actual_state(:,index+1)=actual_state(:,index)+dt*(A1*actual_state(:,index)+B1*u+E1*v);
            actual_output(:,index+1)=C1*actual_state(:,index);
        else
            if counter1==1
                next_state=phi*actual_state(:,index);
                counter1=0;
                counter2=1;
            else
                next_state=actual_state(:,index);
            end
            actual_state(:,index+1)=next_state+dt*(A1*next_state+B1*u+E1*v);
            actual_output(:,index+1)=C1*next_state;
        end
    else
        if counter2==1
            next_state_2=phi*actual_state(:,index);
            counter2=0;
            counter1=1;
        else
             next_state_2=actual_state(:,index);
        end
        actual_state(:,index+1)=next_state_2+dt*(A2*next_state_2+B2*u+E2*v);
        actual_output(:,index+1)=C2*next_state_2;
    end
end

%% Fast Estimator simulation (Peaking phenomenon)
zeta_initial_mode_1=[3;3;3];
zeta_initial_mode_2=[3;3;3];
zeta_mode_1=[zeta_initial_mode_1,zeros(3,tf/dt+1)];
zeta_mode_2=[zeta_initial_mode_2,zeros(3,tf/dt+1)];
x_est_mode_1=zeros(3,tf/dt+1);
x_est_mode_2=zeros(3,tf/dt+1);
x_est_final=zeros(3,tf/dt+1);
counter1_est=0;
counter2_est=0;
mode_est=zeros(1,tf/dt+1);

for index=1:((tf/dt))
    t=(index-1)*dt;
    u=control_vec(index);
    
    zeta_mode_1(:,index+1)=zeta_mode_1(:,index)+dt*(H1_m*zeta_mode_1(:,index)+G1_m*u+L1_m*actual_output(:,index));
    x_est_mode_1(:,index)=zeta_mode_1(:,index)-J1_m*actual_output(:,index);
    r1=norm(C1*x_est_mode_1(:,index)-actual_output(:,index));
    
    zeta_mode_2(:,index+1)=zeta_mode_2(:,index)+dt*(H2_m*zeta_mode_2(:,index)+G2_m*u+L2_m*actual_output(:,index));
    x_est_mode_2(:,index)=zeta_mode_2(:,index)-J2_m*actual_output(:,index);
    r2=norm((C2*x_est_mode_2(:,index)-actual_output(:,index)));
    
    if r2<10^(-10)% mode 2 is active
        
        if counter2_est==1 % switches from mode 1 to mode 2
            next_state_1_to_2=phi*x_est_mode_1(:,index);
            x_est_final(:,index)=next_state_1_to_2;
            counter2_est=0;
            counter1_est=1;
        else  
            x_est_final(:,index)=x_est_mode_2(:,index);
        end
        mode_est(index)=2;     
    else
        if counter1_est==1
            next_state_2_to_1=phi*x_est_mode_2(:,index);
            x_est_final(:,index)=next_state_2_to_1;
            counter2_est=0;
            counter1_est=1;
        else
            x_est_final(:,index)=x_est_mode_1(:,index);
        end   
        mode_est(index)=1;     
    end   
end

subplot(4,1,1)
plot(0:dt:tf-dt,x_est_final(1,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(1,1:end-1),'b--','LineWidth',2)
legend('Estimated state 1','actual state 1')

subplot(4,1,2)
plot(0:dt:tf-dt,x_est_final(2,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(2,1:end-1),'b--','LineWidth',2)
legend('Estimated state 2','actual state 2')

subplot(4,1,3)
plot(0:dt:tf-dt,x_est_final(3,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(3,1:end-1),'b--','LineWidth',2)
legend('Estimated state 3','actual state 3')

subplot(4,1,4)
plot(0:dt:tf-dt,mode_est(1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,mode_vec(1:end-1),'b--','LineWidth',2)
xlabel('Time','FontSize',20)
legend('Estimated mode','actual mode')

%% Slow Estimator simulation 

zeta_initial_mode_1=[3;3;3];
zeta_initial_mode_2=[3;3;3];
zeta_mode_1=[zeta_initial_mode_1,zeros(3,tf/dt+1)];
zeta_mode_2=[zeta_initial_mode_2,zeros(3,tf/dt+1)];
x_est_mode_1=zeros(3,tf/dt+1);
x_est_mode_2=zeros(3,tf/dt+1);
x_est_final=zeros(3,tf/dt+1);
counter1_est=0;
counter2_est=0;
mode_est=zeros(1,tf/dt+1);

for index=1:((tf/dt))
    t=(index-1)*dt;
    u=control_vec(index);
    
    zeta_mode_1(:,index+1)=zeta_mode_1(:,index)+dt*(H1_c*zeta_mode_1(:,index)+G1_c*u+L1_c*actual_output(:,index));
    x_est_mode_1(:,index)=zeta_mode_1(:,index)-J1_c*actual_output(:,index);
    r1=norm(C1*x_est_mode_1(:,index)-actual_output(:,index));
    
    zeta_mode_2(:,index+1)=zeta_mode_2(:,index)+dt*(H2_c*zeta_mode_2(:,index)+G2_c*u+L2_c*actual_output(:,index));
    x_est_mode_2(:,index)=zeta_mode_2(:,index)-J2_c*actual_output(:,index);
    r2=norm((C2*x_est_mode_2(:,index)-actual_output(:,index)));
    
    if r2<10^(-4)  % mode 2 is active

        if counter2_est==1 % switches from mode 1 to mode 2
            next_state_1_to_2=phi*x_est_mode_1(:,index);
            x_est_final(:,index)=next_state_1_to_2;
            counter2_est=0;
            counter1_est=1;
        else  
            x_est_final(:,index)=x_est_mode_2(:,index);
        end
        mode_est(index)=2;     
    else
        if counter1_est==1
            next_state_2_to_1=phi*x_est_mode_2(:,index);
            x_est_final(:,index)=next_state_2_to_1;
            counter2_est=0;
            counter1_est=1;
        else
            x_est_final(:,index)=x_est_mode_1(:,index);
        end   
        mode_est(index)=1;     
    end   
end

subplot(4,1,1)
plot(0:dt:tf-dt,x_est_final(1,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(1,1:end-1),'b--','LineWidth',2)
legend('Estimated state 1','actual state 1')

subplot(4,1,2)
plot(0:dt:tf-dt,x_est_final(2,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(2,1:end-1),'b--','LineWidth',2)
legend('Estimated state 2','actual state 2')

subplot(4,1,3)
plot(0:dt:tf-dt,x_est_final(3,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(3,1:end-1),'b--','LineWidth',2)
legend('Estimated state 3','actual state 3')

subplot(4,1,4)
plot(0:dt:tf-dt,mode_est(1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,mode_vec(1:end-1),'b--','LineWidth',2)
xlabel('Time','FontSize',20)
legend('Estimated mode','actual mode')

%% Hybrid estimator (fast estimator to get mode info, slower to propagate state)
% No peaking phenomenon observed

zeta_initial_mode_1=[3;3;3];
zeta_initial_mode_2=[3;3;3];
zeta_mode_1=[zeta_initial_mode_1,zeros(3,tf/dt+1)];
zeta_mode_2=[zeta_initial_mode_2,zeros(3,tf/dt+1)];
zeta_mode_1_slow=[zeta_initial_mode_1,zeros(3,tf/dt+1)];
zeta_mode_2_slow=[zeta_initial_mode_1,zeros(3,tf/dt+1)];

x_est_mode_1=zeros(3,tf/dt+1);
x_est_mode_2=zeros(3,tf/dt+1);
x_est_mode_1_slow=zeros(3,tf/dt+1);
x_est_mode_2_slow=zeros(3,tf/dt+1);

x_est_final=zeros(3,tf/dt+1);
counter1_est=0;
counter2_est=0;
mode_est=zeros(1,tf/dt+1);

for index=1:((tf/dt))
    t=(index-1)*dt;
    u=control_vec(index);
    
    zeta_mode_1(:,index+1)=zeta_mode_1(:,index)+dt*(H1_m*zeta_mode_1(:,index)+G1_m*u+L1_m*actual_output(:,index));
    x_est_mode_1(:,index)=zeta_mode_1(:,index)-J1_m*actual_output(:,index);
    r1=norm(C1*x_est_mode_1(:,index)-actual_output(:,index));
    
    zeta_mode_2(:,index+1)=zeta_mode_2(:,index)+dt*(H2_m*zeta_mode_2(:,index)+G2_m*u+L2_m*actual_output(:,index));
    x_est_mode_2(:,index)=zeta_mode_2(:,index)-J2_m*actual_output(:,index);
    r2=norm((C2*x_est_mode_2(:,index)-actual_output(:,index)));
    
    zeta_mode_1_slow(:,index+1)=zeta_mode_1_slow(:,index)+dt*(H1_c*zeta_mode_1_slow(:,index)+G1_c*u+L1_c*actual_output(:,index));
    x_est_mode_1_slow(:,index)=zeta_mode_1_slow(:,index)-J1_c*actual_output(:,index);
    
    zeta_mode_2_slow(:,index+1)=zeta_mode_2_slow(:,index)+dt*(H2_c*zeta_mode_2_slow(:,index)+G2_c*u+L2_c*actual_output(:,index));
    x_est_mode_2_slow(:,index)=zeta_mode_2_slow(:,index)-J2_c*actual_output(:,index);
    
    if r2<10^(-10)  % mode 2 is active

        if counter2_est==1 % switches from mode 1 to mode 2
            next_state_1_to_2=phi*x_est_mode_1_slow(:,index);
            x_est_final(:,index)=next_state_1_to_2;
            counter2_est=0;
            counter1_est=1;
        else  
            x_est_final(:,index)=x_est_mode_2_slow(:,index);
        end
        mode_est(index)=2;     
    else
        if counter1_est==1
            next_state_2_to_1=phi*x_est_mode_2_slow(:,index);
            x_est_final(:,index)=next_state_2_to_1;
            counter2_est=0;
            counter1_est=1;
        else
            x_est_final(:,index)=x_est_mode_1_slow(:,index);
        end   
        mode_est(index)=1;     
    end   
end

subplot(4,1,1)
plot(0:dt:tf-dt,x_est_final(1,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(1,1:end-1),'b--','LineWidth',2)
legend('Estimated state 1','actual state 1')

subplot(4,1,2)
plot(0:dt:tf-dt,x_est_final(2,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(2,1:end-1),'b--','LineWidth',2)
legend('Estimated state 2','actual state 2')

subplot(4,1,3)
plot(0:dt:tf-dt,x_est_final(3,1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,actual_state(3,1:end-1),'b--','LineWidth',2)
legend('Estimated state 3','actual state 3')

subplot(4,1,4)
plot(0:dt:tf-dt,mode_est(1:end-1),'r','LineWidth',2)
hold on
grid on
plot(0:dt:tf-dt,mode_vec(1:end-1),'b--','LineWidth',2)
xlabel('Time','FontSize',20)
legend('Estimated mode','actual mode')
