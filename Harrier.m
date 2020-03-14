%% EEE-588/ Design of Multi-variable feedback control system
%Final Project of EEE588
%Topic: AV-8A Harrier Aircraft Longitudinal Dynamicss
%Modelling of plant which includes Actuator and Engine dynamics
%of the motors.
%********************************************************************************
%Model of the plant
%Ap, Bp, Cp, Dp 
%After including actuator there are total 6 states and 2 inputs and 2   
%outputs.
%
clear all;
clc;
close all;
Ap= [0 1 0 0 0 0;
    -1.837 -1.893 1.837 -0.0004 0.0062 -0.1243;
    0.5295 0.0085 -0.5295 0.0006 0.0002 0.0017;
    -34.5 0 2.3 -0.0621 0.4209 -0.0452;
    0 0 0 0 -1.966 0;
    0 0 0 0 0 -12];
Bp= [0 0;
    0 0;
    0 0;
    0 0;
    1.966 0;
    0 12;];
Cp= [0 0 0 1 0 0;
    0 0 57.2958 0 0 0];
Dp= 0*ones(2,2);

P = ss(Ap,Bp,Cp,Dp); %Nominal Plant model
%***************************************************************************
%% Plant Dimensions
%
[ns,nc] = size(Bp);                   % Number of States, Number of Controls
no      = nc;                         % Number of Outputs
%% Modal Analysis - Eigenvalue-EigenVector function
%Natural Modes:

[evec,eval] = eig(Ap); %Poles of nominal plant

%***************************************************************************
%% Transmission Zeros
%
plantzeros = tzero(ss(Ap,Bp,Cp,Dp))       % transmission zeros
%zdir = null([z*eye(4)-Ap -Ap; Cp Dp])   % transmission zero directions
%***************************************************************************
%% SYSTEM TRANSFER FUNCTIONS: From u_i to y_j
%
Plant_zpk = zpk(ss(Ap,Bp,Cp,Dp)) % Zeros, Poles, and Gains fron u_i to x_j

%% Controllability and Observability
%
% Controllability 
%
CM  = [Bp Ap*Bp (Ap^2)*Bp  (Ap^3)*Bp];%   (ap^4)*bp ]  % Controllability Matrix
rcm = rank(CM)                                      % Rank of Controllability Matrix
%
% system is controllable; that is, rank of controllability matrix is 6 
%
%***************************************************************************
%
% Observability
%
OM  = [Cp;Cp*Ap; Cp*(Ap^2);Cp*(Ap^3)]; %Cp*(ap^4) ]          % Observability Matrix
rom = rank(OM)              % Rank of Observability Matrix
%
% system is observable; that is, rank of observability matrix is 6
%
%***************************************************************************
%% Modal analysis
% Visualization of Phugoid Mode (Alpha approximately Constant)
%
% this lightly damped mode takes very long to decay
% the aircraft climbs and drops exchanging kinetic and potential energy
%
t1 = [0:0.05:500];
x = lsim(ss(Ap,Bp, eye(6,6), 0*ones(6,2)), [0*t1; 0*t1] ,t1, real(evec(:,3)) );
figure,
subplot(1,2,1)
plot(t1,x)
grid
title('Visualization of F8 Phugoid Mode (x_o = Re x_{ph})')
xlabel('time (seconds)')
ylabel('x(t)')
%
x = lsim(ss(Ap,Bp, eye(6,6), 0*ones(6,2)), [0*t1; 0*t1] , t1, imag(evec(:,3)) );
subplot(1,2,2)
plot(t1,x)
grid
title('Visualization of F8 Phugoid Mode (x_o = Imag x_{ph})')
xlabel('time (seconds)')
ylabel('x(t)')
%
% Visualization of Short Period Mode (Speed approximately Constant)
%
% this lightly damped mode decays very quickly
%
t1 = [0:0.05:5];
x = lsim(ss(Ap,Bp, eye(6,6), 0*ones(6,2)), [0*t1; 0*t1] , t1, real(evec(:,1)) );
figure,
subplot(1,2,1)
plot(t1,x)
grid
title('Visualization of F8 Short Period Mode (x_o = Re x_{sp})')
xlabel('time (seconds)')
ylabel('x(t)')
%
x = lsim(ss(Ap,Bp, eye(6,6), 0*ones(6,2)), [0*t1; 0*t1] , t1, imag(evec(:,1)) );
subplot(1,2,2)
plot(t1,x)
grid
title('Visualization of F8 Short Period Mode (x_o = Im x_{sp})')
xlabel('time (seconds)')
ylabel('x(t)')
%% Frequency Response: Singular Values of Plant
P = ss(Ap,Bp,Cp,Dp);
winit  = -3;
wfin   =  3;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv_p     = sigma(P,w);
sv_p     = 20*log10(sv_p);
figure; semilogx(w, sv_p)
%clear sv
title('SVD of nominal plant')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%***************************************************************************
%
% SVD ANALYSIS at DC
%
dc =  Cp*inv(-Ap)*Bp + Dp;
[udc,sdc,vdc] = svd(dc);

G_o = udc*sdc*vdc;


figure,
subplot(2,2,1)
bar(sdc(1,1)*udc(1,:))
title('\sigma_11*u1')
grid
xlabel('Output [vel FPA]')
ylabel('\sigma_11*u1')
subplot(2,2,2)
bar(vdc(:,1))
title('v1')
grid
xlabel('Controls [stick_input throttle ]')
ylabel('v1')
%***************************************************************************
subplot(2,2,3)
bar(abs(sdc(2,2)*udc(2,:)))
title('\sigma_22*u2')
grid
xlabel('Output [vel FPA]')
ylabel('\sigma_22*u2')
subplot(2,2,4)
bar(abs(vdc(:,2)))
title('v2')
grid
xlabel('Controls [stick_input throttle ]')
ylabel('v2')


%***************************************************************************
%
% SVD Analysis at a selected frequency
%
s1 = j*0.0975;                                       % At Phugoid Mode
g1 = Cp*inv(s1*eye(6)-Ap)*Bp + Dp;
[u1, sing1, v1 ] = svd(g1);
v1mag = abs(v1);
v1phase = angle(v1)*180/pi;
u1mag = abs(u1);
u1phase = angle(u1)*180/pi;

%***************************************************************************
%
% SVD Analysis at a selected frequency
%
s2 = j*1.1672;                                        % At Short Period Mode
g2 = Cp*inv(s2*eye(6)-Ap)*Bp + Dp;
[u2, sing2, v2 ] = svd(g2);
v2mag = abs(v2);
v2phase = angle(v2)*180/pi;
u2mag = abs(u2);
u2phase = angle(u2)*180/pi;


%***************************************************************************
%% Augment Plant with Integrators at Plant Input and Plot Singular Values
Ad=[Ap Bp 
   0*ones(nc,ns) 0*ones(nc,nc)];
Bd=[0*ones(ns,nc)
    eye(nc)];
Cd=[Cp 0*ones(nc,nc)];
Dd=0*ones(nc,nc);

%***************************************************************************
% Plant analysis after augmentation
Pd=ss(Ad,Bd,Cd,Dd);
pp=eig(Pd);% Eigen values of new plant
ptz=tzero(Pd);%transmission zeros of new plant

figure;
sv_aug=sigma(Pd,w); 
tsv_aug=20*log10(sv_aug);
semilogx(w,tsv_aug);
title('SVD of the plant after Dynamic Augmentation');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');



%***************************************************************************
%% Bilinear transformation
p1=-0.0944; p2=-1e20;
[A_blin,B_blin,C_blin,D_blin]=bilin(Ad,Bd,Cd,Dd,1,'Sft_jw',[p2 p1]);
%***************************************************************************
% Plant analysis after Bilinear
Pd_lin=ss(A_blin,B_blin,C_blin,D_blin);
[evec_alin,eval_alin] = eig(A_blin);%eigenvalues after bilinear trans
auglin_plantzeros = tzero(ss(A_blin,B_blin,C_blin,D_blin));%transmission zeros after bilinear trans
figure;
sv_lin=sigma(Pd_lin,w); 
tsv_lin=20*log10(sv_lin);
semilogx(w,tsv_lin);
title('SVD of the plant after Bilinear transformation');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');
%***************************************************************************
%% LQR Design

Q=C_blin'*C_blin;
rho = 1e-12;                                          
R = rho*eye(nc);
% 
% [G, poles, K] = lqr(A_blin,B_blin,Q,R);
% Q = eye(8)% State Weighting Matrix
% % Q1(1,1) = 2000;
% % Q1(2,2) = 100;
%  Q(3,3) = 10;
%  Q(4,4) = 1;
% rho = 1e-9;                                           % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
R = rho*eye(nc)                                      % Control Weigthing Matrix
[G, K, poles] = lqr(A_blin,B_blin,Q,R); 
Glqr=ss(A_blin,B_blin,G,0);
lqrp=eig(A_blin,B_blin*G);
lqrtz=tzero(Glqr);

[Ak,Bk,Ck,Dk]=bilin(Glqr.A,Glqr.B,Glqr.C,Glqr.D,-1,'Sft_jw',[p2 p1]);
K1=ss(Ak,Bk,Ck,Dk);

So1=inv(eye(size(K1))+K1);
To1=eye(size(K1))-So1;

figure;
sv=sigma(So1,w);
tsv=20*log10(sv);
semilogx(w,tsv);
title('Sensitivity (S) using LQR');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');

figure;
sv=sigma(To1,w);
tsv=20*log10(sv);
semilogx(w,tsv);
title('Complementary Sensitivity (T) using LQR');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');

% t = [0:0.02:5];
% [y, t, x] = step(To1,t);
% figure;
% plot(t,y(:,:,2))
% grid
% title('Pitch & Gamma Response To r = [1  0]  Command')
% ylabel('Theta & Gamma (deg)')
% xlabel('Time (seconds)')
%***************************************************************************
%% Kalman filter
Pint=eye(nc);
mu=0.1;
Mint=mu*eye(nc);
[Kkf, H, sig]=kalman(Pd_lin,Pint,Mint);
Gkf=ss(A_blin-H*C_blin,H,-C_blin,0);

% Plant analysis after Kalman filter
kfp=eig(A_blin-H*C_blin);
kftz=tzero(Gkf);
%***************************************************************************
%% LQG Design
%
% Am=A_blin-B_blin*G-H*C_blin;
% Bm=H;
% Cm=G;
% Dm=0;
% Glqg=ss(Am,Bm,Cm,Dm);

%Form final compensator
Am = [0*ones(2,2)  G;
       0*ones(8,2) A_blin-B_blin*G-H*C_blin];
Bm = [0*ones(2,2);
        H];
Cm = [eye(2,2) 0*ones(2,8)]

Dm = [0*ones(2,2)];
Glqg =ss(Am,Bm,Cm,Dm);
%***************************************************************************
% Analysis after LQG augmentation

lqgp=eig(A_blin-B_blin*G-H*C_blin);
lqgtz=tzero(Glqg);

[Ak,Bk,Ck,Dk]=bilin(Glqg.A,Glqg.B,Glqg.C,Glqg.D,-1,'Sft_jw',[p2 p1]);
K2=ss(Ak,Bk,Ck,Dk);
[Lo2,Li2,So2,Si2,To2,Ti2,KS2,SP2]=f_CLTFM(P,K2);

%**************************************************************************
%
% OPEN LOOP FREQUENCY RESPONSE 
figure; 
sigma(Lo2,w);
title('LQG Loop Singular Values: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sigma(Li2,w);
title('LQG Loop Singular Values: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%**************************************************************************
%
% CLOSED LOOP FREQUENCY RESPONSE 

figure; 
sigma(So2,w);
title('LQG: Sensitivity Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Si2,w);
title('LQG: Sensitivity Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(To2,w);
title('LQG:Complementary Sensitivity Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Ti2,w);
title('LQG:Complementary Sensitivity Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sv=sigma(KS2,w);
tsv=20*log10(sv);
semilogx(w,tsv);
title('K-Sensitivity (KS) using LQG');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');
%***************************************************************************
% Time response
t = [0:0.02:10];
[y, t, x] = step(To2,t);
figure;
subplot(1,2,1)
plot(t,y(:,:,1))
grid
title('Velocity & FPA To r = [1  0]  Command')
ylabel('velocity & FPA (deg)')
xlabel('Time (seconds)')

t = [0:0.02:20];
[y, t, x] = step(To2,t);
subplot(1,2,2)
plot(t,y(:,:,2))
grid
title('Velocity & FPA To r = [0  1]  Command')
ylabel('velocity & FPA (deg)')
xlabel('Time (seconds)')
%***************************************************************************
%% LQG-LTR Design
%
% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(Cp*inv(-Ap)*Bp + Dp);
lh = -inv(Ap)*Bp*ll;
l = [lh;ll]; 

figure; 
sv = sigma(ss(A_blin, l, C_blin, D_blin),w);
sv = 20*log10(sv);
semilogx(w, sv)
title('Filter Open Loop (G_FOL)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
%***************************************************************************
Pint=eye(nc);
mu=0.1;
Mint=mu*eye(nc);
[Kkf1, H1, sig]=kalman(ss(A_blin, [B_blin l], C_blin, [D_blin 0*ones(nc,nc)]),Pint,Mint);

% Recover Target Loop By Solving Cheap LQR Problem
%
Q1 = C_blin'*C_blin;
% Q1 = eye(8)% State Weighting Matrix
% % Q1(1,1) = 2000;
% % Q1(2,2) = 100;
%  Q1(3,3) = 10;
%  Q1(4,4) = 100;
rho = 1e-9;                                           % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
R1 = rho*eye(nc)                                      % Control Weigthing Matrix
[G1, poles, rrr] = lqr(A_blin,B_blin,Q1,R1);                   % Compute Control Gain Matrix G

%***************************************************************************
figure;
sv=sigma(ss(A_blin,H1,C_blin,D_blin),w);
tsv=20*log10(sv);
semilogx(w,tsv);
title('Target Loop (G_KF)');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');

%***************************************************************************
tolpoles = eig(A_blin)                           % Target Open Loop Poles
targzeros = tzero(A_blin,H1,C_blin,0*ones(nc,nc))      % Target Open Loop Zeros
tclpoles = eig(A_blin-H1*C_blin)                       % Target Closed Loop Poles

%***************************************************************************
figure
sv = sigma(ss(A_blin-H1*C_blin, H1, -C_blin, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Target Sensitivity (S_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

sv = sigma(ss(A_blin-H1*C_blin, H1, C_blin, 0*eye(nc)),w);
sv = 20*log10(sv);
figure
semilogx(w, sv, w, 20*log10(10./w))
title('Target Complementary (T_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%***************************************************************************
% Analysis LQG-LTR
% Am1=A_blin-B_blin*G-H1*C_blin;
% Bm1=H1;
% Cm1=G;
% Dm1=0;
% Gltr=ss(Am1,Bm1,Cm1,Dm1);
%Form final compensator
Am1 = [0*ones(2,2)  G1;
       0*ones(8,2) A_blin-B_blin*G1-H1*C_blin];
Bm1 = [0*ones(2,2);
        H1];
Cm1 = [eye(2,2) 0*ones(2,8)]

Dm1 = [0*ones(2,2)];
Gltr=ss(Am1,Bm1,Cm1,Dm1);
ltrp=eig(Am1);
ltrtz=tzero(Gltr);

[ak,bk,ck,dk]=bilin(Gltr.A,Gltr.B,Gltr.C,Gltr.D,-1,'Sft_jw',[p2 p1]);
K3=ss(ak,bk,ck,dk);
[Lo3,Li3,So3,Si3,To3,Ti3,KS3,SP3]=f_CLTFM(P,K3);
%**************************************************************************
%
% OPEN LOOP FREQUENCY RESPONSE 
figure; 
sigma(Lo3,w);
title('LQG-LTR Loop Singular Values: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sigma(Li3,w);
title('LQG-LTR Loop Singular Values: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%**************************************************************************
%
% CLOSED LOOP FREQUENCY RESPONSE 

figure; 
sigma(So3,w);
title('LQG-LTR: Sensitivity Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Si3,w);
title('LQG-LTR: Sensitivity Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(To3,w);
title('LQG-LTR:Complementary Sensitivity Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Ti3,w);
title('LQG-LTR:Complementary Sensitivity Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sv=sigma(KS3,w);
tsv=20*log10(sv);
semilogx(w,tsv);
title('K-Sensitivity (KS) using LQG-LTR');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');
%***************************************************************************
% Time response
t = [0:0.02:10];
[y, t, x] = step(To3,t);
figure;
subplot(1,2,1)
plot(t,y(:,:,1))
grid
title('Velocity & FPA To r = [1  0]  Command')
ylabel('velocity & FPA (deg)')
xlabel('Time (seconds)')

t = [0:0.02:20];
[y, t, x] = step(To3,t);
subplot(1,2,2)
plot(t,y(:,:,2))
grid
title('Velocity & FPA To r = [0  1]  Command')
ylabel('velocity & FPA (deg)')
xlabel('Time (seconds)')
%***************************************************************************
%% H-infinity Robust Constrol Design
% Weights
% Standard first order weights are chosen
% For any non-standard/higher-order weights, form the transfer functions
% appropiately
%***************************************************************************
% Visualizing the Weighting matrix for Robust control
%
M11=7; w11=0.01; Eps11=0.003; M12=7; w12=0.01; Eps12=0.003;
W1 = [tf([1/M11 w11], [1 w11*Eps11]) 0; 0 tf([1/M12 w12], [1 w12*Eps12])];
figure;
sigma(inv(W1),w)
title('Inverse of Weighting matrix W1');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');
M21=0.1; w21=10; Eps21=0.01; M22=0.1; w22=10; Eps22=0.01;
W2 = [tf([1 w21/M21], [Eps21 w21]) 0; 0 tf([1 w22/M22], [Eps22 w22])] ;
figure;
sigma(inv(W2),w)
title('Inverse of Weighting matrix W2');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');
M31=2; w31=20; Eps31=0.001; M32=2; w32=20; Eps32=0.001; 
W3 = [tf([1 w31/M31], [Eps31 w31]) 0; 0 tf([1 w32/M32], [Eps32 w32])] ;
figure;
sigma(inv(W3),w)
title('Inverse of Weighting matrix W3');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');
%*****************************************************************************
% M11=1.2; w11=0.01; Eps11=0.001; M12=1.2; w12=0.05; Eps12=0.001;
% W1 = [tf([1/M11 w11], [1 w11*Eps11]) 0 ; 0 tf([1/M12 w12] , [1 w12*Eps12])];
% figure;
% sigma(inv(W1),w)
% title('Inverse of Weighting matrix W1');
% grid on;
% xlabel('Frequency (rad/sec)');
% ylabel('Singular Values (dB)');
% M21=100; w21=1000; Eps21=0.01; M22=100; w22=1000; Eps22=0.01;
% W2 = [tf([1 w21/M21], [Eps21 w21]) 0 ; 0 tf([1 w22/M22], [Eps22 w22])] ;
% figure;
% sigma(inv(W2),w)
% title('Inverse of Weighting matrix W2');
% grid on;
% xlabel('Frequency (rad/sec)');
% ylabel('Singular Values (dB)');
% M31=1.6; w31=3; Eps31=0.001; M32=1.6; w32=18; Eps32=0.001; 
% W3 = [tf([1 w31/M31], [Eps31 w31]) 0; 0 tf([1 w32/M32], [Eps32 w32])] ;
% figure;
% sigma(inv(W3),w)
% title('Inverse of Weighting matrix W3');
% grid on;
% xlabel('Frequency (rad/sec)');
% ylabel('Singular Values (dB)');

%***************************************************************************
p11=-0.0944; p22=-1e20;
[A_blin1,B_blin1,C_blin1,D_blin1]=bilin(Ad,Bd,Cd,Dd,1,'Sft_jw',[p22 p11]);
Pd_lin1 = ss(A_blin1,B_blin1,C_blin1,D_blin1)
%***************************************************************************
% % Generalized plant
% % Standard weight augmentation using Matlab's augw command
% % For any non-standard augmentation of weights, form the Generalized plant
% % manually using state-space methods
GenP=augw(Pd_lin1,W1,W2,W3);
Kd=hinfsyn(GenP);

[akd,bkd,ckd,dkd] = ssdata(Kd);
[ak4,bk4,ck4,dk4]=bilin(akd,bkd,ckd,dkd,-1,'Sft_jw',[p22 p11]);
Kh=ss(ak4,bk4,ck4,dk4);

[Lo4,Li4,So4,Si4,To4,Ti4,KS4,PS4] = f_CLTFM(Pd,Kh);
%***************************************************************************
%Analysing system after H-infinity design
%
% OPEN LOOP FREQUENCY RESPONSE 
figure; 
sigma(Lo4,w);
title('H-infinity-Open Loop Singular Values: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sigma(Li4,w);
title('H-infinity-Open Loop Singular Values: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%**************************************************************************
%
% CLOSED LOOP FREQUENCY RESPONSE 

figure; 
sigma(So4,w);
title('H-infinity-Sensitivity: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Si4,w);
title('H-infinity-Sensitivity: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(To4,w);
title('H-infinity-Complementary Sensitivity: Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Ti4,w);
title('H-infinity-Complementary Sensitivity: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sv=sigma(KS4,w);
tsv=20*log10(sv);
semilogx(w,tsv);
title('K-Sensitivity (KS) using H-infinity');
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Singular Values (dB)');

%%
%**************************************************************************
%
% CLOSED LOOP TIME RESPONSE 
t = [0:0.02:150];
[y, t, x] = step(To4,t);

% PITCH COMMAND
% Pitch: r = [ 1  0 ] Pitch Command

figure;
subplot(1,2,1)
%step(To4)
plot(t,y(:,:,1))
grid
title('Pitch & Gamma Response To r = [1  0]  Command')
ylabel('Theta & Gamma (deg)')
xlabel('Time (seconds)')
subplot(1,2,2)
plot(t,y(:,:,2))
grid
title('Pitch & Gamma Response To r = [0  1]  Command')
ylabel('Theta & Gamma (deg)')
xlabel('Time (seconds)')
%
%end
