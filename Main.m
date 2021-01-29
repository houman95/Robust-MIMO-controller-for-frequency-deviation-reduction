% clear all
close all
global Tg
Tg = 0.1;
Td = 5;
M = 0.15;
D = 0.008;
R = 3;
Tb = 0.1;
Rw = 14;
J = 62.993;
rho = 1.225;
Pwg = 350;
V = 692.82;
R1 = 0.00397;
R2 = 0.00443;
R3 = 0.0534;
N_meas = 3;
N_ctrl = 2;
clc
[A_P,B_P,C_P,D_P] = linmod('untitled1');
P = ss(A_P,B_P,C_P,D_P);
Iz = [1:2]'; 
Ie = [3:7]';
Iy = [8:10]'; 
Iv = [1:2]';
Iw = [3:6]'; 
Iu = [7:8]';
nz=numel(Iz);
nv=numel(Iv);
ne=numel(Ie);
nw=numel(Iw);
ny=numel(Iy); 
nu=numel(Iu); 

%% --------------------------design [e;y] <- [w;u]-------------------------
Pnomdesign = P([Ie;Iy],[Iw;Iu]); % select [e;y] <- [w;u]
Pnomdesign = minreal(Pnomdesign); % remove non-minimal states.

%% --------------------H_{\infty} controller K_nom ------------------------
nmeas=ny; nctrl=nu;
[K_nom,G_nom, gamma, info] = hinfsyn(Pnomdesign, nmeas, nctrl, 'METHOD','maxe');
G_nom_cl = lft(P, K_nom);
RS_blk = [1 1;1 1];
RP_blk = [RS_blk; numel(Iw) numel(Ie)];
omega = logspace(-6,4,250); % frequency vector
G_robust_f = frd(G_nom_cl,omega); % frequency response
%% --------------------------RS of the closed loop system with nominal controller --------------------------
mu_RS = mussv(G_robust_f(Iz,Iv),RS_blk); % robust stability
figure
[mu_RP,muinfo0] = mussv(G_robust_f,RP_blk);
A01=mu_RS.ResponseData(1,1,:);
A01=A01(:);
A02=mu_RP.ResponseData(1,1,:);
A02=A02(:);
semilogx(omega,A01,'LineWidth',1.5)

hold on
semilogx(omega,A02,'LineWidth',1.5)
a = svd(G_robust_f(Ie,Iw)).ResponseData(1,:);
semilogx(omega,a,'LineWidth',1.5,'color','k')
legend('Robust Stability','Robust Performance', 'Nominal Performance')
%% -------------------------------DK iteration-----------------------------
[D_L_0,~]=mussvunwrap(muinfo0);
D_L_0_perf=D_L_0(3,3);
D_L_0_1=D_L_0(1,1)/D_L_0_perf;
D_L_0_2=D_L_0(2,2)/D_L_0_perf;
D_L_0_1=D_L_0(1,1);
D_L_0_2=D_L_0(2,2);
DL0_1a=fitfrd(genphase(D_L_0_1),3);
DL0_2a=fitfrd(genphase(D_L_0_2),3);
sysDL_0 = [DL0_1a 0;0 DL0_2a];
%% ------------------------------- D-scaling transfer function fit -----------------------------

% figure
% h4=bodeplot(D_L_0_1,omega);
% hold on 
% h1=bodeplot(DL0_1a,omega);
% h3=bodeplot(D_L_0_2,omega);
% h2=bodeplot(DL0_2a,omega);
% setoptions(h1,'PhaseVisible','off');
% setoptions(h2,'PhaseVisible','off');
% setoptions(h3,'PhaseVisible','off');
% setoptions(h4,'PhaseVisible','off');
% legend('D1','D1 fit','D2','D2 fit')
%% 
sysDR_0 = inv(sysDL_0);
Pmu1design=[sysDL_0 zeros(nz,ne+nmeas);zeros(ne,nz) eye(ne) zeros(ne,nmeas);...
        zeros(nmeas,nz+ne),eye(nmeas)]*P*[sysDR_0 zeros(nv,nw+nctrl);...
        zeros(nw,nv),eye(nw),zeros(nw,nctrl);zeros(nctrl,nv+nw) eye(nctrl)];

%[K_mu_1,Gmu1,gamma1,info1] = hinfsyn(Pmu1design,nmeas,nctrl)

[K_mu_1,Gmu1,gamma1,info1] = hinfsyn(Pmu1design,nmeas,nctrl,...
'METHOD','lmi',... % LMI solution
'TOLGAM',0.1); % gamma tolerance;
G_cl_Robust = lft(P, K_mu_1);
Gmu1_f = frd(G_cl_Robust,omega); % frequency response
muRS1 = mussv(Gmu1_f(Iz,Iv),RS_blk); % robust stability
muNP1= svd(Gmu1_f(Ie,Iw)); % nominal performance
[muRP1,muinfo1] = mussv(Gmu1_f,RP_blk); % robust performance
A01=muRS1.ResponseData(1,1,:);
A01=A01(:);
A02=muRP1.ResponseData(1,1,:);
A02=A02(:);
figure
semilogx(omega,muNP1.ResponseData(1,:),'LineWidth',2)
hold on
semilogx(omega,A01,'LineWidth',2)
semilogx(omega,A02,'LineWidth',2)
xlabel('Frequency [rad/sec]')
grid on
title('Closed-loop robustness analysis with K1_{\mu}')
legend('Nominal Performance','Robust Stability','Robust Performance')  
%% --------------------2nd iteration---------------------
[D_L_2,~]=mussvunwrap(muinfo1);
D_L_2_perf=D_L_2(3,3);
D_L_2_1=D_L_2(1,1)/D_L_2_perf;
D_L_2_2=D_L_2(2,2)/D_L_2_perf;
DL2_1a=fitfrd(genphase(D_L_2_1),3);
DL2_2a=fitfrd(genphase(D_L_2_2),3);
sysDL_2 = [DL2_1a 0;0 DL2_2a];
sysDR_2 = inv(sysDL_2);
Pmu2design=[sysDL_2 zeros(nz,ne+nmeas);zeros(ne,nz) eye(ne) zeros(ne,nmeas);...
        zeros(nmeas,nz+ne),eye(nmeas)]*P*[sysDR_2 zeros(nv,nw+nctrl);...
        zeros(nw,nv),eye(nw),zeros(nw,nctrl);zeros(nctrl,nv+nw) eye(nctrl)];
 
[K_mu_2,~,gamma2,info2] = hinfsyn(Pmu2design,nmeas,nctrl); % gamma tolerance;
G_cl_Robust = lft(P, K_mu_2);
Gmu1_f = frd(G_cl_Robust,omega); % frequency response
muRS2 = mussv(Gmu1_f(Iz,Iv),RS_blk); % robust stability
muNP2= svd(Gmu1_f(Ie,Iw)); % nominal performance
[muRP2,muinfo2] = mussv(Gmu1_f,RP_blk); % robust performance
A01=muRS2.ResponseData(1,1,:);
A01=A01(:);
A02=muRP2.ResponseData(1,1,:);
A02=A02(:);
figure
semilogx(omega,muNP2.ResponseData(1,:),'LineWidth',2)
hold on

semilogx(omega,A01,'LineWidth',2)
semilogx(omega,A02,'LineWidth',2)
xlabel('Frequency [rad/sec]')
grid on
title('Closed-loop robustness analysis with K2_{\mu}')
legend('Nominal Performance','Robust Stability','Robust Performance')  
%% --------------------3rd iteration---------------------
[D_L_3,~]=mussvunwrap(muinfo2);
D_L_3_perf=D_L_3(3,3);
D_L_3_1=D_L_3(1,1)/D_L_3_perf;
D_L_3_2=D_L_3(2,2)/D_L_3_perf;
DL3_1a=fitfrd(genphase(D_L_3_1),3);
DL3_2a=fitfrd(genphase(D_L_3_2),3);
sysDL_3 = [DL3_1a 0;0 DL3_2a];
sysDR_3 = inv(sysDL_3);
Pmu3design=[sysDL_3 zeros(nz,ne+nmeas);zeros(ne,nz) eye(ne) zeros(ne,nmeas);...
        zeros(nmeas,nz+ne),eye(nmeas)]*P*[sysDR_3 zeros(nv,nw+nctrl);...
        zeros(nw,nv),eye(nw),zeros(nw,nctrl);zeros(nctrl,nv+nw) eye(nctrl)];
 
[K_mu_3,~,gamma3,info3] = hinfsyn(Pmu3design,nmeas,nctrl);
G_cl_Robust = lft(P, K_mu_3);
Gmu3_f = frd(G_cl_Robust,omega); % frequency response
muRS3 = mussv(Gmu3_f(Iz,Iv),RS_blk); % robust stability
muNP3= svd(Gmu3_f(Ie,Iw)); % nominal performance
[muRP3,muinfo2] = mussv(Gmu3_f,RP_blk); % robust performance
A01=muRS3.ResponseData(1,1,:);
A01=A01(:);
A02=muRP3.ResponseData(1,1,:);
A02=A02(:);
figure
semilogx(omega,muNP3.ResponseData(1,:),'LineWidth',2)
hold on
semilogx(omega,A01,'LineWidth',2)
semilogx(omega,A02,'LineWidth',2)
xlabel('Frequency [rad/sec]')
grid on
title('Closed-loop robustness analysis with K3_{\mu}')
legend('Nominal Performance','Robust Stability','Robust Performance') 
% ----------------------------- worst case ------------------------------
 mudata = frdata(muRP3);
 maxmu = max(mudata);
 
 maxidx = find(maxmu == max(maxmu));
 maxidx = maxidx(1);
  Delta0 = mussvunwrap(muinfo0);
 Delta0data = frdata(Delta0);
 Delta0_data_w = Delta0data(:,:,maxidx);
 
 %% --------Perturbed plant a transfer function to the perturbation---------
 close all
s=tf('s');
Delta0_wc = ss(zeros(nv,nz));
for i = 1:2
    delta_i = Delta0_data_w(i,i);
    gamma = abs(delta_i);
    if imag(delta_i) > 0 
        delta_i = -1*delta_i;
         gamma = -1*gamma;
    end
    x = real(delta_i)/abs(gamma);                      % fit a Pade with the same phase
    tau = 2*omega(maxidx)*(sqrt((1+x)/(1-x))); 
    Delta0_wc(i,i) = gamma*(-s+tau/2)/(s+tau/2); 
end
nDelta = norm(Delta0_data_w);                          % the size should be 1/mu.
Delta0_wc = Delta0_wc/nDelta;                          % scale the perturbation to be of size 1
%% --------Perturbed plant: close the loop around \delta-------------------
Ppert = lft(Delta0_wc,P);
Ppert_f = frd(Ppert,omega);  
hinf = lft(Ppert,K_nom);
mu = lft(Ppert,K_mu_2);
step(hinf(2,3))
hold on 
step(mu(2,3));
legend('H\infty Res.','\mu Res.')
title('Ref1 to error1')
figure
bodemag(hinf(2,3))
hold on
bodemag(mu(2,3))


%%  
type = 'pulse';
tau = 2;
Tf = 80000;
Ts = 0.1;
[u,t] = gensig(type,tau,Tf,Ts);
%u = 0.0001*u;
out2 = wgn(1,length(t),45);
%out = sin(1e-4*t);
out3 = wgn(1,length(t),20);

out = wgn(1,length(t),35);
figure
plot(t,[out2])
hold on
plot(t,[out])
plot(t,[out3])

figure
y = lsim(hinf(2,3),out,t) + lsim(hinf(2,2),out3,t)+ lsim(hinf(2,4),out2,t);
y2 = lsim(mu(2,3),out,t) + lsim(mu(2,2),out3,t) + lsim(mu(2,3),out2,t);

plot(t,60*[y],'LineWidth',3)
hold on 
% 
 plot(t,60*[y2],'LineWidth',3)
legend('Hinfinity controller','mu controller')


