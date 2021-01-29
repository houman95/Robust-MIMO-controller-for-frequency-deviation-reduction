close all
P = bodeoptions
P.FreqUnits = 'Hz';
Ww2 = tf([0.0005],[1 0.00001]);

Wbe = tf([20 100],[1 0.001]);
Wse = tf([50 0.001],[0.5 0.1]);

bodemag(Ww2,Wbe,Wse,P)
ylim([-100,120])
% Wcg = tf([1 0.0001],[5 1]);
% Wcb = tf([0.2 0.1],[100 0.1]);
% 
% bodemag(Wcg,Wcb,P)
legend
set(findall(gcf,'type','line'),'linewidth',2)
