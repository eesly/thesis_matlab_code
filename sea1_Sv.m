clear all;
close all;
clc;
load('.\sea1_data\sea1_EI.mat');
load('.\sea1_data\sea1_ADCPpara.mat')
load('.\sea1_data\sea1_AD.mat')
load('.\sea1_data\sea1_vz.mat');
load('.\sea1_data\sea1_depth');
load('.\sea1_data\sea1_edepth');
load('.\sea1_data\sea1_distance');
%% 平均EI--------------------------------------------------------%
%EL
[pN dN] = size(EI1);
ensemble_N = 5;
avEI = zeros(floor(pN/ensemble_N),dN);
for i = 1:floor(pN/ensemble_N)
    avEI(i,:) = mean(EI1((i-1)*ensemble_N+1:i*ensemble_N,:));
end
avEI(avEI == Inf | avEI == -Inf) = 114;
ADCPpara.blinddepth = ADCPpara.ranges(1);

%depth
maxN = floor(pN/ensemble_N)*ensemble_N;
depth_m = mean(reshape(depth(1:maxN),ensemble_N,maxN/ensemble_N));
depth_e = mean(reshape(edepth(1:maxN),ensemble_N,maxN/ensemble_N));

%% 计算Sv--------------------------------------------------------%
NSL = 31;
W = 63e3;
theta = 4/180*pi;
solid_angle = 2*pi*(1-cos(theta/2));
D = 4*pi/solid_angle;
DI = 10*log10(D);
NL = NSL + 10*log10(W) - DI;

T = 10;         %温度10摄氏度
S = 35;         %盐度35ppt
Z = 0;          %深度0km
PH = 8;
f = 300e3;
f1 = 780*exp(T/29);
f2 = 42000*exp(T/18);
A = 0.083*(S/35)*exp(T/31 - Z/91 + 1.8*(PH-8));
B = 22*(S/35)*exp(T/14 - Z/6);
C = 4.9*10^(-10)*exp(-T/26 - Z/25);
rw = (A ./ (f1^2 + f.^2) + B ./ (f2^2 + f.^2) + C).*f.^2./1000;

pulseW = 0.0012648;
pulseL = pulseW*ADCPpara.soundspeed;

SL = 216;
P = 10^(SL - DI - 170.8)/10;

% Con = -85;
% Kc = 0.9;
Con = 20;

% theta = 4/180*pi;
% solid_angle = 2*pi*(1-cos(theta/2));
% layerH = 2;
% Vrec = ((ADCPpara.ranges./10+layerH/2).^3 - (ADCPpara.ranges./10-layerH/2).^3)*solid_angle/3;
% 
% Vrec1 = layerH*solid_angle;

Sv = zeros(size(avEI));
for i = 1:floor(pN/avN)
      Sv(i,:) = avEI(i,:) + 2*rw*depth' + 20*log10(depth') - 216 - 10*log10(Vrec1') - Con;
%       Sv(i,:) = 0.9*(avEI(i,:) - NL) + 2*rw*depth' + 20*log10(depth') - 216  - 10*log10(Vrec1');
end

% figure;
% set(gca,'FontSize',14);
% [X Y] = meshgrid(depth,(1:floor(pN/avN))*avN/2);
% surf(X,Y,Sv,'LineStyle','none',...
%     'FaceColor','interp',...
%     'DisplayName','EI_effect');colorbar;view([90 90]);
% xlabel('Depth(m)');ylabel('Time(s)');zlabel('S_v(dB ref 1uPa)');title('S_v(dB ref 1uPa)');
% xlim([depth(1) depth(end)]);ylim([1 floor(pN/avN)]*avN/2);

% figure
% set(gca,'FontSize',14);
% step = floor(pN/avN/4);
% plot(ADCPpara.ranges./10,Sv(step*4,:),'linewidth',1,...
%     'MarkerFaceColor',[0.16 0.38 0.27],'Marker','square','Color',[1 0 0]);hold on;
% plot(ADCPpara.ranges./10,Sv(step*3,:),'linewidth',1,...
%     'MarkerFaceColor',[1 1 0],'Marker','o','Color',[0 0 1]);hold on;
% plot(ADCPpara.ranges./10,Sv(step*2,:),'linewidth',1,...
%     'MarkerFaceColor',[0 0 1],'Marker','v','Color',[0 0 0]);hold on;
% plot(ADCPpara.ranges./10,Sv(step*1,:),'linewidth',1,...
%     'MarkerFaceColor',[1 0 0],'MarkerSize',10,'Marker','*',...
%     'LineStyle','--','Color',[0.043 0.51 0.78]);hold off;
% xlim([min(ADCPpara.ranges./10) max(ADCPpara.ranges./10)]);
% xlabel('Depth(m)');ylabel('S_v(dB ref 1uPa)');
% legend(['t = ' num2str(step*4*avN*2) 's'],...
%        ['t = ' num2str(step*3*avN*2) 's'],...
%        ['t = ' num2str(step*2*avN*2) 's'],...
%        ['t = ' num2str(step*1*avN*2) 's'])
% grid on;





