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


depth_range = ADCPpara.ranges./10;
distance = mean(reshape(distance(1:maxN),ensemble_N,maxN/ensemble_N));

Sv = zeros(floor(pN/ensemble_N),size(avEI,2));
for i = 1:floor(pN/ensemble_N)
    Sv(i,:) = 0.45*(avEI(i,:) - NL) + 2*rw*depth_range' + 20*log10(depth_range') - 216;
end

figure;
set(gca,'FontSize',14);

[X Y] = meshgrid(depth_range,distance);
surf(X,Y,Sv,'LineStyle','none',...
    'FaceColor','interp',...
    'DisplayName','EI_effect');colorbar;view([90 90]);

xlabel('Depth(m)');
ylabel('Time(s)');
zlabel('S_v(dB ref 1uPa)');
title('S_v(dB ref 1uPa)');

xlim([depth_range(1) depth_range(end)]);
ylim([distance(1) distance(end)]);


