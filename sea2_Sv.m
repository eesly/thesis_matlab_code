clear all;
close all;
clc;
load('.\sea2_data\sea2_EI.mat');
load('.\sea2_data\sea2_ADCPpara.mat')
load('.\sea2_data\sea2_AD.mat')
load('.\sea2_data\sea2_vz.mat');
load('.\sea2_data\sea2_depth');
load('.\sea2_data\sea2_edepth');
load('.\sea2_data\sea2_distance');
%% 平均EI--------------------------------------------------------%

%EL
[pN dN] = size(EI1);
ensemble_N = 5;
avEI = zeros(floor(pN/ensemble_N),dN);
for i = 1:floor(pN/ensemble_N)
    avEI(i,:) = mean(EI1((i-1)*ensemble_N+1:i*ensemble_N,:));
end
avEI(avEI == Inf | avEI == -Inf) = 115;
ADCPpara.ranges = ADCPpara.ranges(1:47);
ADCPpara.blinddepth = ADCPpara.ranges(1);

%depth distance
maxN = floor(pN/ensemble_N)*ensemble_N;
depth_m = mean(reshape(depth(1:maxN),ensemble_N,maxN/ensemble_N));
depth_e = mean(reshape(edepth(1:maxN),ensemble_N,maxN/ensemble_N));

depth_range = ADCPpara.ranges./10;
distance = mean(reshape(distance(1:maxN),ensemble_N,maxN/ensemble_N));

%% 计算Sv--------------------------------------------------------%

%NL
NSL = 31;
W = 63e3;
theta = 4/180*pi;
solid_angle = 2*pi*(1-cos(theta/2));
D = 4*pi/solid_angle;
DI = 10*log10(D);
NL = NSL + 10*log10(W) - DI;
NL

%rw
T = 20;          %温度
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
rw

%pulse
pulseW = 0.0012648;
pulseL = pulseW*ADCPpara.soundspeed;

%SL
SL = 216;

%体积
layerH = 2;
Vrec = ((ADCPpara.ranges./10+layerH/2).^3 - (ADCPpara.ranges./10-layerH/2).^3)*solid_angle/3;
Vrec1 = (ADCPpara.ranges./10).^2*solid_angle;

%Sv
Sv = zeros(floor(pN/ensemble_N),size(avEI,2));
for i = 1:floor(pN/ensemble_N)
    Sv(i,:) = avEI(i,:) + 2*rw*depth_range' + 20*log10(depth_range') - 216 - 10*log10(Vrec1');
%     Sv(i,:) = 0.45*(avEI(i,:) - NL) + 2*rw*depth_range' + 20*log10(depth_range') - 216;
end

%%无效区域去除
earea = floor(-depth_e/2);
for i = 1:length(earea)
    Sv(i,earea(i)+1:end) = -80;
end

%% figure
figure;
set(gca,'FontSize',14);

[X Y] = meshgrid(depth_range,distance);
contourf(X,Y,Sv,'LevelStep',2,'LineStyle','-');
% surf(X,Y,Sv,'LineStyle','none',...
%     'FaceColor','interp',...
%     'DisplayName','EI_effect');colorbar;view([90 90]);
hold on;
colorbar;
caxis([-80 -55]);
view([90 90]);

% plot(-depth_m,distance,'+black','linewidth',3);hold on;
plot(-depth_e,distance,'*blue','linewidth',3);hold on;

xlabel('深度(m)');
ylabel('距离(m)');
title('单位体积散射强度(dB ref 1uPa)')

xlim([depth_range(1) depth_range(end-20)]);
ylim([distance(1) distance(end)]);


