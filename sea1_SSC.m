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
%% ƽ��EI--------------------------------------------------------%

%EL
[pN dN] = size(EI1);
ensemble_N = 5;
avEI = zeros(floor(pN/ensemble_N),dN);
for i = 1:floor(pN/ensemble_N)
    avEI(i,:) = mean(EI1((i-1)*ensemble_N+1:i*ensemble_N,:));
end
avEI(avEI == Inf | avEI == -Inf) = 114;
ADCPpara.blinddepth = ADCPpara.ranges(1);

%depth distance
maxN = floor(pN/ensemble_N)*ensemble_N;
depth_m = mean(reshape(depth(1:maxN),ensemble_N,maxN/ensemble_N));
depth_e = mean(reshape(edepth(1:maxN),ensemble_N,maxN/ensemble_N));

depth_range = ADCPpara.ranges./10;
distance = mean(reshape(distance(1:maxN),ensemble_N,maxN/ensemble_N));

%% ����Sv--------------------------------------------------------%

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
T = 27.5;       %�¶�
S = 35;         %�ζ�35ppt
Z = 0;          %���0km
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

%���
layerH = 2;
Vrec = ((ADCPpara.ranges./10+layerH/2).^3 - (ADCPpara.ranges./10-layerH/2).^3)*solid_angle/3;
Vrec1 = (ADCPpara.ranges./10).^2*solid_angle;

%Sv
Sv = zeros(floor(pN/ensemble_N),size(avEI,2));
for i = 1:floor(pN/ensemble_N)
    Sv(i,:) = avEI(i,:) + 2*rw*depth_range' + 20*log10(depth_range') - 216 - 10*log10(Vrec1');
%     Sv(i,:) = 0.45*(avEI(i,:) - NL) + 2*rw*depth_range' + 20*log10(depth_range') - 216;
end

%% Sv->SSC-------------------------------------------------------%
fc = 300e3;
band = 90e3;
frange = [fc-band/2 fc fc+band/2];

density = 2650;     %kg/m^3;
as = 60e-6
% SSC = 10.^(Sv/10)*100e5/2;
SSC = (10.^(Sv./10))*3*density*ADCPpara.soundspeed.^4 /(as).^3/frange(2).^4/1.082/4/pi^3;

%% ��Ч����ȥ��
earea = floor(-depth_e/2);
for i = 1:length(earea)
    SSC(i,earea(i):end) = 5;
end

%% figure
figure;
set(gca,'FontSize',14);

[X Y] = meshgrid(depth_range,distance);
contourf(X,Y,SSC,'LevelStep',0.1,'LineStyle','-');
% surf(X,Y,SSC,'LineStyle','none',...
%     'FaceColor','interp',...
%     'DisplayName','EI_effect');
hold on;
colorbar;
caxis([0 2]);
view([90 90]);

% plot(-depth_m,distance,'+black','linewidth',3);hold on;
plot(-depth_e,distance,'*blue','linewidth',3);hold on;

xlabel('���(m)');
ylabel('����(m)');
title('��ɳŨ��()')

xlim([depth_range(1) depth_range(end)]);
ylim([distance(1) distance(end)]);
% figure
% set(gca,'FontSize',18);
% step = floor(pN/avN/4);
% plot(ADCPpara.ranges(1:26)./10,SSC(step*4,1:26),'linewidth',2,...
%     'MarkerFaceColor',[0.16 0.38 0.27],'Marker','square','Color',[1 0 0]);hold on;
% plot(ADCPpara.ranges(1:26)./10,SSC(step*3,1:26),'linewidth',2,...
%     'MarkerFaceColor',[1 1 0],'Marker','o','Color',[0 0 1]);hold on;
% plot(ADCPpara.ranges(1:29)./10,SSC(step*2,1:29),'linewidth',2,...
%     'MarkerFaceColor',[0 0 1],'Marker','v','Color',[0 0 0]);hold on;
% plot(ADCPpara.ranges(1:36)./10,SSC(step*1,1:36),'linewidth',2,...
%     'MarkerFaceColor',[1 0 0],'MarkerSize',10,'Marker','*',...
%     'LineStyle','--','Color',[0.043 0.51 0.78]);hold off;
% xlim([min(ADCPpara.ranges./10) max(ADCPpara.ranges(1:36)./10)]);
% xlabel('Depth(m)');ylabel('SSC(kg m^-3)');
% legend(['t = ' num2str(step*4*avN*2) 's'],...
%        ['t = ' num2str(step*3*avN*2) 's'],...
%        ['t = ' num2str(step*2*avN*2) 's'],...
%        ['t = ' num2str(step*1*avN*2) 's'])
% grid on;
