clear all;
close all;
clc;
load('EI1.mat');
load('ADCPpara.mat')
load('AD.mat')
load('vz.mat');
load('Vrms1.mat');
vr = 1;
%% 平均EI--------------------------------------------------------%
close all
[pN dN] = size(EI1);
avN= 4;
avEI = zeros(floor(pN/avN),dN);
avV = zeros(floor(pN/avN),dN);
for i = 1:floor(pN/avN)
    avEI(i,:) = mean(EI1((i-1)*avN+1:i*avN,:));
    avV(i,:) = mean(Vrms((i-1)*avN+1:i*avN,:));
end
avEI(avEI == Inf | avEI == -Inf) = 114;
ADCPpara.blinddepth = ADCPpara.ranges(1);
depth = ADCPpara.ranges./10;

figure;
set(gca,'FontSize',14);
[X Y] = meshgrid(depth,(1:floor(pN/avN))*avN*2);
% contourf(X,Y,avEI,'LineStyle','none',...
%     'FaceColor','interp',...
%     'DisplayName','EI_effect');
contourf(X,Y,avEI,'LevelStep',2)
colorbar;view([90 90]);
xlabel('Depth(m)');ylabel('Time(s)');zlabel('Backscatter Strength(dB ref 1uPa)')
xlim([depth(1) depth(end)]);ylim([1 floor(pN/avN)]*avN*2);
title('EL (dB ref 1\mupa)')

% 
% figure;
% step = floor(pN/avN/4);
% set(gca,'FontSize',18);
% plot(ADCPpara.ranges(1:26)./10,avEI(step*4,1:26),'linewidth',2,...
%     'MarkerFaceColor',[0.16 0.38 0.27],'Marker','square','Color',[1 0 0]);hold on;
% plot(ADCPpara.ranges(1:26)./10,avEI(step*3,1:26),'linewidth',2,...
%     'MarkerFaceColor',[1 1 0],'Marker','o','Color',[0 0 1]);hold on;
% plot(ADCPpara.ranges(1:29)./10,avEI(step*2,1:29),'linewidth',2,...
%     'MarkerFaceColor',[0 0 1],'Marker','v','Color',[0 0 0]);hold on;
% plot(ADCPpara.ranges(1:36)./10,avEI(step*1,1:36),'linewidth',2,...
%     'MarkerFaceColor',[1 0 0],'MarkerSize',10,'Marker','*',...
%     'LineStyle','--','Color',[0.043 0.51 0.78]);hold off;
% xlim([min(ADCPpara.ranges(1)./10) max(ADCPpara.ranges(36)./10)]);
% xlabel('Depth(m)');
% set(gca,'XTickLabel',[]);
% ylabel('EL(dB ref 1uPa)');ylabel('Backscatter Strength(dB ref 1uPa)');
% legend(['t = ' num2str(step*4*avN*2) 's'],...
%        ['t = ' num2str(step*3*avN*2) 's'],...
%        ['t = ' num2str(step*2*avN*2) 's'],...
%        ['t = ' num2str(step*1*avN*2) 's'])
% grid on;
% hold on;

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

theta = 4/180*pi;
solid_angle = 2*pi*(1-cos(theta/2));
layerH = 2;
Vrec = ((ADCPpara.ranges./10+layerH/2).^3 - (ADCPpara.ranges./10-layerH/2).^3)*solid_angle/3;

Vrec1 = layerH*solid_angle;

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
%% Sv->SSC-------------------------------------------------------%
fc = 300e3;
band = 90e3;
frange = [fc-band/2 fc fc+band/2];

density = 2650;     %kg/m^3;
if vr == 1
    as = 8e-6;
    SSC = 10.^(Sv/10)*100e3/2;
else
    as = 400e-6
    SSC = 10.^(Sv/10)*100e3/2;
%     SSC = (10.^(Sv/10))*3*density*ADCPpara.soundspeed.^4 /(as).^3/frange(2).^4/1.082/4/pi^3;
end

figure
set(gca,'FontSize',18);
step = floor(pN/avN/4);
plot(ADCPpara.ranges(1:26)./10,SSC(step*4,1:26),'linewidth',2,...
    'MarkerFaceColor',[0.16 0.38 0.27],'Marker','square','Color',[1 0 0]);hold on;
plot(ADCPpara.ranges(1:26)./10,SSC(step*3,1:26),'linewidth',2,...
    'MarkerFaceColor',[1 1 0],'Marker','o','Color',[0 0 1]);hold on;
plot(ADCPpara.ranges(1:29)./10,SSC(step*2,1:29),'linewidth',2,...
    'MarkerFaceColor',[0 0 1],'Marker','v','Color',[0 0 0]);hold on;
plot(ADCPpara.ranges(1:36)./10,SSC(step*1,1:36),'linewidth',2,...
    'MarkerFaceColor',[1 0 0],'MarkerSize',10,'Marker','*',...
    'LineStyle','--','Color',[0.043 0.51 0.78]);hold off;
xlim([min(ADCPpara.ranges./10) max(ADCPpara.ranges(1:36)./10)]);
xlabel('Depth(m)');ylabel('SSC(kg m^-3)');
legend(['t = ' num2str(step*4*avN*2) 's'],...
       ['t = ' num2str(step*3*avN*2) 's'],...
       ['t = ' num2str(step*2*avN*2) 's'],...
       ['t = ' num2str(step*1*avN*2) 's'])
grid on;

% 模型
speed = ADCPpara.soundspeed;
blindL = ADCPpara.blinddepth/10;
layerH = ADCPpara.celldepth/10;
layerN = ADCPpara.cellnumber;
vcu = ones(1,ADCPpara.cellnumber)*7.616328346675831;        %+表示靠近ADCP


pR = (ones(1,layerN))*as;
step = floor(pN/avN/4);
vC = vcu;
SLv = 10^1.3;
EL = zeros(4,3,layerN);
for j = 1:4
    for i = 1:3
        [ELv EL(j,i,:)]= micscat(layerH,layerN,pR,SSC(step*j,:),vz(step*j,:).*cos(20*pi/180)./100,frange(i),SLv,blindL); 
    end
end
% EL =EL+20;
t1 = 1:36;
t2 = 1:29;
t3 = 1:26;
t4 = 1:26;
%相关系数
s11 = reshape(EL(1,1,t1),1,length(t1));
s12 = reshape(EL(1,2,t1),1,length(t1));
s13 = reshape(EL(1,3,t1),1,length(t1));
R1 = CR(s12,avEI(step*1,t1));

s21 = reshape(EL(2,1,t2),1,length(t2));
s22 = reshape(EL(2,2,t2),1,length(t2));
s23 = reshape(EL(2,3,t2),1,length(t2));
R2 = CR(s22,avEI(step*2,t2));

s31 = reshape(EL(3,1,t3),1,length(t3));
s32 = reshape(EL(3,2,t3),1,length(t3));
s33 = reshape(EL(3,3,t3),1,length(t3));
R3 = CR(s32,avEI(step*3,t3));

s41 = reshape(EL(4,1,t3),1,length(t4));
s42 = reshape(EL(4,2,t3),1,length(t4));
s43 = reshape(EL(4,3,t3),1,length(t4));
R4 = CR(s42,avEI(step*4,t4));

figure1 = figure
boardW = 0.1;
gap = 0.01;
sizeL = (1 - boardW*2 - gap)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
subplot('Position',[boardW boardW+sizeL + gap  sizeL sizeL])
set(gca,'FontSize',14);
plot(avEI(step*1,t1),'linewidth',2,...
    'MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','*',...
    'LineStyle','--','Color',[1 0 0]);hold on;
plot(s11,'-black','linewidth',2);hold on;
plot(s13,'-black','linewidth',2);hold on;
plot(s12,'linewidth',1,...
    'MarkerFaceColor',[0.34 0.20 0.32],'MarkerSize',3,'Marker','o',...
    'LineStyle','-','Color',[0 0 0]);hold on;
highlow(s13', s11',ones(1,length(t1))'.*max(s13)*2,ones(1,length(t1))','black');
hold off
xlim([0 length(t1)+1]);grid on;
ylim([min(s11)-1 max(s13)+2]);
ylabel('Backscatter Strength');
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2
subplot('Position',[boardW+sizeL + gap boardW+sizeL + gap  sizeL sizeL]) 
set(gca,'FontSize',14);
plot(avEI(step*2,t2),'linewidth',2,...
    'MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','*',...
    'LineStyle','--','Color',[1 0 0]);hold on;
plot(s21,'-black','linewidth',2);hold on;
plot(s23,'-black','linewidth',2);hold on;
plot(s22,'linewidth',1,...
    'MarkerFaceColor',[0.34 0.20 0.32],'MarkerSize',3,'Marker','o',...
    'LineStyle','-','Color',[0 0 0]);hold on;
highlow(s23', s21',ones(1,length(t2))'.*max(s23)*2,ones(1,length(t2))','black');
hold off
xlim([0 length(t2)+1]);grid on;
ylim([min(s21)-1 max(s23)+2]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3
subplot('Position',[boardW boardW  sizeL sizeL]) 
set(gca,'FontSize',14);
plot(avEI(step*3,t3),'linewidth',2,...
    'MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','*',...
    'LineStyle','--','Color',[1 0 0]);hold on;
plot(s31,'-black','linewidth',2);hold on;
plot(s33,'-black','linewidth',2);hold on;
plot(s32,'linewidth',1,...
    'MarkerFaceColor',[0.34 0.20 0.32],'MarkerSize',3,'Marker','o',...
    'LineStyle','-','Color',[0 0 0]);hold on;
highlow(s33', s31',ones(1,length(t3))'.*max(s33)*2,ones(1,length(t3))','black');hold off
xlim([0 length(t3)]);grid on;
ylim([min(s31)-1 max(s33)+2]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4
subplot('Position',[boardW+sizeL + gap boardW  sizeL sizeL]) 
set(gca,'FontSize',14);
plot(avEI(step*4,t4),'linewidth',2,...
    'MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','*',...
    'LineStyle','--','Color',[1 0 0]);hold on;
plot(s41,'-black','linewidth',2);hold on;
plot(s43,'-black','linewidth',2);hold on;
plot(s42,'linewidth',1,...
    'MarkerFaceColor',[0.34 0.20 0.32],'MarkerSize',3,'Marker','o',...
    'LineStyle','-','Color',[0 0 0]);hold on;
highlow(s43', s41',ones(1,length(t4))'.*max(s43)*2,ones(1,length(t4))','black');hold off
xlim([0 length(t4)+1]);grid on; xlabel('Depth');
ylim([min(s41)-1 max(s43)+2]);
set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
annotation(figure1,'textbox',...
    [0.25 0.55 0.15 0.05],...
    'FontSize',12,...
    'FontName','Arial Rounded MT Bold',...
    'String',{['R =' num2str(R1)]},...
    'FitBoxToText','on');

% Create textbox
annotation(figure1,'textbox',...
    [0.66 0.55 0.15 0.05],...
    'FontSize',12,...
    'FontName','Arial Rounded MT Bold',...
    'String',{['R =' num2str(R2)]},...
    'FitBoxToText','on');

% Create textbox
annotation(figure1,'textbox',...
    [0.66 0.15 0.15 0.05],...
    'FontSize',12,...
    'FontName','Arial Rounded MT Bold',...
    'String',{['R =' num2str(R3)]},...
    'FitBoxToText','on');

% Create textbox
annotation(figure1,'textbox',...
    [0.25 0.15 0.15 0.05],...
    'FontSize',12,...
    'FontName','Arial Rounded MT Bold',...
    'String',{['R =' num2str(R4)]},...
    'FitBoxToText','on');
% 
% fig(avEI,SSC,t1,t2,t3,t4,step,ADCPpara);
% %%
E1 = avEI(step*1,t1);
E2 = avEI(step*2,t2);
E3 = avEI(step*3,t3);
E4 = avEI(step*4,t4);
SSC1 = SSC(step*1,t1);
SSC2 = SSC(step*2,t2);
SSC3 = SSC(step*3,t3);
SSC4 = SSC(step*4,t4);

