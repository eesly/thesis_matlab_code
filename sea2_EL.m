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
%% 平均EI与depth--------------------------------------------------------%

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

%depth
maxN = floor(pN/ensemble_N)*ensemble_N;
depth_m = mean(reshape(depth(1:maxN),ensemble_N,maxN/ensemble_N));
depth_e = mean(reshape(edepth(1:maxN),ensemble_N,maxN/ensemble_N));

%figure
depth_range = ADCPpara.ranges./10;
distance = mean(reshape(distance(1:maxN),ensemble_N,maxN/ensemble_N));

figure;
set(gca,'FontSize',14);

[X Y] = meshgrid(depth_range,distance);
contourf(X,Y,avEI,'LevelStep',1.5,'LineStyle','-');
hold on;
colorbar;
view([90 90]);

plot(-depth_m,distance,'+black','linewidth',3);hold on;
plot(-depth_e,distance,'*blue','linewidth',3);hold on;

xlabel('深度(m)');
ylabel('距离(m)');
title('回波强度(dB ref 1uPa)')

xlim([depth_range(1) depth_range(end)]);
ylim([distance(1) distance(end)]);


