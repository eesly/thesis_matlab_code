clear all;
close all;
clc;

load('.\sea1_data\sea1_EI.mat');
load('.\sea1_data\sea1_ADCPpara.mat')
load('.\sea1_data\sea1_AD.mat')
addpath '.\function';
%% -----------------------------------------------------%
pN = length(cur_AD);
dN = length(cur_AD(1).AD);

fs = 98039;
c = ADCPpara.soundspeed;
bL = ADCPpara.blinddepth/10;
dL = 0;
bdN = floor((bL+dL)*2/c*fs);
layerN = floor(2*ADCPpara.celldepth/10*fs/1500);
blindL = ADCPpara.blinddepth/10 ;
layerH = ADCPpara.celldepth/10;
[sendS R I]= sig_s(fs);
ADCPpara.cellnumber = 47; 

%%
i = 10;
ADs = double(cur_AD(i).AD);
AD1 = ADs(1:8:end) + ADs(5:8:end)*1i;   
AD1 = AD1 - mean(AD1);

fftN = 4096;
fr = linspace(-fs/2,fs/2,fftN)./1000;
ADf = fftshift(abs(fft(AD1(1:fftN))));
ADf = ADf./max(ADf);

adr = real(AD1(1:layerN*ADCPpara.cellnumber+ bdN))./max(abs(real(AD1(1:layerN*ADCPpara.cellnumber+ bdN))));
adi = imag(AD1(1:layerN*ADCPpara.cellnumber+ bdN))./max(abs(imag(AD1(1:layerN*ADCPpara.cellnumber+ bdN))));
range = (1:length(adr))/fs/2*c;

figure(1)
subplot(211);
set(gca,'FontSize',14);

plot(range,adr,'LineStyle','-.','Color',[0 0 1]);hold on;
plot(range,adi,'LineStyle',':','Color',[0.87 0.49 0]);

xlim([min(range) max(range)]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid on;

legend('实部','虚部','FontSize',6)
xlabel('时间');
ylabel('幅度');

subplot(212);
set(gca,'FontSize',14);

plot(fr,ADf);
grid on;
set(gca,'YTickLabel',[]);

xlabel('频率(kHz)');
ylabel('幅度')
      
%% 
t_temp = [10*layerN:layerN*15] + bdN;
range = t_temp./fs/2*c;
anysis1 = adr(t_temp);
anysis2 = adi(t_temp);
anysis3 = abs(anysis1+1i*anysis2)/max(abs(anysis1+1i*anysis2));

figure(2)
fighist(anysis1,anysis2);     

%%
sig_r = adr + 1i*adi;
sig = R + 1i*I;
Rx_f = abs(conv(adr,R(end:-1:1),'same')); 

t = [1: ADCPpara.cellnumber * layerN] + bdN;
r = t./fs/2*c;  

% %对数正太分布估计
% m = mean(log(anysis3))
% s = std(log(anysis3))
% x = linspace(0.0001,1,100);
% width = 0.9999/103;
% f_lognormal = exp(-(log(x)-m).^2/2/s^2)./x./s./sqrt(2*pi)*width;  
% 
% %k分布参数估计（log-I 估计方法）
% m2 = sum(anysis3.^2)/length(anysis3);
% M1 = sum(log(anysis3.^2))/length(anysis3);
% M2 = log(sum(anysis3.^2)/length(anysis3));
% C = 0.57721566490;
% left = M1 - M2 + C;
% v_shape_array = 0.01:0.01:10;
% right = psi(v_shape_array) + log(v_shape_array);
% result = abs(left - right);
% v_shape = v_shape_array(find(result == min(result)));
% a_scale = sqrt(m2/v_shape)/2;
% f_k=2 * ((x/2/a_scale).^v_shape) .* besselk((v_shape-1),x/a_scale) ./ (a_scale*gamma(v_shape))*width; 
% plot(x,f_lognormal,'markersize',3,...
%     'MarkerFaceColor',[0 1 0],...
%     'Marker','square',...
%     'Color',[0 1 0]);
% hold on;
% plot(x,f_k,'MarkerFaceColor',[0.043 0.52 0.78],...
%     'Marker','o',...
%     'Color',[0.043 0.52 0.78],...
%     'markersize',3);
% hold off;

figure(3)
subplot(211)
set(gca,'FontSize',14);
plot(r,adr(t),'-blue');hold on;
plot(range,adr(t_temp),'-black');hold off;
 
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ylabel('幅度');
xlabel('时间')

subplot(212)
set(gca,'FontSize',14);
histfit2(anysis3);

grid on;  
xlim([0 1])

xlabel('归一化包络')
ylabel('PDF')
set(gca,'YTickLabel',[]);
legend('数据','rayleigh','lognormal','k')

pd3 = fitdist(anysis3,'rayleigh');
disp(['rayleith parameter:' num2str(pd3.Params)]);  
pd1 = fitdist(anysis1,'normal')

