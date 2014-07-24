function fighist(anysis1,anysis2)
subplot(211);
set(gca,'FontSize',14);
% span = linspace(min(anysis1),max(anysis1),100);
% [n,xout] = hist(anysis1,span);
% bar(xout,n./sum(n),'hist');
h = histfit1(anysis1);
% hold on;
% area(span,n./sum(n),'FaceColor',[0.168627455830574 0.505882382392883 0.337254911661148]);
xlim([-1,1]);
ylim([0 h*1.1])
ylabel('PDF');
set(gca,'XTickLabel',[]);
grid on;
% title('Probability density distribution of the instantaneous signal')

subplot(212);
set(gca,'FontSize',14);
% span = linspace(min(anysis2),max(anysis2),100);
% [n1,xout1] = hist(anysis2,span);
% bar(xout1,n1./sum(n1),'hist');
h = histfit1(anysis2);
% hold on;
% area(span,n1./sum(n1),'FaceColor',[0.168627455830574 0.505882382392883 0.337254911661148]);
% xlim([min(anysis2),max(anysis2)]);
xlim([-1,1]);
ylim([0 h*1.1])
ylabel('PDF');
xlabel('归一化瞬时值')
grid on;
% title('Probability density distribution of the signal envelop')