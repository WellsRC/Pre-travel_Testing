close all;
clear;
clc;

pA=0.308; % Proportion of asymptomatic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_NoTest_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RNoTest=pA.*RA+(1-pA).*RS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RTPCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_RTPCR_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RRTPCR=pA.*RA+(1-pA).*RS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BinaxNOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_BinaxNow_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RBinax=pA.*RA+(1-pA).*RS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDVeritor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_BDVeritor_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RBDVeritor=pA.*RA+(1-pA).*RS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sofia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_Sofia_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RSofia=pA.*RA+(1-pA).*RS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lumarax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_Lumarax_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RLumarax=pA.*RA+(1-pA).*RS;

Lcolour=[hex2rgb('#9A9EAB');hex2rgb('#000B29'); hex2rgb('#807dba'); hex2rgb('#CF3721'); hex2rgb('#4CB5F5'); hex2rgb('#F5BE41');];
Lstyle={':','-','-.','-.','-.','-.'};
figure('units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
subplot('Position',[0.218510638297872,0.136779071547693,0.762872340425532,0.834686557505484]);
pdt=plot(w,RNoTest,w,RRTPCR,w,RBinax,w,RBDVeritor,w,RSofia,w,RLumarax,'LineWidth',2);
box off;
for ii=1:6
   pdt(ii).Color=Lcolour(ii,:); 
   pdt(ii).LineStyle=Lstyle{ii}; 
end
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:7],'XminorTick','on','YTick',[0:0.05:1],'Yminortick','on','Fontsize',22, 'xdir', 'reverse' );
ylim([0 0.9])
xlim([0 7]);
xlabel('Days pre-departure test conducted','Fontsize',28);
ylabel({'Post-travel','transmission'},'Fontsize',28);
% text(8.3,0.9,'A','Fontsize',34','FontWeight','bold');
legend({'No test','RT-PCR','BinaxNOW','BD Veritor','Sofia','LumiraDx'},'Fontsize',22,'Location','SouthWest')
legend boxoff;
print(gcf,['PreDeparture_Expected_Transmission_Hellewell.png'],'-dpng','-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Risk=1;

figure('units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
subplot('Position',[0.218510638297872,0.136779071547693,0.762872340425532,0.834686557505484]);
pdt=plot(w,Probability_Onward(RNoTest,Risk),w,Probability_Onward(RRTPCR,Risk),w,Probability_Onward(RBinax,Risk),w,Probability_Onward(RBDVeritor,Risk),w,Probability_Onward(RSofia,Risk),w,Probability_Onward(RLumarax,Risk),'LineWidth',2);
box off;
for ii=1:6
   pdt(ii).Color=Lcolour(ii,:); 
   pdt(ii).LineStyle=Lstyle{ii}; 
end
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:7],'XminorTick','on','YTick',[0.0:0.02:0.32],'Yminortick','on','Fontsize',22, 'xdir', 'reverse' );
ylim([0 0.32])
xlim([0 7]);
xlabel('Days pre-departure test conducted','Fontsize',28);
ylabel({'Probability of','post-travel transmission'},'Fontsize',28);
legend({'No test','RT-PCR','BinaxNOW','BD Veritor','Sofia','LumiraDx'},'Fontsize',22,'Location','SouthWest')
legend boxoff;
% text(8.86,0.32,'B','Fontsize',34','FontWeight','bold');
print(gcf,['PreDeparture_Probability_Transmission_Hellewell.png'],'-dpng','-r600');