close all;
clear;
clc;
VOCv={'Omicron','Delta','Original'};
tsvB=[3.1 4.4 5.723];

load('PCR_mapping.mat','VmI','VmS','t1','t2');
PCR_Map.VmI=VmI;
PCR_Map.VmS=VmS;
PCR_Map.t1=t1;
PCR_Map.t2=t2;
figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Omicron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.061491596638655,0.592256332320174,0.424600840336134,0.390000000000001]);
vv=1;
ts1=tsvB(vv);

load('Peak_Infection.mat','mmv','tsv');
% opts=optimset('TolX',10^(-16));
%   tsv=[0.1:0.1:28];
%   mmv=zeros(length(tsv),1);
%   for ii=1:280
%   mmv(ii)=fminbnd(@(x)-Infectivity_Profile(x,tsv(ii),inf),0,tsv(ii),opts);
%   end
%  save('Peak_Infection.mat');
mm=pchip(tsv,mmv,ts1);

load(['MLE-Estimate-RTPCR.mat'],'beta');
t=linspace(0,50,1001);
S1 = TestSensitivity(t,tsvB(vv),[],beta,mm,PCR_Map);
load(['RTPCR_Parameter_Sample_Uncertainty'],'betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),length(t));
for ii=1:length(betaRTPCRv(:,1))
    SU(ii,:)=  TestSensitivity(t,tsvB(vv),[],betaRTPCRv(ii,:),mm,PCR_Map);
end
patch([t flip(t)],[prctile(SU,2.5) flip(prctile(SU,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on
plot(t,S1,'k','LineWidth',2);
plot(tsvB(vv).*ones(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);

box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days post-infection','Fontsize',22);

text(-6.30407911001236,0.995,'A','Fontsize',32,'FontWeight','bold');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Non-Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.061491596638655,0.114,0.424600840336134,0.390000000000001]);

vv=3;
ts3=tsvB(vv);

mm=pchip(tsv,mmv,ts3);
load(['MLE-Estimate-RTPCR.mat'],'TPtID','TDate','TResult','PtID','par','beta');
TI=par(1:end-2);
[~,b]=ismember(TPtID,PtID);
dt=round(TDate'-TI(b(b>0)));
t=linspace(0,50,1001);

S3 = TestSensitivity(t,tsvB(vv),[],beta,mm,PCR_Map);

load(['RTPCR_Parameter_Sample_Uncertainty'],'betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),length(t));
for ii=1:length(betaRTPCRv(:,1))
   SU(ii,:)=  TestSensitivity(t,tsvB(vv),[],betaRTPCRv(ii,:),mm,PCR_Map);
end
patch([t flip(t)],[prctile(SU,2.5) flip(prctile(SU,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on
LL1=plot(t,S3,'k','LineWidth',2);
plot(tsvB(vv).*ones(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);

L=[1:40];
S=zeros(1,40);

for ii=1:40
   S(ii)=sum(TResult(dt<=(ii+0.5)&dt>(ii-0.5)))./length(dt(dt<=(ii+0.5)&dt>(ii-0.5))); 
end
LL2=scatter(L,S,40,'k','filled');
legend([LL2,LL1],{'Avg Test Result:Binned','Estimated RT-PCR'},'Fontsize',20);

box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days post-infection','Fontsize',22);

text(-6.30407911001236,0.995,'C','Fontsize',32,'FontWeight','bold');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
subplot('Position',[0.561491596638655,0.592256332320174,0.424600840336134,0.390000000000001]);
vv=2;
ts2=tsvB(vv);

mm=pchip(tsv,mmv,ts2);
load(['MLE-Estimate-RTPCR.mat'],'beta');

S2 = TestSensitivity(t,tsvB(vv),[],beta,mm,PCR_Map);

load(['RTPCR_Parameter_Sample_Uncertainty'],'betaRTPCRv');
SU=zeros(length(betaRTPCRv(:,1)),length(t));

for ii=1:length(betaRTPCRv(:,1))
   SU(ii,:)=  TestSensitivity(t,tsvB(vv),[],betaRTPCRv(ii,:),mm,PCR_Map);
end
patch([t flip(t)],[prctile(SU,2.5) flip(prctile(SU,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on
plot(t,S2,'k','LineWidth',2);
plot(tsvB(vv).*ones(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);

box off;
set(gca,'LineWidth',2,'tickdir','out','XTick',[0:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days post-infection','Fontsize',22);
text(-6.30407911001236,0.995,'B','Fontsize',32,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% All together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.561491596638655,0.114,0.424600840336134,0.390000000000001]);

plot(t-ts1,S1,'k',t-ts2,S2,'b',t-ts3,S3,'r','LineWidth',2); hold on
plot(zeros(101,1),linspace(0,1,101),':','color',[0.75 0.75 0.75],'LineWidth',2);
legend({'3.1 day incubation period','4.4 day incubation period','5.72 day incubation period'},'Fontsize',19);
legend boxoff;
box off;
xlim([-5.5 40]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[-5:5:50],'Xminortick','on','YTick',[0:0.1:1],'YminorTick','on','Fontsize',20);
ylabel('Diagnostic sensitivity','Fontsize',22);
xlabel('Days since symptom onset','Fontsize',22);

text(-11.2533992,0.995,'D','Fontsize',32,'FontWeight','bold');
print(gcf,['RT-PCR_Curves.png'],'-dpng','-r600');