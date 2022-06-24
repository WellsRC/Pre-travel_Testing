clear;
close all;
clc;
R0v=[5.08*(1+0.365) 5.08 2.79];
VOCv={'Omicron','Delta','Original'};
rng default;
pA=0.351; % Proportion of asymptomatic
par=fminbnd(@(x)((betainv(0.025,x,x*(1-pA)./pA)-0.307).^2+(betainv(0.975,x,x*(1-pA)./pA)-0.399).^2),[0],[400]);
pAv=betarnd(par,par*(1-pA)./pA,1000,1);
pAv=repmat(pAv,1,6);

load('RAgTest_Name_Top.mat','testName');
testName=sort(testName);
NumTests=length(testName);

MLE_Red=zeros(NumTests+1,6);
LB_95_Red=zeros(NumTests+1,6);
UB_95_Red=zeros(NumTests+1,6);

Lcolour=[hex2rgb('#000B29'); hex2rgb('#807dba'); hex2rgb('#CF3721'); hex2rgb('#4CB5F5'); hex2rgb('#F5BE41');  hex2rgb('#A2C523');];

for vv=1:3
    R0VOC=R0v(vv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_No_Testing_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RNoTest=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RTPCR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_RTPCR_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RRTPCR=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);
    
    
    Risk=1;
    MLE_Red(1,:)=1-RRTPCR(w==0 | w==0.5 | w==1 | w==2 | w==3 | w==7)./RNoTest(w==0 | w==0.5 | w==1 | w==2 | w==3 | w==7);  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_No_Testing_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);
    
    RA=repmat(RA(w==0 | w==0.5 | w==1 | w==2 | w==3 | w==7),1000,1);
    RS=repmat(RS(w==0 | w==0.5 | w==1 | w==2 | w==3 | w==7),1000,1);
    
    RNoTestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RTPCR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_RTPCR_' VOCv{vv} '_Uncertainty.mat'],'PosttestISv','PosttestIAv','PretestISv','PretestIAv');
    w=repmat([0 0.5 1 2 3 7],1000,1);
    RS=(w./ts).*PosttestISv+((ts-w)./ts).*PretestISv;
    RA=(w./td).*PosttestIAv+((td-w)./td).*PretestIAv;

    RS(w==0)=PretestISv(w==0);
    RA(w==0)=PretestIAv(w==0);
    
    
    RS(w==ts)=PosttestISv(w==ts);
    RA(w==td)=PosttestIAv(w==td);
    
    RS=reshape(RS,1000,6);
    RA=reshape(RA,1000,6);
    
    RRTPCRv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    R_temp=1-RRTPCRv./RNoTestv; 
    
           
    LB_95_Red(1,:)=prctile(R_temp,2.5,1);
    UB_95_Red(1,:)=prctile(R_temp,97.5,1);
    
    for jj=1:NumTests
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RA TEST
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(['Pre_Testing_' testName{jj} '_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
        RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
        RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

        RS(w==0)=PretestIS(w==0);
        RA(w==0)=PretestIA(w==0);


        RS(w==ts)=PosttestIS(w==ts);
        RA(w==td)=PosttestIA(w==td);

        R_RATest=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);


        Risk=1;
        MLE_Red(jj+1,:)=1-R_RATest(w==0 | w==0.5 | w==1 | w==2 | w==3 | w==7)./RNoTest(w==0 | w==0.5 | w==1 | w==2 | w==3 | w==7);   
                
        load(['Pre_Testing_' testName{jj} '_' VOCv{vv} '_Uncertainty.mat'],'PosttestISv','PosttestIAv','PretestISv','PretestIAv');
        w=repmat([0 0.5 1 2 3 7],1000,1);
        RS=(w./ts).*PosttestISv+((ts-w)./ts).*PretestISv;
        RA=(w./td).*PosttestIAv+((td-w)./td).*PretestIAv;

        RS(w==0)=PretestISv(w==0);
        RA(w==0)=PretestIAv(w==0);


        RS(w==ts)=PosttestISv(w==ts);
        RA(w==td)=PosttestIAv(w==td);

        R_RATestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);

        R_temp=1-R_RATestv./RNoTestv; 


        LB_95_Red(1+jj,:)=prctile(R_temp,2.5,1);
        UB_95_Red(1+jj,:)=prctile(R_temp,97.5,1);
    end
     
  if(vv==1)  
        figure('units','normalized','outerposition',[0. 0. 0.6 1]);
        subplot('Position',[0.1285,0.7125,0.8614788732,0.28]);
        xx=[-0.465 -0.305 -0.15 0.005 0.16 0.315];
        dxx=0.145;
        for nt=1:6
            for ii=1:6
                patch(ii+[xx(nt) xx(nt)+dxx xx(nt)+dxx xx(nt)], 100.*[LB_95_Red(nt,ii) LB_95_Red(nt,ii) UB_95_Red(nt,ii) UB_95_Red(nt,ii)],Lcolour(nt,:),'LineStyle','none','facealpha',0.3); hold on;
                plot(ii+linspace(xx(nt),xx(nt)+dxx,101),100.*MLE_Red(nt,ii).*ones(101,1),'LineWidth',2,'color',Lcolour(nt,:));
            end
            if(nt==1)
                text(4.15,58,'RT-PCR','Color',Lcolour(nt,:),'Fontsize',20);
            else
                text(4.15,58-5.*(nt-1),AdjustedNames_Plotting(testName{nt-1}),'Color',Lcolour(nt,:),'Fontsize',20);                
            end
        end
        set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTickLabel','','XTick',[1:6],'YTick',[0:10:60]);
    grid on
        xlim([0.5 6.5]);
        ylim([0 60]);
        ytickformat('percentage');
        text(-0.355464759959142,58.0505415162455,'A','Fontsize',32,'FontWeight','bold');
  elseif(vv==2)
    subplot('Position',[0.1285,0.4087,0.8614788732,0.28]);
    
        for nt=1:6
            for ii=1:6
                patch(ii+[xx(nt) xx(nt)+dxx xx(nt)+dxx xx(nt)], 100.*[LB_95_Red(nt,ii) LB_95_Red(nt,ii) UB_95_Red(nt,ii) UB_95_Red(nt,ii)],Lcolour(nt,:),'LineStyle','none','facealpha',0.3); hold on;
                plot(ii+linspace(xx(nt),xx(nt)+dxx,101),100.*MLE_Red(nt,ii).*ones(101,1),'LineWidth',2,'color',Lcolour(nt,:));
            end
        end
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTickLabel','','XTick',[1:6],'YTick',[0:10:60]);
    grid on
    xlim([0.5 6.5]);
    ylim([0 60]);
    ylabel({'Reduction in post-arrival transmission relative to no test'},'Fontsize',20);
    ytickformat('percentage');
    text(-0.355464759959142,58.0505415162455,'B','Fontsize',32,'FontWeight','bold');
  else
    subplot('Position',[0.1285,0.4087-(0.7125-0.4087),0.8614788732,0.28]);
    
        for nt=1:6
            for ii=1:6
                patch(ii+[xx(nt) xx(nt)+dxx xx(nt)+dxx xx(nt)], 100.*[LB_95_Red(nt,ii) LB_95_Red(nt,ii) UB_95_Red(nt,ii) UB_95_Red(nt,ii)],Lcolour(nt,:),'LineStyle','none','facealpha',0.3); hold on;
                plot(ii+linspace(xx(nt),xx(nt)+dxx,101),100.*MLE_Red(nt,ii).*ones(101,1),'LineWidth',2,'color',Lcolour(nt,:));
            end
        end
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTickLabel',{'0 h','12 h','24 h','48 h','72 h','7 days'},'XTick',[1:6],'YTick',[0:10:60]);
    grid on
    xlabel('Time test taken prior to arrival','Fontsize',20);
    xlim([0.5 6.5]);
    ylim([0 60]);
    ytickformat('percentage');
    text(-0.355464759959142,58.0505415162455,'C','Fontsize',32,'FontWeight','bold');
  end

end
print(gcf,['Uncertainty reduction.png'],'-dpng','-r300');


