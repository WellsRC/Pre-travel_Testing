clear;
clc;
R0v=[5.08*(1+0.365) 5.08 2.79];
VOCv={'Omicron','Delta','Original'};
rng default;
pA=0.351; % Proportion of asymptomatic
par=fminbnd(@(x)((betainv(0.025,x,x*(1-pA)./pA)-0.307).^2+(betainv(0.975,x,x*(1-pA)./pA)-0.399).^2),[0],[400]);
pAv=betarnd(par,par*(1-pA)./pA,1000,1);

vv=1;


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
    RNoTest=RNoTest(w==0);
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
    
    RA=repmat(RA(w==0),1000,1);
    RS=repmat(RS(w==0),1000,1);

    RNoTestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    
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
    Red_RTPCR=1-RRTPCR(w==7)./RNoTest;   
    
    
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
    
    RS=RS(w==7);
    RA=RA(w==7);
    
    RRTPCRv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    R_temp=1-RRTPCRv./RNoTestv; 
    
 fprintf(['RT-PCR testing a week prior: ' num2str(100.*Red_RTPCR,'%3.2f') '%% (' num2str(prctile(100.*R_temp,2.5),'%3.2f') '%%' char(8211) num2str(prctile(100.*R_temp,97.5),'%3.2f') '%%) \n']);      