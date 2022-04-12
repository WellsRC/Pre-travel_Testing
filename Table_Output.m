clear;
clc;
R0v=[5.08*(1+0.365) 5.08 2.79];
VOCv={'Omicron','Delta','Original'};
rng default;
pA=0.351; % Proportion of asymptomatic
par=fminbnd(@(x)((betainv(0.025,x,x*(1-pA)./pA)-0.307).^2+(betainv(0.975,x,x*(1-pA)./pA)-0.399).^2),[0],[400]);
pAv=betarnd(par,par*(1-pA)./pA,1000,1);
pAv=repmat(pAv,1,5);

load('RAgTest_Name_Final.mat','testName');
NumTests=length(testName);

MLE_Red=zeros(NumTests+1,5);
LB_95_Red=zeros(NumTests+1,5);
UB_95_Red=zeros(NumTests+1,5);

MLE_PT=zeros(NumTests+1,5);
LB_95_PT=zeros(NumTests+1,5);
UB_95_PT=zeros(NumTests+1,5);


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
    MLE_Red(1,:)=1-RRTPCR(w==0 | w==0.5 | w==1 | w==2 | w==3)./RNoTest(w==0 | w==0.5 | w==1 | w==2 | w==3);   
    
    MLE_PT_NT=100.*Probability_Onward(RNoTest(w==0 | w==0.5 | w==1 | w==2 | w==3),Risk);
    MLE_PT(1,:)=[ 100.*Probability_Onward(RRTPCR(w==0 | w==0.5 | w==1 | w==2 | w==3),Risk)];
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
    
    RA=repmat(RA(w==0 | w==0.5 | w==1 | w==2 | w==3),1000,1);
    RS=repmat(RS(w==0 | w==0.5 | w==1 | w==2 | w==3),1000,1);
    
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
    
    RS=reshape(RS(w~=7),1000,5);
    RA=reshape(RA(w~=7),1000,5);
    
    RRTPCRv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    R_temp=1-RRTPCRv./RNoTestv; 
    
           
    LB_95_Red(1,:)=prctile(R_temp,2.5,1);
    UB_95_Red(1,:)=prctile(R_temp,97.5,1);
    
    P_temp=100.*Probability_Onward(RNoTestv,Risk);
    LB_95_PT_NT=prctile(P_temp,2.5,1);
    UB_95_PT_NT=prctile(P_temp,97.5,1);
    
    P_temp=100.*Probability_Onward(RRTPCRv,Risk);
    
    LB_95_PT(1,:)=prctile(P_temp,2.5,1);
    UB_95_PT(1,:)=prctile(P_temp,97.5,1);
    
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
        MLE_Red(jj+1,:)=1-R_RATest(w==0 | w==0.5 | w==1 | w==2 | w==3)./RNoTest(w==0 | w==0.5 | w==1 | w==2 | w==3);   
        MLE_PT(jj+1,:)=[100.*Probability_Onward(R_RATest(w==0 | w==0.5 | w==1 | w==2 | w==3),Risk);];
        
        load(['Pre_Testing_' testName{jj} '_' VOCv{vv} '_Uncertainty.mat'],'PosttestISv','PosttestIAv','PretestISv','PretestIAv');
        w=repmat([0 0.5 1 2 3 7],1000,1);
        RS=(w./ts).*PosttestISv+((ts-w)./ts).*PretestISv;
        RA=(w./td).*PosttestIAv+((td-w)./td).*PretestIAv;

        RS(w==0)=PretestISv(w==0);
        RA(w==0)=PretestIAv(w==0);


        RS(w==ts)=PosttestISv(w==ts);
        RA(w==td)=PosttestIAv(w==td);

        
        RS=reshape(RS(w~=7),1000,5);
        RA=reshape(RA(w~=7),1000,5);

        R_RATestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);

        R_temp=1-R_RATestv./RNoTestv; 


        LB_95_Red(1+jj,:)=prctile(R_temp,2.5,1);
        UB_95_Red(1+jj,:)=prctile(R_temp,97.5,1);

        P_temp=100.*Probability_Onward(R_RATestv,Risk);

        LB_95_PT(jj+1,:)=prctile(P_temp,2.5,1);
        UB_95_PT(jj+1,:)=prctile(P_temp,97.5,1);
    end
    
    Time_Test_Red=cell(NumTests+1,5);
    Time_Test_PT=cell(NumTests+1,6);
    Name_Red=cell(NumTests+1,1);
    Name_PT=cell(NumTests+1,1);
    for jj=1:NumTests+1
           switch  jj
               case 1
                    Name_Red{jj}='RT-PCR';
                    Name_PT{jj}='RT-PCR';
               otherwise
                    Name_PT{jj}=AdjustedNames_Plotting(testName{jj-1});
                    Name_Red{jj}=AdjustedNames_Plotting(testName{jj-1});
           end
           Time_Test_PT{jj,1}=[num2str(MLE_PT_NT(1),'%4.2f') '% (' num2str(LB_95_PT_NT(1),'%4.2f') '%' char(8211) num2str(UB_95_PT_NT(1),'%4.2f') '%)'];               
           for ww=1:5
               Time_Test_PT{jj,7-ww}=[num2str(MLE_PT(jj,ww),'%4.2f') '% (' num2str(LB_95_PT(jj,ww),'%4.2f') '%' char(8211) num2str(UB_95_PT(jj,ww),'%4.2f') '%)'];               
               Time_Test_Red{jj,6-ww}=[num2str(100.*MLE_Red(jj,ww),'%4.2f') '% (' num2str(100.*LB_95_Red(jj,ww),'%4.2f') '%' char(8211) num2str(100.*UB_95_Red(jj,ww),'%4.2f') '%)'];               
           end
    end
    
    
    writetable(table(Name_Red,Time_Test_Red),[VOCv{vv} '_Reduction_Transmission.csv']);
    writetable(table(Name_PT,Time_Test_PT),[VOCv{vv} '_Probability_Transmission.csv']);  
end




