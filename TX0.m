% Random entry in qiaratine with testing on exit
clear;

pobj=parpool(32); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w=sort([0:0.04:7  0.5]);
% 
SelfIsolate=1; % Self-isolation
% 
% % Allcoate memory for output
% PretestIS=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
% PosttestIS=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
% 
% PretestIA=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
% PosttestIA=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
% 
% Get Basline parameters
R0=1; % To scale easier after
R0S=R0; % Set R0 for symptomatic
R0A=R0; % Set R0 for asymptomatic

VOC={'Omicron','Delta','Original'};
tsvB=[3.1 4.4 5.723];

load('PCR_mapping.mat','VmI','VmS','t1','t2');
PCR_Map.VmI=VmI;
PCR_Map.VmS=VmS;
PCR_Map.t1=t1;
PCR_Map.t2=t2;

load('Peak_Infection.mat','mmv','tsv');
% 
% for vv=1:3
%     ts=tsvB(vv);
%     
%     mm=pchip(tsv,mmv,ts);
%     td=ts+20; % Asymptomatic increase 20 days from symptom onset
%     betaRTPCR=ParameterRTPCRTest(1);
% 
%     testtype=cell(1,1);
%     testtype{1}=[];
% 
% 
%     parfor jj=1:177
%     % Infected after the test was adminsted 
%         PosttestIS(jj)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,[],testtype,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR,mm,PCR_Map),0,w(jj),@(u)(w(jj)-u),ts));
%         PosttestIA(jj)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,[],testtype,R0S,R0A,1,ts,td,0,betaRTPCR,mm,PCR_Map),0,w(jj),@(u)(w(jj)-u),td));
% 
%          % Infected before the test was adminsted 
%         PretestIS(jj)=((1./(ts-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,0,testtype,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR,mm,PCR_Map),0,ts-w(jj),@(u)(u+w(jj)),ts));
%         PretestIA(jj)=((1./(td-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,0,testtype,R0S,R0A,1,ts,td,0,betaRTPCR,mm,PCR_Map),0,td-w(jj),@(u)(u+w(jj)),td));
% 
%     end
% 
%     save(['Pre_Testing_RTPCR_' VOC{vv} '.mat']);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% RT-PCR (Uncertainty)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% betaRTPCRv=zeros(1000,3);
% for jj=1:1000
%         [betaRTPCRv(jj,:)]=ParameterRTPCRTest(0);
% end
load(['RTPCR_Parameter_Sample_Uncertainty.mat'],'betaRTPCRv');

w=[0 0.5 1 2 3 7];
for vv=1:3
    ts=tsvB(vv);
    td=ts+20; % Asymptomatic increase 20 days from symptom onset
    
    mm=pchip(tsv,mmv,ts);
    testtype=cell(1,1);
    testtype{1}=[];

    % Allcoate memory for output
    PretestISv=zeros(1000,length(w)); % Assumes asymptomatic enter over infections period 
    PosttestISv=zeros(1000,length(w)); % Assumes asymptomatic enter over infections period 

    PretestIAv=zeros(1000,length(w)); % Assumes asymptomatic enter over infections period 
    PosttestIAv=zeros(1000,length(w)); % Assumes asymptomatic enter over infections period 
    
    
    % Allcoate memory for output
    PretestIS=zeros(1000,1); % Assumes asymptomatic enter over infections period 
    PosttestIS=zeros(1000,1); % Assumes asymptomatic enter over infections period 

    PretestIA=zeros(1000,1); % Assumes asymptomatic enter over infections period 
    PosttestIA=zeros(1000,1); % Assumes asymptomatic enter over infections period 
    
    for jj=1:length(w)
        parfor ss=1:1000
        betaRTPCR=betaRTPCRv(ss,:);
        % Infected after the test was adminsted 
            PosttestIS(ss)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,[],testtype,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR,mm,PCR_Map),0,w(jj),@(u)(w(jj)-u),ts));
            PosttestIA(ss)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,[],testtype,R0S,R0A,1,ts,td,0,betaRTPCR,mm,PCR_Map),0,w(jj),@(u)(w(jj)-u),td));

             % Infected before the test was adminsted 
            PretestIS(ss)=((1./(ts-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,0,testtype,R0S,R0A,0,ts,td,SelfIsolate,betaRTPCR,mm,PCR_Map),0,ts-w(jj),@(u)(u+w(jj)),ts));
            PretestIA(ss)=((1./(td-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,0,testtype,R0S,R0A,1,ts,td,0,betaRTPCR,mm,PCR_Map),0,td-w(jj),@(u)(u+w(jj)),td));

        end
        PretestISv(:,jj)=PretestIS;
        PosttestISv(:,jj)=PosttestIS;
        PretestIAv(:,jj)=PretestIA;
        PosttestIAv(:,jj)=PosttestIA;
    end

    save(['Pre_Testing_RTPCR_' VOC{vv} '_Uncertainty.mat'],'PretestISv','PosttestISv','PretestIAv','PosttestIAv','ts','td','w','R0S','R0A');
end

