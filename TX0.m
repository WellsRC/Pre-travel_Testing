% Random entry in qiaratine with testing on exit
clear;

pobj=parpool(20); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=linspace(0,7,1001);
SelfIsolate=1; % Self-isolation
tL=[2.9]; % vecotor for the incbation periods to be integrated over

% Allcoate memory for output
PretestIS=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
PosttestIS=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 

PretestIA=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
PosttestIA=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 

% Get Basline parameters
[~,~,R0,ts] = BaselineParameters(tL); % Does not matter here what ts is fed in 


td=ts+20; % Asymptomatic increase 20 days from symptom onset

R0S=R0; % Set R0 for symptomatic
R0A=R0; % Set R0 for asymptomatic


load('MLE-Estimate-RTPCR-Hill.mat','beta')
betaRTPCR=beta;

testtype=cell(1,1);
testtype{1}=[];


parfor jj=1:1001 
    % Infected after the test was adminsted 
    PosttestIS(jj)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[],testtype,R0S,R0A,0,ts,tL,td,SelfIsolate,betaRTPCR),0,w(jj),0,inf));
    PosttestIA(jj)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,[],testtype,R0S,R0A,1,ts,tL,td,0,betaRTPCR),0,w(jj),0,inf));
    
     % Infected before the test was adminsted 
    PretestIS(jj)=((1./(ts-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,0,testtype,R0S,R0A,0,ts,tL,td,SelfIsolate,betaRTPCR),0,ts-w(jj),w(jj),inf));
    PretestIA(jj)=((1./(td-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,0,testtype,R0S,R0A,1,ts,tL,td,0,betaRTPCR),0,td-w(jj),w(jj),inf));
    
end

save('Testing_RTPCR_Hellewell.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR (OLD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=linspace(0,7,1001);
SelfIsolate=1; % Self-isolation
tL=[2.9]; % vecotor for the incbation periods to be integrated over

% Allcoate memory for output
PretestIS=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
PosttestIS=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 

PretestIA=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 
PosttestIA=zeros(1,length(w)); % Assumes asymptomatic enter over infections period 

% Get Basline parameters
[~,~,R0,ts] = BaselineParameters(tL); % Does not matter here what ts is fed in 


td=ts+20; % Asymptomatic increase 20 days from symptom onset

R0S=R0; % Set R0 for symptomatic
R0A=R0; % Set R0 for asymptomatic

testtype=cell(1,1);
testtype{1}=[];

parfor jj=1:1001 
    % Infected after the test was adminsted 
    PosttestIS(jj)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTestingOLD(t+u,u,[],testtype,R0S,R0A,0,ts,tL,td,SelfIsolate),0,w(jj),0,inf));
    PosttestIA(jj)=((1./w(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTestingOLD(t+u,u,[],testtype,R0S,R0A,1,ts,tL,td,0),0,w(jj),0,inf));
    
     % Infected before the test was adminsted 
    PretestIS(jj)=((1./(ts-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTestingOLD(t+u,u,0,testtype,R0S,R0A,0,ts,tL,td,SelfIsolate),0,ts-w(jj),w(jj),inf));
    PretestIA(jj)=((1./(td-w(jj))).*integral2(@(u,t)InfectiousnessfromInfectionTestingOLD(t+u,u,0,testtype,R0S,R0A,1,ts,tL,td,0),0,td-w(jj),w(jj),inf));
    
end

save('Testing_RTPCR_NatComm.mat');

clear;

