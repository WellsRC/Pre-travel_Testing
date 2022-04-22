 
clear;

VOCv={'Omicron','Delta','Original'};
rng default;
pA=0.351; % Proportion of asymptomatic
par=fminbnd(@(x)((betainv(0.025,x,x*(1-pA)./pA)-0.307).^2+(betainv(0.975,x,x*(1-pA)./pA)-0.399).^2),[0],[400]);
pAv=betarnd(par,par*(1-pA)./pA,1000,1);


Expected_Transmission=zeros(3,1);
Expected_Transmission_V=zeros(3,1000);
for vv=1:3
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

    temp=pA.*RA.*(1./R0A)+(1-pA).*RS.*(1./R0S);
    Expected_Transmission(vv)=temp(w==0);
    
    
    load(['Pre_Testing_No_Testing_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);
    
    RA=repmat(RA(w==0),1000,1);
    RS=repmat(RS(w==0),1000,1);
    
    Expected_Transmission_V(vv,:)=pAv.*RA.*(1./R0A)+(1-pAv).*RS.*(1./R0S);
    
end

fprintf(['No Testing post-arrival transmission: Delta relative to Omicron: ' num2str(100.*(Expected_Transmission(2)./Expected_Transmission(1)-1),'%3.1f') '%% (' num2str(prctile(100.*(Expected_Transmission_V(2,:)./Expected_Transmission_V(1,:)-1),2.5),'%3.1f') '%%' char(8211) num2str(prctile(100.*(Expected_Transmission_V(2,:)./Expected_Transmission_V(1,:)-1),97.5),'%3.1f') '%%) \n']);      
fprintf(['No Testing post-arrival transmission: Original relative to Omicron: ' num2str(100.*(Expected_Transmission(3)./Expected_Transmission(1)-1),'%3.1f') '%% (' num2str(prctile(100.*(Expected_Transmission_V(3,:)./Expected_Transmission_V(1,:)-1),2.5),'%3.1f') '%%' char(8211) num2str(prctile(100.*(Expected_Transmission_V(3,:)./Expected_Transmission_V(1,:)-1),97.5),'%3.1f') '%%) \n']); 