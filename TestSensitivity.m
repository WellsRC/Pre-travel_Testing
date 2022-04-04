function S = TestSensitivity(t,ts,testtype,beta,mm,PCR_Map)
%SensitivityvsViralLoad(V,asym) - Returns the sensitivity for a given viral
%load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ts- incubation period
% testtype - test type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S - Probability of that it is a true positive based on the viral load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if(strcmp(VOC,'Original'))
%     Gx =[6.127 3.277];
%     mm=Gx(1).*((Gx(2)-1)/Gx(2))^(1/Gx(2));
%     tt_max=46.2;
%     [S]=PCRSens(t,beta);
% elseif(strcmp(VOC,'Delta'))
%        
%     Gx =[4.230765931945100   0.959892250476516];
%     mm=(Gx(1)-1).*Gx(2);
%     tt_max=90;
%     load('PCR_Mapping.mat','VmI','VmS','t1','t2');
%     S(t<=mm)=pchip(VmI,PCRSens(t1,beta),Relative_Infection_PCR(t(t<=mm),VOC));
%     S(t>mm)=pchip(VmS,PCRSens(t2,beta),Relative_Infection_PCR(t(t>mm),VOC));
% elseif(strcmp(VOC,'Omicron'))
%     m=3.3;
%     v=3.5^2;
%     Gx=[log(m^2/sqrt(v+m^2)) sqrt(log(v/m^2+1))];
%     mm=exp(Gx(1)-Gx(2)^2);
%     tt_max=90;
%     
%     load('PCR_Mapping.mat','VmI','VmS','t1','t2');
%     S(t<=mm)=pchip(VmI,PCRSens(t1,beta),Relative_Infection_PCR(t(t<=mm),VOC));
%     S(t>mm)=pchip(VmS,PCRSens(t2,beta),Relative_Infection_PCR(t(t>mm),VOC));
% end


% PCR-Map
% load('MLE-Estimate-RTPCR.mat','beta')
% opts=optimset('TolX',10^(-16));
% ts=5.723;
% mm=fminbnd(@(x)-Infectivity_Profile(x,ts,inf),0,ts,opts);
% t1=linspace(0,mm,1001);
% t2=linspace(mm,50,1001);
% [SmI]=PCRSens(t1,beta);
% [SmS]=PCRSens(t2,beta);
% VmI=Relative_Infection_PCR(t1,5.723,mm);
% VmS=Relative_Infection_PCR(t2,5.723,mm);
 
 S(t<=mm)=pchip(PCR_Map.VmI,PCRSens(PCR_Map.t1,beta),Relative_Infection_PCR(t(t<=mm),ts,mm));
 S(t>mm)=pchip(PCR_Map.VmS,PCRSens(PCR_Map.t2,beta),Relative_Infection_PCR(t(t>mm),ts,mm));
 
if(~isempty(testtype))    
    V=Relative_Infection_PCR(t,ts,mm); % Use inf as we need to construct the mapping
    tt=[mm:0.1:90]; % made coarsered to improve the search tt=[ts:0.001:90]; 
    Vx=Relative_Infection_PCR(tt,ts,mm); % Use inf as we need to construct the mapping
    PPA=LR(tt-ts,testtype);
    S=pchip(Vx,PPA,V).*S;
end

end

