function V = Relative_Infection_PCR(t,ts,mm)
%ViralShedding_Symptomatic(t,tL) reutrns the level of infectiousness for an symptomatic individual at
%time t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t - time post-infection
%tL- duration of the latent period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V - the amount of viral shedding at time t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(strcmp(VOC,'Original'))
%     Gx =[6.127 3.277];
%     mm=Gx(1).*((Gx(2)-1)/Gx(2))^(1/Gx(2));
% elseif(strcmp(VOC,'Delta'))
%     Gx =[4.230765931945100   0.959892250476516];
%     mm=(Gx(1)-1).*Gx(2);
% elseif(strcmp(VOC,'Omicron'))
%     m=3.3;
%     v=3.5^2;
%     Gx=[log(m^2/sqrt(v+m^2)) sqrt(log(v/m^2+1))];
%     mm=exp(Gx(1)-Gx(2)^2);
% end
V=Infectivity_Profile(t,ts,inf)./Infectivity_Profile(mm,ts,inf);

end

