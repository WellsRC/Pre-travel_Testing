function V = ViralShedding_Symptomatic(t,ts,td)

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
%     V=wblpdf(t,Gx(1),Gx(2))./wblcdf(td,Gx(1),Gx(2));
% elseif(strcmp(VOC,'Delta'))
%     Gx =[4.230765931945100   0.959892250476516];
%     V=gampdf(t,Gx(1),Gx(2))./gamcdf(td,Gx(1),Gx(2));
% elseif(strcmp(VOC,'Omicron'))
%     m=3.3;
%     v=3.5^2;
%     Gx=[log(m^2/sqrt(v+m^2)) sqrt(log(v/m^2+1))];
%     V=lognpdf(t,Gx(1),Gx(2))./logncdf(td,Gx(1),Gx(2));
% end
% V(t>td)=0;

V=Infectivity_Profile(t,td,ts);



VC=integral(@(x)Infectivity_Profile(x,ts,td),0,td);

V=V./VC;
end

