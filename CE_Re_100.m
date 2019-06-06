%    clear all
 %   load ROMtestSV16KN20dt002_Re1
% % % 
%load snapshotData35Kdt002SV_Re100
   %load ROMtestSV16KN20dt002_Re10 
% % % %Snapshots MassROM StiffROM GlobalV PhiR MassMatrix

%load ROMtestSV35KN20dt002_Re100

%load Gsnap_SV35K_r6_d7_N16_166.mat
%load DNSProjMat16k_r2_Re_100
 %load Gsnap_SV16K_r4_d7_N20_Re_100
 %load Gsnap_SV16K_r4_d7_N20_Re_1

 
% load the Re=500 stuff
%load Gsnap_SV103K_r6_d7_N30_Re_500_219.mat
 
% d=7;
% r=3;

%problem = 'With_CE';
%problem = 'Without_CE';

% dfilter = 1;
% 
% data_driven_rom = 0;
% 
% my_factor = 0;


% endTimestep = 166;
% snapIndex =  1000;
% endTime =   0.332;


% skp=1;
% 
% delta = 1e-4;
% % % % 
%     endTimestep = 166;
% % % % % 
%      r = 16;
%      d = 16;
%     delta = 0.15;

MassROMdd = MassROM(1:d,1:d);
MassROMrr = MassROM(1:r,1:r);
MassROMrd = MassROM(1:r,1:d);
StiffROMrr = StiffROM(1:r,1:r);
StiffROMrd = StiffROM(1:r,1:d);



for ts=1:endTimestep
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We project snapshot on to the space X^d first to get a_d, the
    % coefficient of u_d
    vvv = Snapshots(:,ts+1000);
    
    RHS = zeros(d,1);
    for j=1:d
        RHS(j) = vvv' * (MassMatrix * PhiR(:,j) );
    end
    
    A = MassROMdd;
    RHS = RHS - A(:,1)*1;
    A(1,:)=0;
    A(:,1)=0;
    A(1,1)=1;
    RHS(1)=1;
    a_d = A  \ RHS;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute b_r for both differential and projection filters
    RHS1 = zeros(r,1);
    
    RHS1 = MassROMrd*a_d;
    
    system_mat_dfilter = MassROMrr+delta*delta*StiffROMrr;
    
    dfilter_a_r = system_mat_dfilter\RHS1;   %solving equation 11
    pfilter_a_r = MassROMrr\RHS1;            %solving equation 14
    
    
    RHS_pfilter = zeros(r,1);
    RHS_dfilter = zeros(r,1);
    
    RHS_pfilter = -StiffROMrr*pfilter_a_r;
    RHS_dfilter = -StiffROMrr*dfilter_a_r;
    
  
    pfilter_b_r = MassROMrr\RHS_pfilter;
    dfilter_b_r = MassROMrr\RHS_dfilter;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  Next we compute c_r for both differential and projection filters
    
    RHS2 = zeros(r,1);
    
    RHS2 = -StiffROMrd*a_d;
    
    pfilter_c_r = MassROMrr\RHS2;
    dfilter_c_r = system_mat_dfilter\RHS2;
    
    pfilter_CE(:,ts) = pfilter_b_r - pfilter_c_r;
    dfilter_CE(:,ts) = dfilter_b_r - dfilter_c_r;
    
    weak_pfilter_CE(:,ts) = (pfilter_CE(:,ts)'* MassROMrr)';
    weak_dfilter_CE(:,ts) = (dfilter_CE(:,ts)'* MassROMrr)';
end

% pfilter_error = 0.0;
% dfilter_error = 0.0;
% 
% for ts=1:endTimestep
%     pfilter_error = pfilter_error + pfilter_CE(:,ts)'*MassROMrr*pfilter_CE(:,ts);
%     dfilter_error = dfilter_error + dfilter_CE(:,ts)'*MassROMrr*dfilter_CE(:,ts);
% end
% 
% pfilter_CE_error = nu*sqrt(pfilter_error/endTimestep)
% dfilter_CE_error = nu*sqrt(dfilter_error/endTimestep)
    
% nu*sqrt(sum(sum(pfilter_CE.^2))/endTimestep)
% nu*sqrt(sum(sum(Gsnap.^2))/endTimestep)
clear RHS;
