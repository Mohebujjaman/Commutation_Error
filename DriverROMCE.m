clear all;

 load ROMtestSV35KN20dt002_Re100  
 load DNSProjMat35k_r13_Re_100
 %load Gsnap_SV35K_r5_d9_N20_Re_100_166
 %load Gsnap_SV16K_r4_d7_N20_Re_1
 
 %Snapshots MassROM StiffROM TriLinROM NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode


d=14;
r=13;

problem = 'With_CE';
%problem = 'Without_CE';
dfilter = 0;

data_driven_rom = 1;

%Do you want to add commutation error in Gsnap?
my_factor = 1;


endTimestep = 166;
snapIndex =  1000;
endTime =   0.332;


skp=2;

delta = 0.00;
createGsnap_9;

%if(data_driven_rom)
 
CE_Re_100;

  % unconstrained optimization (DDF-ROM)

  My_tol_DDROM = 5e-3;%  best 1.17e-3<=tol<9.9e-3 gives 4.6792e-04
%  5.5124e-05 5.2493e-05
% very good #3 3e-15 2.98e-15<=tol<=3.2e-15

%For projection filter ROM error
%r=3 cetol=3.2e-16
%r=4 cetol=3e-16
%r=5 cetol= 7e-15

%Differential filter ROM error

%r=5 celtol=4.4e-15 dtol all \delta=1e-4
tic
createABtilde_noconstraints;
toc
%end

ROMDriverPlotCE;

error_CE_DDC = my_error;

CEdataTableDDC = dataTableDDC;

%Without CE 

my_factor = 0;
data_driven_rom = 1;

%My_tol_DDROM = 6e-3;% 7e-3 (any value work) 	5.5124e-05

createABtilde_noconstraints;

ROMDriverPlotCE;

error_DDC = my_error;

WCEdataTableDDC = dataTableDDC;

% a=[1];
% b=[0.525];
% c=[0.5247];
% ,a,b,'w',a,c,'w';
% 
% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,4),'b-.',CEdataTableDDC(:,1),CEdataTableDDC(:,4),'r',WCEdataTableDDC(:,1),WCEdataTableDDC(:,4),'k',a,b,'w',a,c,'w','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% title(['r=' num2str(r)],'FontSize',20)
% K = legend('DNS','CE-DDC-ROM','DDC-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

error_DDC
error_CE_DDC
% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,3),'b-.',CEdataTableDDC(:,1),CEdataTableDDC(:,3),'r',WCEdataTableDDC(:,1),WCEdataTableDDC(:,3),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['r=' num2str(r)],'FontSize',20)
% K = legend('DNS','CEDDC-ROM','WCEDDC-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,2),'b-.',CEdataTableDDC(:,1),CEdataTableDDC(:,2),'r',WCEdataTableDDC(:,1),WCEdataTableDDC(:,2),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['r=' num2str(r)],'FontSize',20)
% K = legend('DNS','CEDDC-ROM','WCEDDC-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight


