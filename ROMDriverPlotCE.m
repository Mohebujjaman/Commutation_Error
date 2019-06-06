%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 7/20/2018                                    %
% Author: Jaman Mohebujjaman                   %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%load ROMtestSV35K_N16_166  Snapshots MassROM StiffROM  TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode
 
% load ROMCyl500_zt_103KSV_N30

%endTime = 0.3320; %
%endTime = 10;
GlobalV = Snapshots(:,snapIndex);

MassROM = MassROM(1:r,1:r);
StiffROM = StiffROM(1:r,1:r);

TriLinROM = TriLinROM(1:r,1:r,1:r);
NLlift = NLlift(1:r,1:r);
NLdrag = NLdrag(1:r,1:r);
vdmass = vdmass(1:r);
vdstiff = vdstiff(1:r);
vlmass = vlmass(1:r);
vlstiff = vlstiff(1:r);



% L2 project the initial condition (held in GlobalV) into ROM basis - coeff
% vector is put into "velInit"
vvv = GlobalV;
RHS = zeros(r,1);
 for i=1:r
     RHS(i) = vvv' * (MassMatrix * PhiR(:,i) );
 end
A = MassROM;
RHS = RHS - A(:,1)*1;
A(1,:)=0;
A(:,1)=0;
A(1,1)=1;
RHS(1)=1;
velInit = A  \ RHS;

velPrevPrev=velInit;


dataTableDDC=[];
dataTableDNS = [];
errorTable = [];

% simulation parameters

%T=endTime;
numTimeSteps =  round(endTime/dt);


solns = zeros(r,numTimeSteps+1);
solns(:,1) = velInit; 

%initial condition

 ts  = 1;
% 
 b = 1.0/dt * MassROM * velPrevPrev;
 %-my_factor*nu*weak_pfilter_CE(:,ts)-Gsnap(:,ts)
 % b = 1.0/dt * MassROM * velPrevPrev;
% % build matrix -0*weak_pfilter_CE(:,ts) 
%-data_driven_rom*nu*weak_dfilter_CE(:,ts)
 NLmat = 0*MassROM;
 for k=1:r
     %NLmat = NLmat + velPrevPrev(k)*(TriLinROM(:,:,k));
     NLmat = NLmat + velPrevPrev(k)*(TriLinROM(:,:,k)+data_driven_rom*BtildeDDC(:,:,k));
 end
  % A = 1.0/dt * MassROM + nu*StiffROM + NLmat +myAnew;
 %A = 1.0/dt * MassROM + nu*StiffROM + NLmat + data_driven_rom*AtildeDDC; 
 A = 1.0/dt * MassROM + nu*StiffROM + NLmat;     
% 
 A(1,:) = 0;
 A(1,1) = 1;
 b(1) = 1;
%     
% % solve the linear system
 velSoln = A \ b;
 


vt = 1/dt*(velSoln - velPrevPrev);
lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );

vtDNS = 1/dt * (DNSProjectionMatrix(:,snapIndex+ts)-DNSProjectionMatrix(:,snapIndex));
liftDNS = -20*( vtDNS'*vlmass' + nu*DNSProjectionMatrix(:,snapIndex+ts)'*vlstiff' + DNSProjectionMatrix(:,snapIndex+ts)' * NLlift * DNSProjectionMatrix(:,snapIndex+ts) );
energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,snapIndex+ts)' *(MassROM * DNSProjectionMatrix(:,snapIndex+ts)) ));
dragDNS = -20*( vtDNS'*vdmass' + nu*DNSProjectionMatrix(:,snapIndex+ts)'*vdstiff' + DNSProjectionMatrix(:,snapIndex+ts)' * NLdrag * DNSProjectionMatrix(:,snapIndex+ts) );


drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );
energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));

error = (velSoln-DNSProjectionMatrix(:,snapIndex+ts))' * (MassROM* (velSoln-DNSProjectionMatrix(:,snapIndex+ts)));

display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '   ' num2str(energy) ])
dataTableDDC = [dataTableDDC; ts*dt, lift, drag,   energy];
dataTableDNS = [dataTableDNS; ts*dt, liftDNS, dragDNS, energyDNS];
errorTable = [errorTable; error];
%      
%     
 solns(:,ts+1) = velSoln;
 velPrev = velSoln;



%BDF2 with linear scheme
   
for ts=2:numTimeSteps
    %RHS
    %-pfilter_CE(:,ts)
       
    b = 2/dt * MassROM * velPrev -0.5/dt*MassROM*velPrevPrev;
    % build matrix -0*weak_pfilter_CE(:,ts) -0*nu*weak_pfilter_CE(:,ts)
    % -my_factor*nu*weak_pfilter_CE(:,ts)-Gsnap(:,ts)
    NLmat = 0*MassROM;
    for k=1:r
        %NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM(:,:,k) );
       NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM(:,:,k)+data_driven_rom*BtildeDDC(:,:,k));
    end
    %A = 1.5/dt * MassROM + nu*StiffROM + NLmat;
    A = 1.5/dt * MassROM + nu*StiffROM + NLmat + data_driven_rom*AtildeDDC;
   % A = 1.5/dt * MassROM + nu*StiffROM + NLmat +myAnew;

    A(1,:)=0;
    A(1,1)=1;
    b(1)=1;
    
    % solve the linear system
    velSoln = A \ b;

       
     vt = 1/dt*(velSoln - velPrev);
     lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );
     drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );
     energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
     
     
     
  if(snapIndex==1000)
     
     if(ts<=1500)
     vtDNS = 1/dt * (DNSProjectionMatrix(:,snapIndex+ts)-DNSProjectionMatrix(:,snapIndex+ts-1));
     liftDNS = -20*( vtDNS'*vlmass' + nu*DNSProjectionMatrix(:,snapIndex+ts)'*vlstiff' + DNSProjectionMatrix(:,snapIndex+ts)' * NLlift * DNSProjectionMatrix(:,snapIndex+ts) );
     energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,snapIndex+ts)' *(MassROM * DNSProjectionMatrix(:,snapIndex+ts)) ));
     dragDNS = -20*( vtDNS'*vdmass' + nu*DNSProjectionMatrix(:,snapIndex+ts)'*vdstiff' + DNSProjectionMatrix(:,snapIndex+ts)' * NLdrag * DNSProjectionMatrix(:,snapIndex+ts) );
     dataTableDNS = [dataTableDNS; ts*dt,liftDNS, dragDNS, energyDNS];
     error = (velSoln-DNSProjectionMatrix(:,snapIndex+ts))' * (MassROM* (velSoln-DNSProjectionMatrix(:,snapIndex+ts)));
     errorTable = [errorTable; error];
     end
  else
      if(ts<=5000)
     vtDNS = 1/dt * (DNSProjectionMatrix(:,snapIndex+ts)-DNSProjectionMatrix(:,snapIndex+ts-1));
     liftDNS = -20*( vtDNS'*vlmass' + nu*DNSProjectionMatrix(:,snapIndex+ts)'*vlstiff' + DNSProjectionMatrix(:,snapIndex+ts)' * NLlift * DNSProjectionMatrix(:,snapIndex+ts) );
     energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,snapIndex+ts)' *(MassROM * DNSProjectionMatrix(:,snapIndex+ts)) ));
     dragDNS = -20*( vtDNS'*vdmass' + nu*DNSProjectionMatrix(:,snapIndex+ts)'*vdstiff' + DNSProjectionMatrix(:,snapIndex+ts)' * NLdrag * DNSProjectionMatrix(:,snapIndex+ts) );
     dataTableDNS = [dataTableDNS; ts*dt,liftDNS, dragDNS, energyDNS];
     error = (velSoln-DNSProjectionMatrix(:,snapIndex+ts))' * (MassROM * (velSoln-DNSProjectionMatrix(:,snapIndex+ts)));
     errorTable = [errorTable; error];
      end
  end
    
     display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '     ' num2str(energy) ])
     dataTableDDC = [dataTableDDC; ts*dt, lift, drag,   energy];
     
        
    solns(:,ts+1)=velSoln;
    velPrevPrev=velPrev;
    velPrev = velSoln;
    
end

my_error = sqrt(sum(errorTable)*dt/0.332);
% figure
% plot(dataTable1(:,1),dataTable1(:,3),'g', dataTable1(:,1),dataTable2(:,3),'r',dataTable1(:,1),dataTable3(:,3),'k', DNStable10(1:numTimeSteps,1)-7, DNStable10(1:numTimeSteps,4),'b-.','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% I = legend('G-ROM','iDD-ROM','DD-ROM', 'DNS')
% set(I,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% 
% 
% figure
% plot(dataTable1(:,1),dataTable1(:,2),'g',dataTable1(:,1),dataTable2(:,2),'r',dataTable1(:,1),dataTable3(:,2),'k',DNStable10(1:numTimeSteps,1)-7, DNStable10(1:numTimeSteps,5),'b-.','LineWidth',2)
% %plot(dataTable1(:,1),dataTable1(:,2)- BalanceTable1(end-1499:end,5),'k-',dataTable2(:,1),dataTable2(:,2)- BalanceTable2(end-1499:end,5),'r-.','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% J = legend('G-ROM','iDD-ROM','DD-ROM', 'DNS')
% set(J,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

%,



% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,2),'b-.',dataTableGROM(:,1),dataTableGROM(:,2),'g',dataTableDDF(:,1),dataTableDDF(:,2),'r',dataTableCDDF(:,1),dataTableCDDF(:,2),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% K = legend('DNS','G-ROM','DDF-ROM-quadratic','DDF-ROM-linear')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,3),'b-.',dataTableGROM(:,1),dataTableGROM(:,3),'g',dataTableDDF(:,1),dataTableDDF(:,3),'r',dataTableCDDF(:,1),dataTableCDDF(:,3),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% K = legend('DNS','G-ROM','DDF-ROM-quadratic','DDF-ROM-linear')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,4),'b-.',dataTableGROM(:,1),dataTableGROM(:,4),'g',dataTableDDF(:,1),dataTableDDF(:,4),'r',dataTableCDDF(:,1),dataTableCDDF(:,4),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% K = legend('DNS','G-ROM','DDF-ROM-quadratic','DDF-ROM-linear')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aa=3.002:0.002:5.974;
% appendd=[aa',dataTableDNS(14:end,2),dataTableDNS(14:end,3),dataTableDNS(14:end,4)];
% 
% dataTableDNS=[dataTableDNS;appendd];
% 
% aa=5.976:0.002:9.942;
% appendd=[aa',dataTableDNS(1004:end,2),dataTableDNS(1004:end,3),dataTableDNS(1004:end,4)];
% dataTableDNS=[dataTableDNS;appendd];
% 
% aa=9.944:0.002:10;
% 
% appendd=[aa',dataTableDNS(1501:1529,2),dataTableDNS(1501:1529,3),dataTableDNS(1501:1529,4)];
% dataTableDNS=[dataTableDNS;appendd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% aa=10.002:0.002:20;
% appendd=[aa',dataTableDNS(1:end,2),dataTableDNS(1:end,3),dataTableDNS(1:end,4)];
% dataTableDNS=[dataTableDNS;appendd];

% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,3),'b-.',dataTableCDDF(:,1),dataTableCDDF(:,3),'r',dataTableDDF(:,1),dataTableDDF(:,3),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['r=' num2str(r)],'FontSize',20)
% K = legend('DNS','CDDC-ROM','DDC-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,4),'b-.',dataTableDDC(:,1),dataTableDDC(:,4),'r','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% title(['r=' num2str(r)],'FontSize',20)
% K = legend('DNS','DDC-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

% figure
% plot(dataTableDNS(1:end,1), dataTableDNS(1:end,2),'b-.',dataTableCDDF(:,1),dataTableCDDF(:,2),'r',dataTableDDF(:,1),dataTableDDF(:,2),'k','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['r=' num2str(r)],'FontSize',20)
% K = legend('DNS','CDDC-ROM','DDC-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

%l=4500;

% display('For CDDC-ROM')
% dragaveCDD = sum(dataTableDDC(l:end,3))/length(dataTableDDC(l:end,3))
% diffdragCDD = max(dataTableDDC(l:end,3))-min(dataTableDDC(l:end,3))
% liftavgCDD = sum(dataTableDDC(l:end,2))/length(dataTableDDC(l:end,2))
% diffliftCDD = max(dataTableDDC(l:end,2))-min(dataTableDDC(l:end,2))
% 
% dragmaxCDD = max(dataTableDDC(l:end,3))
% dragminCDD = min(dataTableDDC(l:end,3))
% liftmaxCDD = max(dataTableDDC(l:end,2))
% liftminCDD = min(dataTableDDC(l:end,2))
% 
% 
% display('For DDC-ROM')
% dragaveDD = sum(dataTableDDF(l:end,3))/length(dataTableDDF(l:end,3))
% diffdragDD = max(dataTableDDF(l:end,3))-min(dataTableDDF(l:end,3))
% liftavgDD = sum(dataTableDDF(l:end,2))/length(dataTableDDF(l:end,2))
% diffliftDD = max(dataTableDDF(l:end,2))-min(dataTableDDF(l:end,2))
% dragmaxDD = max(dataTableDDF(l:end,3))
% dragminDD = min(dataTableDDF(l:end,3))
% liftmaxDD = max(dataTableDDF(l:end,2))
% liftminDD = min(dataTableDDF(l:end,2))
% 
% 
% display('For DNS')
% dragaveDNS = sum(dataTableDNS(l:end,3))/length(dataTableDNS(l:end,3))
% diffdragDNS = max(dataTableDNS(l:end,3))-min(dataTableDNS(l:end,3))
% liftavgDNS = sum(dataTableDNS(l:end,2))/length(dataTableDNS(l:end,2))
% diffliftDNS = max(dataTableDNS(l:end,2))-min(dataTableDNS(l:end,2))
% dragmaxDNS = max(dataTableDNS(l:end,3))
% dragminDNS = min(dataTableDNS(l:end,3))
% liftmaxDNS = max(dataTableDNS(l:end,2))
% liftminDNS = min(dataTableDNS(l:end,2))
% 
% 
% display('ERROR-CDDC')
% er_drag_CDD = abs(dragaveCDD-dragaveDNS)/abs(dragaveDNS)*100
% er_drag_range_CDD = abs(diffdragCDD-diffdragDNS)/abs(diffdragDNS)*100
% er_lift_CDD = abs(liftavgCDD-liftavgDNS)/abs(liftavgDNS)*100
% er_lift_range_CDD = abs(diffliftCDD-diffliftDNS)/abs(diffliftDNS)*100
% 
% er_max_drag_cdd= abs(dragmaxCDD-dragmaxDNS)/abs(dragmaxDNS)*100
% er_min_drag_cdd= abs(dragminCDD-dragminDNS)/abs(dragminDNS)*100
% er_max_lift_cdd= abs(liftmaxCDD-liftmaxDNS)/abs(liftmaxDNS)*100
% er_min_lift_cdd= abs(liftminCDD-liftminDNS)/abs(liftminDNS)*100
% 
% display('ERROR-DDC')
% er_drag_DD = abs(dragaveDD-dragaveDNS)/abs(dragaveDNS)*100
% er_drag_range_DD = abs(diffdragDD-diffdragDNS)/abs(diffdragDNS)*100
% er_lift_DD = abs(liftavgDD-liftavgDNS)/abs(liftavgDNS)*100
% er_lift_range_DD = abs(diffliftDD-diffliftDNS)/abs(diffliftDNS)*100
% 
% er_max_drag_dd= abs(dragmaxDD-dragmaxDNS)/abs(dragmaxDNS)*100
% er_min_drag_dd= abs(dragminDD-dragminDNS)/abs(dragminDNS)*100
% er_max_lift_dd= abs(liftmaxDD-liftmaxDNS)/abs(liftmaxDNS)*100
% er_min_lift_dd= abs(liftminDD-liftminDNS)/abs(liftminDNS)*100
% 


% figure
% plot(dataTable1(:,1),dataTable1(:,4),'g', dataTable1(:,1), dataTableDNS(:,2),'b',dataTable1(:,1),dataTable2(:,4),'r',a,bmm,'k','LineWidth',2)
% %plot(dataTable1(:,1), dataTableDNS(:,2),'b',dataTable1(:,1),dataTable2(:,4),'r',dataTable1(:,1),dataTable3(:,4),'k','LineWidth',2)
% %plot(dataTable1(:,1),dataTable1(:,5)- BalanceTable1(end-1499:end,2),'k-',dataTable2(:,1),dataTable2(:,5)- BalanceTable2(end-1499:end,2),'r-.','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% %title(['N=' num2str(N)],'FontSize',20)
% K = legend('G-ROM','DNS','DDF-ROM + Constraints','DDF-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

