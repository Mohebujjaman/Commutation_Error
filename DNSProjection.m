clear all

load ROMCyl500_zt_103KSV_N30
load snapshots103K_SV_Re500_zt_dt002_u0_restart
%load snapshots103K_SV_Re500_zt 


%load snapshotData35Kdt002SV_Re100

%load ROMtestSV35K_N32_1500  Snapshots MassROM StiffROM  TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode

num_time_steps = size(Snapshots,2);

%num_time_steps = 10;
%size(GlobalV)
%GlobalVPrev = Snapshots(:,snapIndex-1);


N = 1;
DNSProjectionMatrix = zeros(N,num_time_steps);
MassROM = MassROM(1:N,1:N);
StiffROM = StiffROM(1:N,1:N);
%GradDivROM = GradDivROM(1:N,1:N);
%TriLinROM2 = TriLinROM2(1:N,1:N,1:N);
NLlift = NLlift(1:N,1:N);
NLdrag = NLdrag(1:N,1:N);
vdmass = vdmass(1:N);
vdstiff = vdstiff(1:N);
vlmass = vlmass(1:N);
vlstiff = vlstiff(1:N);



for ts=1:num_time_steps
    GlobalV = Snapshots(:,ts);

    vvv = GlobalV;
    RHS = zeros(N,1);
    for i=1:N
        RHS(i) = vvv' * (MassMatrix * PhiR(:,i) );
    end
    A = MassROM;
    RHS = RHS - A(:,1)*1;
    A(1,:)=0;
    A(:,1)=0;
    A(1,1)=1;
    RHS(1)=1;
    velts = A  \ RHS;

    DNSProjectionMatrix(:,ts) = velts;
end

save ('DNSProjectionMatrix_Re500_r1.mat','DNSProjectionMatrix')