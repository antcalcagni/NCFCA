function [data] = ALS_NCFPCA(M0,H,I,L,R,MU0,MU1,MU2,p,tolMax,maxIter,scaling,varimax_cond,plotFN,plotRES,axisPlot,k_rect)

m = size(M0,2);
n = size(M0,1);
Im = eye(m,m);
iter=0;
fail=0;

%% Pre-processing procedure
M1 = M0-H; 
M2 = M0+I;
switch scaling
    case 1
        MU0 = MU0*(diag((std(M0))))^-1;
        MU1 = MU1*(diag((std(M1))))^-1;
        MU2 = MU2*(diag((std(M2))))^-1;
        H = H*(diag((std(M0))))^-1;
        I = I*(diag((std(M0))))^-1;
        L = L*(diag((std(M1))))^-1;
        R = R*(diag((std(M2))))^-1;
        M0 = (M0-ones(n,m)*diag(mean(M0)))*(diag((std(M0))))^-1;
    case 2
        MU0 = MU0*(diag(sqrt(std(M0))))^-1;
        MU1 = MU1*(diag(sqrt(std(M1))))^-1;
        MU2 = MU2*(diag(sqrt(std(M2))))^-1;
        H = H*(diag(sqrt(std(M0))))^-1;
        I = I*(diag(sqrt(std(M0))))^-1;
        L = L*(diag(sqrt(std(M1))))^-1;
        R = R*(diag(sqrt(std(M2))))^-1;
        M0 = (M0-ones(n,m)*diag(mean(M0)))*(diag(sqrt(std(M0))))^-1;
end 

%% Generating Starting Values
[Um Sm Vm] = svd(M0);
Am = Um*Sm(:,1:p);

[Uh Sh Vh] = svd(H);
Ah = Uh*Sh(:,1:p);

[Ui Si Vi] = svd(I);
Ai = Ui*Si(:,1:p);

[Ul Sl Vl] = svd(L);
Al = Ul*Sl(:,1:p);

[Ur Sr Vr] = svd(R);
Ar = Ur*Sr(:,1:p);

[Umu0 Smu0 Vmu0] = svd(MU0);
Amu0 = Umu0*Smu0(:,1:p);

[Umu1 Smu1 Vmu1] = svd(MU1);
Amu1 = Umu1*Smu1(:,1:p);

[Umu2 Smu2 Vmu2] = svd(MU2);
Amu2 = Umu2*Smu2(:,1:p);

f = (Vm+Vh+Vi+Vl+Vr+Vmu0+Vmu1+Vmu2)/8;
F = f(:,1:p);

%% ALS loop
for i=1:maxIter
    
    Al_old=Al;Ar_old=Ar;Ah_old=Ah;Ai_old=Ai;F_old=F;
    Am_old=Am;Amu0_old=Amu0;Amu1_old=Amu1;Amu2_old=Amu2;
    
    %%% Computing alternating estimations
    Al = (2^m)*(H*F - Ah*F'*F + L*F)*inv((2^m)*F'*F);
    Ar = (2^m)*(I*F - Ai*F'*F + R*F)*inv((2^m)*F'*F);
    Ai = (2^m)*(R*F - Ar*F'*F + I*F)*inv((2^m)*F'*F);
    Ah = (2^m)*(L*F - Al*F'*F + H*F)*inv((2^m)*F'*F);
    Am = (M0*F)*inv(F'*F);
    Amu0 = (MU0*F)*inv(F'*F);
    Amu1 = (MU1*F)*inv(F'*F);
    Amu2 = (MU2*F)*inv(F'*F);
    F = inv( ((2^m)*(kron(Ah'*Ah,Im) + kron(Al'*Al,Im) + kron(Ai'*Ai,Im) + kron(Ar'*Ar,Im))) + ((2^m)*(kron(Am'*Am,Im) + kron(Amu0'*Amu0,Im) + kron(Amu1'*Amu1,Im) + kron(Amu2'*Amu2,Im))) + ((2^(m))*(kron(Ah'*Al,Im) + kron(Al'*Ah,Im) + kron(Ai'*Ar,Im) + kron(Ar'*Ai,Im))))*vec( ((2^m)*(H'*Ah + L'*Al  + I'*Ai + R'*Ar)) + ((2^m)*(M0'*Am + MU0'*Amu0 + MU1'*Amu1 + MU2'*Amu2)) + ((2^m)*(H'*Al + L'*Ah + I'*Ar + R'*Ai)));
    F = reshape(F,m,p);
    
    %%% Computing convergence criteria 
    k1(i) = norm(Al_old-Al)^2;
    k2(i) = norm(Ar_old-Ar)^2;
    k3(i) = norm(Am_old-Am)^2;
    k4(i) = norm(Amu1_old-Amu1)^2;
    k5(i) = norm(Amu2_old-Amu2)^2;
    k6(i) = norm(Amu0_old-Amu0)^2;
    k7(i) = norm(F_old-F)^2;
    k8(i) = norm(Ah_old-Ah)^2;
    k9(i) = norm(Ai_old-Ai)^2;
    CP(i,:) = [k1(i) k2(i) k3(i) k4(i) k5(i) k6(i) k7(i) k8(i) k9(i)];
    
    %%% Computing estimated components 
    M0star = Am*F';
    Lstar = Al*F';
    Hstar = Ah*F';
    Istar = Ai*F';
    Rstar = Ar*F';
    MU0star = Amu0*F';
    MU1star = Amu1*F';
    MU2star = Amu2*F';
    
    %%% Stopping rule 
    foVal(i) = (2^(m-1)*norm(M0-M0star)^2) + (2^(m-1)*norm(I-Istar)^2) + (2^(m-1)*norm(H-Hstar)^2) + (2^(m-1)*norm(L-Lstar)^2) + (2^(m-1)*norm(R-Rstar)^2) + (2^(m)*trace((I-Istar)'*(R-Rstar))) + (2^(m)*trace((H-Hstar)'*(L-Lstar))) + (2^(m-1)*(norm(MU0-MU0star)^2 + norm(MU1-MU1star)^2 + norm(MU2-MU2star)^2));
    if i > 2
        if abs(foVal(i-1)-foVal(i)) < tolMax
            iter=i;
            break
        end
    end

end

if iter==0
    fail=1;
    disp('ALS did not converge at a global/local minimum')
else
    disp('ALS converged at a global/local minimum')
end

%% Goodness-of-fit of the model
data.fit.RSS = (norm(M0-M0star)^2 + norm(I-Istar)^2 + norm(H-Hstar)^2 + norm(L-Lstar)^2 + norm(R-Rstar)^2 + 2*trace((I-Istar)'*(R-Rstar)) + 2*trace((H-Hstar)'*(L-Lstar)) + norm(MU0-MU0star)^2 + norm(MU1-MU1star)^2 + norm(MU2-MU2star)^2);
data.fit.TSS = (norm(M0)^2 + norm(I)^2 + norm(H)^2 + norm(L)^2 + norm(R)^2 + 2*trace((I)'*(R)) + 2*trace((H)'*(L)) + norm(MU0)^2 + norm(MU1)^2 + norm(MU2)^2);
data.fit.ESS = data.fit.TSS - data.fit.RSS;
data.fit.R2 = 1-(data.fit.RSS/data.fit.TSS);

%% Saving final data
data.pars.Am = Am;
data.pars.Ah = Ah;
data.pars.Ai = Ai;
data.pars.Al = Al;
data.pars.Ar = Ar;
data.pars.Amu0 = Amu0;
data.pars.Amu1 = Amu1;
data.pars.Amu2 = Amu2;
data.pars.F = F;

data.estim.M0 = M0star;
data.estim.H = Hstar;
data.estim.I = Istar;
data.estim.L = Lstar;
data.estim.R = Rstar;
data.estim.MU0 = MU0star;
data.estim.MU1 = MU1star;
data.estim.MU2 = MU2star;

data.CP = CP;
data.FoVal = foVal;
data.TotIter = iter;
data.fail = fail;

%% Deleting null values from estimated components
data.estim.L(data.estim.L<0)=0;
data.estim.H(data.estim.H<0)=0;
data.estim.I(data.estim.I<0)=0;
data.estim.R(data.estim.R<0)=0;

data.estim.MU1(data.estim.H==0)=data.estim.MU0(data.estim.H==0);
data.estim.MU2(data.estim.I==0)=data.estim.MU0(data.estim.I==0);


%% Normalizing MU0,MU1,MU2
[data.estim.MU0,data.estim.MU1,data.estim.MU2] = normalizeMU(data.estim.MU0,data.estim.MU1,data.estim.MU2);

%% Modified Gram-Schmidt Orthonormalization for estimated loadings and scores (with Varimax Rotation)
[F_ort T]= MGMa(data.pars.F);

data.ort.F = F_ort;
data.ort.Am = data.pars.Am*inv(T');
data.ort.Ah = data.pars.Ah*inv(T');
data.ort.Ai = data.pars.Ai*inv(T');
data.ort.Al = data.pars.Al*inv(T');
data.ort.Ar = data.pars.Ar*inv(T');
data.ort.Amu0 = data.pars.Amu0*inv(T');
data.ort.Amu1 = data.pars.Amu1*inv(T');
data.ort.Amu2 = data.pars.Amu2*inv(T');

switch varimax_cond
    case 1
    [F_ort T]= rotatefactors(F_ort);
    data.ort.F = F_ort;
    data.ort.Am = data.ort.Am*inv(T');
    data.ort.Ah = data.ort.Ah*inv(T');
    data.ort.Ai = data.ort.Ai*inv(T');
    data.ort.Al = data.ort.Al*inv(T');
    data.ort.Ar = data.ort.Ar*inv(T');
    data.ort.Amu0 = data.ort.Amu0*inv(T');
    data.ort.Amu1 = data.ort.Amu1*inv(T');
    data.ort.Amu2 = data.ort.Amu2*inv(T');
end


%% Graphical output for convergence criteria
close all
switch plotFN
    case 1
        figure();
        subplot(2,5,1);plot(foVal,'r');title('FoVal');
        subplot(2,5,2);plot(k3);title('Am');axis([0 i+3 min(k3)-(2*mean(k3)) max(k3)])
        subplot(2,5,3);plot(k8);title('Ah');axis([0 i+3 min(k8)-(2*mean(k8)) max(k8)])
        subplot(2,5,4);plot(k1);title('Al');axis([0 i+3 min(k1)-(2*mean(k1)) max(k1)])
        subplot(2,5,5);plot(k9);title('Ai');axis([0 i+3 min(k9)-(2*mean(k9)) max(k9)])
        subplot(2,5,6);plot(k2);title('Ar');axis([0 i+3 min(k2)-(2*mean(k2)) max(k2)])
        subplot(2,5,7);plot(k6);title('Amu0');axis([0 i+3 min(k6)-(2*mean(k6)) max(k6)])
        subplot(2,5,8);plot(k4);title('Amu1');axis([0 i+3 min(k4)-(2*mean(k4)) max(k4)])
        subplot(2,5,9);plot(k5);title('Amu2');axis([0 i+3 min(k5)-(2*mean(k5)) max(k5)])
        subplot(2,5,10);plot(k7);title('F');axis([0 i+3 min(k7)-(2*mean(k7)) max(k7)])
end

%c = {'Piedmont', 'Aosta_Valley', 'Lombardy', 'Veneto', 'Friuli', 'Liguria', 'Emilia_Romagna', 'Tuscany', 'Umbria', 'Marche', 'Lazio','Abruzzo', 'Molise', 'Campania', 'Apulia', 'Basilicata', 'Calabria', 'Sicily', 'Sardinia', 'Trentino'};
%c = {'FRA1', 'FRA2', 'FRA3', 'HUS1', 'HUS2', 'HUS3', 'INC1', 'INC2', 'INC3', 'ISA1', 'ISA2', 'ISA3', 'JPL1', 'JPL2', 'JPL3', 'KHA1', 'KHA2', 'KHA3', 'LOT1', 'LOT2', 'LOT3', 'PHI1', 'PHI2', 'PHI3', 'ROM1', 'ROM2', 'ROM3'};
%c = {'LIN', 'PER', 'COT', 'SES', 'CAM', 'OLI', 'BEE', 'HOG'};
%c = {'NewEng', 'MidAtl', 'MidWestE', 'MidWestW', 'SouthAtl', 'SouthE', 'SouthW', 'Mountain', 'Pacific'};

%% Plot of results
switch plotRES
    case 1
    %%% Score Plot
    figure();hold on
    rect_colors = rand(n,3);
    for i=1:size(data.ort.Am)
        x_M0 = data.ort.Am(i,axisPlot(1));
        y_M0 = data.ort.Am(i,axisPlot(2));
        %plot(x_M0, y_M0,'r.')

        x_LB = (x_M0-k_rect*abs(data.ort.Ah(i,axisPlot(1)))-k_rect*abs(data.ort.Al(i,axisPlot(1))));
        y_LB = (y_M0-k_rect*abs(data.ort.Ah(i,axisPlot(2)))-k_rect*abs(data.ort.Al(i,axisPlot(2))));
        %plot(x_LB, y_LB,'k.')

        x_UB = (x_M0+k_rect*abs(data.ort.Ai(i,axisPlot(1)))+k_rect*abs(data.ort.Ar(i,axisPlot(1))));
        y_UB = (y_M0+k_rect*abs(data.ort.Ai(i,axisPlot(2)))+k_rect*abs(data.ort.Ar(i,axisPlot(2))));
        %plot(x_UB, y_UB,'k.')
        
        rectangle('Position',[x_LB y_LB abs(x_LB-x_UB) abs(y_LB-y_UB)],'edgecolor',rect_colors(i,:)) %,'FaceColor',rgb('White')
        rectangle('Position',[x_M0 y_M0 abs(x_UB-x_M0) abs(y_UB-y_M0)],'LineStyle','--','edgecolor',rect_colors(i,:))
        rectangle('Position',[x_LB y_LB abs(x_LB-x_M0) abs(y_LB-y_M0)],'LineStyle','--','edgecolor',rect_colors(i,:))
        colormap(rect_colors);colorbar

        text(x_LB,y_LB,['n.' num2str(i,'%d')],'FontSize',14,'Color',rect_colors(i,:))
        %text(x_LB,y_LB,c(i),'FontSize',12,'Color',rect_colors(i,:))
        
    end
    axescenter()

    %%% Loadings plot
    figure();hold on
    scatter(data.ort.F(:,axisPlot(1)),data.ort.F(:,axisPlot(2)),'.k')
    for j=1:m
        a=data.ort.F(j,axisPlot(1));
        b=data.ort.F(j,axisPlot(2));
        text(a,b,['var.' num2str(j,'%d')],'FontSize',15)
    end
    axescenter()
    
end



end