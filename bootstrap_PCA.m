function [boot] = bootstrap_PCA(M0,H,I,L,R,MU0,MU1,MU2,F,q,k,p,scaling)

count=0;
for i=1:q
    i
    
    bootIndices = randsample(1:k,k,true);
    
    M0boot = M0(bootIndices,:);
    Hboot = H(bootIndices,:);
    Iboot = I(bootIndices,:);
    Lboot = L(bootIndices,:);
    Rboot = R(bootIndices,:);
    MU0boot = MU0(bootIndices,:);
    MU1boot = MU1(bootIndices,:);
    MU2boot = MU2(bootIndices,:);
    
    
    [data] = ALS_NCFPCA(M0boot,Hboot,Iboot,Lboot,Rboot,MU0boot,MU1boot,MU2boot,p,1e-6,500,scaling,0,0,0);
    
    if data.fail==0
        Fp = data.pars.F;
        T = inv(Fp'*Fp)*Fp'*F;
        Fp_star = Fp*T;
        PA(i) = 1-((norm(Fp_star-F)^2)/norm(F)^2);
        Fboot(i,:,:)=Fp_star';
        count=count+1;
        %input('')
    end
    
end


for b=1:size(Fboot,2)
    for c=1:size(Fboot,3)
        Fstar(c,b) = mean(Fboot(:,b,c));
        F_se(c,b) = std(Fboot(:,b,c));
        F_CI_lb(c,b) = quantile(Fboot(:,b,c),[.05])';
        F_CI_ub(c,b) = quantile(Fboot(:,b,c),[.95])';
    end
end

boot.F = Fboot;
boot.PA.dist = PA;
boot.PA.mean = mean(PA);
boot.PA.std= std(PA);
boot.Fstar = Fstar;
boot.se = F_se;
boot.CI.lb = F_CI_lb;
boot.CI.ub = F_CI_ub;
boot.conv=count;

end

