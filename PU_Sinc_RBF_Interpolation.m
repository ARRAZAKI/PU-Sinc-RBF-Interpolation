clear; close all; clc;
Nlist = [1250,2500,5000]; % 10000, 20000 % ajouter 5000 pour etude de convergence
nRepeats = 1;               
kernels = {'sinc'};% ,'mq','gaussian'
modes = {'loocv'};%,'loocv_cond','LOOCV-GA'
cond_limit = 1e12;         % seuil d'etude 
eps_bounds = [5,15];       % domaine historique; tester une plage plus large si necessaire
s_d = 60;                  % 3600 noeuds de validation fixes
add_data_perturbation = false; 
Results = table();

d = 2;
bar_x = MakeSDGrid(d,s_d);
truebar = f2(bar_x(:,1),bar_x(:,2)); % utilise uniquement pour EVALUER l'erreur
w = @(supp,r) max(1-supp*r,0).^4.*(4*supp*r+1);
for irep = 1:nRepeats
    for iN = 1:numel(Nlist)
        N = Nlist(iN);
        hs = haltonset(d);
        if irep>1
            hs.Skip = (irep-1)*N;
            hs = scramble(hs,'RR2');
        end
        x = net(hs,N);
        f = f2(x(:,1),x(:,2));
        m_d = max(2,floor(sqrt(N)/4));
        for ik = 1:numel(kernels)
            kernel = kernels{ik};
            switch kernel
                case 'sinc'
                    phi = @(ep,X,Y) sinc(ep*DistanceMatrix(X,Y));
                case 'mq'
                    phi = @(ep,X,Y) sqrt(1+(ep*DistanceMatrix(X,Y)).^2);
                case 'gaussian'
                    phi = @(ep,X,Y) exp(-(ep*DistanceMatrix(X,Y)).^2);
            end
            for im = 1:numel(modes)
                mode = modes{im};
                if ~strcmp(kernel,'sinc') && ~strcmp(mode,'loocv')
                    continue;
                end
                t0=tic;
                [Pf,diagPU] = PU_SINC_RBF( ...
                    d,x,bar_x,m_d,phi,w,f,[],mode,cond_limit,eps_bounds);
                tsec=toc(t0);
                err=Pf-truebar;
                rmse=sqrt(mean(err.^2)); linf=max(abs(err));
                cv=diagPU.cond_patch; cv=cv(~isnan(cv));
                if isempty(cv)
                    condMean=NaN; condMax=NaN; badFraction=NaN;
                else
                    condMean=mean(cv); condMax=max(cv);
                    badFraction=mean(cv>cond_limit);
                end
                ev=diagPU.eps_patch;
                ev=ev(isfinite(ev));
                if isempty(ev)
                    epsMean=NaN; epsMin=NaN; epsMax=NaN;
                else
                    epsMean=mean(ev); epsMin=min(ev); epsMax=max(ev);
                end
                dataAmplification=NaN;
                if add_data_perturbation
                    df=1e-8*max(1,norm(f,inf))*sin((1:N)');
                    [Pfp,~] = PU_SINC_RBF( ...
                        d,x,bar_x,m_d,phi,w,f+df,diagPU.eps_patch,mode,cond_limit,eps_bounds);
                    relU=norm(Pfp-Pf,2)/max(norm(Pf,2),eps);
                    relD=norm(df,2)/max(norm(f,2),eps);
                    dataAmplification=relU/max(relD,eps);
                end
                fallback=sum(diagPU.fallback_patch);
                row=table({kernel},{mode},irep,N,rmse,linf,tsec,epsMean, ...
                    epsMin,epsMax,condMean,condMax,badFraction,fallback,dataAmplification, ...
                    'VariableNames',{'Kernel','Mode','Repeat','N','RMSE','Linf','Time_s', ...
                    'Eps_mean','Eps_min','Eps_max','Cond_mean','Cond_max', ...
                    'Bad_fraction','Fallback_count','Data_amplification'});
                Results=[Results;row];
                fprintf('%s %-11s repeat=%d N=%d RMSE=%9.2e LINF=%.3e condMax=%9.2e CPU=%.2f s\n', ...
                    kernel,mode,irep,N,rmse,linf,condMax,tsec);
            end
        end
    end
end

S = groupsummary(Results,{'Kernel','Mode','N'}, ...
    {'mean','std','max'}, ...
    {'Eps_mean','RMSE','Linf','Time_s','Cond_mean','Cond_max','Data_amplification'});

T = S(:,{'Kernel','Mode','N','mean_Eps_mean','mean_RMSE', ...
    'std_RMSE','mean_Linf','mean_Time_s','mean_Cond_mean', ...
    'max_Cond_max','mean_Data_amplification'});

T.Properties.VariableNames = {'RBF','Mode','N','Eps_mean', ...
    'RMSE','RMSE_std','LINF','Time_s','Cond_mean', ...
    'Cond_max','Data_amplification'};

disp(T);
writetable(T,'Interpolation_Results.csv');
