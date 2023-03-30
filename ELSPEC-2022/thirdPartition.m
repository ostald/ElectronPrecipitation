    idxPartition3 = floor(size(pp,2)/3);
    if idxPartition3 == 0
      AICthirds  = inf; % just to avoid splitting periods of 2 profiles
    else
    [AIC1_3,polycoefs1_3,b_o1_3,n_params1_3,ne1_3,neEnd1_3,Ie1_3,exitflag1_3] = AICcFitParSeq(pp(:,1:idxPartition3),...
                                                      ppstd(:,1:idxPartition3),...
                                                      alpha(:,1:idxPartition3),...
                                                      dt(1:idxPartition3),...
                                                      ne00,...
                                                      A,...
                                                      polycoefs(:,1:idxPartition3),...
                                                      [],... % Was: ieprior,...
                                                      [], ... % was: stdprior,...
                                                      Ie0(:,end),...
                                                      Directives);
    [AIC2_3,polycoefs2_3,b_o2_3,n_params2_3,ne2_3,neEnd2_3,Ie2_3,exitflag2_3] = AICcFitParSeq(pp(:,(idxPartition3+1):(2*idxPartition3)),...
                                                      ppstd(:,(idxPartition3+1):(2*idxPartition3)),...
                                                      alpha(:,(idxPartition3+1):(2*idxPartition3)),...
                                                      dt((idxPartition3+1):(2*idxPartition3)),...
                                                      neEnd1_3(:,end),...
                                                      A,...
                                                      polycoefs(:,(idxPartition3+1):(2*idxPartition3)),...
                                                      [],... % was: ieprior,...
                                                      [],... % was: stdprior(),...
                                                      Ie1_3(:,end),...
                                                      Directives);
    [AIC3_3,polycoefs3_3,b_o3_3,n_params3_3,ne3_3,neEnd3_3,Ie3_3,exitflag3_3] = AICcFitParSeq(pp(:,(2*idxPartition3+1):end),...
                                                      ppstd(:,(2*idxPartition3+1):end),...
                                                      alpha(:,(2*idxPartition3+1):end),...
                                                      dt((2*idxPartition3+1):end),...
                                                      neEnd2_3(:,end),...
                                                      A,...
                                                      polycoefs(:,(2*idxPartition3+1):end),...
                                                      [],... % was: ieprior,...
                                                      [],... % was: stdprior(),...
                                                      Ie2_3(:,end),...
                                                      Directives);
    AICthirds = AICc( pp, ...
                      ppstd.^2, ...
                      [ne1_3,ne2_3,ne3_3], ...
                      n_params1_3(1) + n_params2_3(1)+ n_params3_3(1), ...
                      numel(ne1_3) + numel(ne2_3) + numel(ne3_3), ...
                      Directives.ErrType, ...
                      Directives.ErrWidth );
    end