    idxPartition2 = floor(size(pp,2)/2);
    idxPartition3 = floor(size(pp,2)/3);
    [AIC1_2,polycoefs1_2,b_o1_2,n_params1_2,ne1_2,neEnd1_2,Ie1_2,exitflag1_2] = AICcFitParSeq(pp(:,1:idxPartition2),...
                                                      ppstd(:,1:idxPartition2),...
                                                      alpha(:,1:idxPartition2),...
                                                      dt(1:idxPartition2),...
                                                      ne00,...
                                                      A,...
                                                      polycoefs(:,1:idxPartition2),...
                                                      [],... % Was: ieprior,...
                                                      [], ... % was: stdprior,...
                                                      Ie0(:,end),...
                                                      Directives);
    [AIC2_2,polycoefs2_2,b_o2_2,n_params2_2,ne2_2,neEnd2_2,Ie2_2,exitflag2_2] = AICcFitParSeq(pp(:,(idxPartition2+1):end),...
                                                      ppstd(:,(idxPartition2+1):end),...
                                                      alpha(:,(idxPartition2+1):end),...
                                                      dt((idxPartition2+1):end),...
                                                      neEnd1_2(:,end),...
                                                      A,...
                                                      polycoefs(:,(idxPartition2+1):end),...
                                                      [],... % was: ieprior,...
                                                      [],... % was: stdprior(),...
                                                      Ie1_2(:,end),...
                                                      Directives);
    AIChalves = AICc( pp, ...
                      ppstd.^2, ...
                      [ne1_2,ne2_2], ...
                      n_params1_2(1) + n_params2_2(1), ...
                      numel(ne1_2) + numel(ne2_2), ...
                      Directives.ErrType, ...
                      Directives.ErrWidth );