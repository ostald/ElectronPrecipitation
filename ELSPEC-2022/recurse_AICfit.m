function [ne,neEnd,Ie,polycoefs,best_order,n_params,exitflag,nSteps,idx_out] = recurse_AICfit(ne,neEnd,pp,ppstd,alpha,dt,ne00,A,polycoefs,best_order,n_params,Ie,idx_in,exitflag,ieprior,stdprior,Ie0,Directives)
%                                                                                                                        ne00 a priori electron density
% RECURSE_AICFIT - maybe better to scrap AICFull in the outputs.
%   
%  savename = sprintf('in-%03i-%03i.mat',idx_in(1),numel(idx_in));
%  save(savename,'ne00')
  nSteps = numel(dt).*ones(size(dt(:)'));
  % disp(['Entering recurse_AICfit, nt: ',num2str(numel(dt))])
  AICFull = AICc( pp, ...
                  ppstd.^2, ...
                  ne, ...
                  n_params(1), ...
                  numel(ne), ...
                  Directives.ErrType, ...
                  Directives.ErrWidth );
  if numel(dt) == 1
    % Nothing more to be done just abandon mission, we have reached the
    % end of the branch
    idx_out = idx_in;
  else
    % Then there splitting-work to be tested
    % First we try breaking into halves

    
    j1 = batch('halfPartition');

%     idxPartition2 = floor(size(pp,2)/2);
%     idxPartition3 = floor(size(pp,2)/3);
%     [AIC1_2,polycoefs1_2,b_o1_2,n_params1_2,ne1_2,neEnd1_2,Ie1_2,exitflag1_2] = AICcFitParSeq(pp(:,1:idxPartition2),...
%                                                       ppstd(:,1:idxPartition2),...
%                                                       alpha(:,1:idxPartition2),...
%                                                       dt(1:idxPartition2),...
%                                                       ne00,...
%                                                       A,...
%                                                       polycoefs(:,1:idxPartition2),...
%                                                       [],... % Was: ieprior,...
%                                                       [], ... % was: stdprior,...
%                                                       Ie0(:,end),...
%                                                       Directives);
%     [AIC2_2,polycoefs2_2,b_o2_2,n_params2_2,ne2_2,neEnd2_2,Ie2_2,exitflag2_2] = AICcFitParSeq(pp(:,(idxPartition2+1):end),...
%                                                       ppstd(:,(idxPartition2+1):end),...
%                                                       alpha(:,(idxPartition2+1):end),...
%                                                       dt((idxPartition2+1):end),...
%                                                       neEnd1_2(:,end),...
%                                                       A,...
%                                                       polycoefs(:,(idxPartition2+1):end),...
%                                                       [],... % was: ieprior,...
%                                                       [],... % was: stdprior(),...
%                                                       Ie1_2(:,end),...
%                                                       Directives);
%     AIChalves = AICc( pp, ...
%                       ppstd.^2, ...
%                       [ne1_2,ne2_2], ...
%                       n_params1_2(1) + n_params2_2(1), ...
%                       numel(ne1_2) + numel(ne2_2), ...
%                       Directives.ErrType, ...
%                       Directives.ErrWidth );


    % Then we try to split the intervall into thirds

    j2 = batch('thirdPartition');

%     idxPartition3 = floor(size(pp,2)/3);
%     if idxPartition3 == 0
%       AICthirds  = inf; % just to avoid splitting periods of 2 profiles
%     else
%     [AIC1_3,polycoefs1_3,b_o1_3,n_params1_3,ne1_3,neEnd1_3,Ie1_3,exitflag1_3] = AICcFitParSeq(pp(:,1:idxPartition3),...
%                                                       ppstd(:,1:idxPartition3),...
%                                                       alpha(:,1:idxPartition3),...
%                                                       dt(1:idxPartition3),...
%                                                       ne00,...
%                                                       A,...
%                                                       polycoefs(:,1:idxPartition3),...
%                                                       [],... % Was: ieprior,...
%                                                       [], ... % was: stdprior,...
%                                                       Ie0(:,end),...
%                                                       Directives);
%     [AIC2_3,polycoefs2_3,b_o2_3,n_params2_3,ne2_3,neEnd2_3,Ie2_3,exitflag2_3] = AICcFitParSeq(pp(:,(idxPartition3+1):(2*idxPartition3)),...
%                                                       ppstd(:,(idxPartition3+1):(2*idxPartition3)),...
%                                                       alpha(:,(idxPartition3+1):(2*idxPartition3)),...
%                                                       dt((idxPartition3+1):(2*idxPartition3)),...
%                                                       neEnd1_3(:,end),...
%                                                       A,...
%                                                       polycoefs(:,(idxPartition3+1):(2*idxPartition3)),...
%                                                       [],... % was: ieprior,...
%                                                       [],... % was: stdprior(),...
%                                                       Ie1_3(:,end),...
%                                                       Directives);
%     [AIC3_3,polycoefs3_3,b_o3_3,n_params3_3,ne3_3,neEnd3_3,Ie3_3,exitflag3_3] = AICcFitParSeq(pp(:,(2*idxPartition3+1):end),...
%                                                       ppstd(:,(2*idxPartition3+1):end),...
%                                                       alpha(:,(2*idxPartition3+1):end),...
%                                                       dt((2*idxPartition3+1):end),...
%                                                       neEnd2_3(:,end),...
%                                                       A,...
%                                                       polycoefs(:,(2*idxPartition3+1):end),...
%                                                       [],... % was: ieprior,...
%                                                       [],... % was: stdprior(),...
%                                                       Ie2_3(:,end),...
%                                                       Directives);
%     AICthirds = AICc( pp, ...
%                       ppstd.^2, ...
%                       [ne1_3,ne2_3,ne3_3], ...
%                       n_params1_3(1) + n_params2_3(1)+ n_params3_3(1), ...
%                       numel(ne1_3) + numel(ne2_3) + numel(ne3_3), ...
%                       Directives.ErrType, ...
%                       Directives.ErrWidth );

    wait(j1);
    wait(j2);
    load(j1);
    load(j2);
    %end
    
    if AIChalves < AICFull & AIChalves <= AICthirds
      % Divide into 2 halves
      [ne1_2,neEnd1_2,Ie1_2,polycoefs1_2,best_order1_2,n_params1_2,exitflag1_2,nSteps1_2,idx_out1_2] = recurse_AICfit(ne1_2,...
                                                        neEnd1_2,...% was:(:,1:idxPartition2),...
                                                        pp(:,1:idxPartition2),...
                                                        ppstd(:,1:idxPartition2),...
                                                        alpha(:,1:idxPartition2),...
                                                        dt(1:idxPartition2),...
                                                        ne00,...
                                                        A,...
                                                        polycoefs1_2,...%was:(:,:,1:idxPartition2),...
                                                        b_o1_2,...%was:best_order(1:idxPartition2),...
                                                        n_params1_2,...%was(1:idxPartition2),...
                                                        Ie1_2,...% was(:,1:idxPartition2),...
                                                        idx_in(1:idxPartition2),...
                                                        exitflag1_2,...%was: (1:idxPartition2),...
                                                        [], ... % was: ieprior,...
                                                        [], ... % was: stdprior,...
                                                        Ie0(:,end),...
                                                        Directives);
      % Then we need to re-search for the best single-parameters for the
      % second half
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
      % and continue on with this into the recursion for the second half
      [ne2_2,neEnd2_2,Ie2_2,polycoefs2_2,best_order2_2,n_params2_2,exitflag2_2,nSteps2_2,idx_out2_2] = recurse_AICfit(ne2_2,...
                                                        neEnd2_2,...% was: (:,(idxPartition2+1):end),...
                                                        pp(:,(idxPartition2+1):end),...
                                                        ppstd(:,(idxPartition2+1):end),...
                                                        alpha(:,(idxPartition2+1):end),...
                                                        dt((idxPartition2+1):end),...
                                                        neEnd1_2(:,end),... % 3## FIX!
                                                        A,...
                                                        polycoefs2_2,...%was:(:,:,(idxPartition2+1):end),...
                                                        b_o2_2,...%was:best_order((idxPartition2+1):end),...
                                                        n_params2_2,...%was:((idxPartition2+1):end),...
                                                        Ie2_2,...% was:(:,(idxPartition2+1):end),...
                                                        idx_in(idxPartition2+1:end),...
                                                        exitflag2_2,...%was: (1:idxPartition2),...
                                                        [], ... % was: ieprior,...
                                                        [], ... % was: stdprior,...
                                                        Ie1_2(:,end),...
                                                        Directives);
      % and conquer!
      % AICFull = [AICc1_2,AICc2_2];
      polycoefs = cat(2,polycoefs1_2,polycoefs2_2);
      best_order = [best_order1_2,best_order2_2];
      n_params = [n_params1_2,n_params2_2];
      ne = [ne1_2,ne2_2];
      neEnd = [neEnd1_2,neEnd2_2];
      Ie = [Ie1_2,Ie2_2];
      exitflag = [exitflag1_2,exitflag2_2];
      nSteps = [nSteps1_2,nSteps2_2];
      idx_out = [idx_out1_2,idx_out2_2];
      %load recurse_idx.mat r_idx
      %r_idx(4) = r_idx(4) + 1;
      %save('recurse_idx.mat','r_idx')
      %save(['recurse-checks-4-',sprintf('-%04i',r_idx)])
    elseif AICthirds < AICFull & AICthirds <= AIChalves
      % Divide into thirds, same as above but branching into the recursion
      % in the 3 thirds 
      [ne1_3,neEnd1_3,Ie1_3,polycoefs1_3,best_order1_3,n_params1_3,exitflag1_3,nSteps1_3,idx_out1_3] = recurse_AICfit(ne1_3,...
                                                        neEnd1_3,...% was:(:,1:idxPartition2),...
                                                        pp(:,1:idxPartition3),...
                                                        ppstd(:,1:idxPartition3),...
                                                        alpha(:,1:idxPartition3),...
                                                        dt(1:idxPartition3),...
                                                        ne00,...
                                                        A,...
                                                        polycoefs1_3,...%was:(:,:,1:idxPartition2),...
                                                        b_o1_3,...%was:best_order(1:idxPartition2),...
                                                        n_params1_3,...%was(1:idxPartition2),...
                                                        Ie1_3,...% was(:,1:idxPartition2),...
                                                        idx_in(1:idxPartition3),...
                                                        exitflag1_3,...%was: (1:idxPartition2),...
                                                        [], ... % was: ieprior,...
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

      [ne2_3,neEnd2_3,Ie2_3,polycoefs2_3,best_order2_3,n_params2_3,exitflag2_3,nSteps2_3,idx_out2_3] = recurse_AICfit(ne2_3,...
                                                        neEnd2_3,...% was: (:,(idxPartition2+1):end),...
                                                        pp(:,(idxPartition3+1):(2*idxPartition3)),...
                                                        ppstd(:,(idxPartition3+1):(2*idxPartition3)),...
                                                        alpha(:,(idxPartition3+1):(2*idxPartition3)),...
                                                        dt((idxPartition3+1):(2*idxPartition3)),...
                                                        neEnd1_3(:,end),... % 3## FIX!
                                                        A,...
                                                        polycoefs2_3,...%was:(:,:,(idxPartition2+1):end),...
                                                        b_o2_3,...%was:best_order((idxPartition2+1):end),...
                                                        n_params2_3,...%was:((idxPartition2+1):end),...
                                                        Ie2_3,...% was:(:,(idxPartition2+1):end),...
                                                        idx_in((idxPartition3+1):(2*idxPartition3)),...
                                                        exitflag2_3,...%was: (1:idxPartition2),...
                                                        [], ... % was: ieprior,...
                                                        [], ... % was: stdprior,...
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
      [ne3_3,neEnd3_3,Ie3_3,polycoefs3_3,best_order3_3,n_params3_3,exitflag3_3,nSteps3_3,idx_out3_3] = recurse_AICfit(ne3_3,...
                                                        neEnd3_3,...% was: (:,(idxPartition2+1):end),...
                                                        pp(:,(2*idxPartition3+1):end),...
                                                        ppstd(:,(2*idxPartition3+1):end),...
                                                        alpha(:,(2*idxPartition3+1):end),...
                                                        dt((2*idxPartition3+1):end),...
                                                        neEnd2_3(:,end),... % 3## FIX!
                                                        A,...
                                                        polycoefs3_3,...%was:(:,:,(idxPartition2+1):end),...
                                                        b_o3_3,...%was:best_order((idxPartition2+1):end),...
                                                        n_params3_3,...%was:((idxPartition2+1):end),...
                                                        Ie3_3,...% was:(:,(idxPartition2+1):end),...
                                                        idx_in((2*idxPartition3+1):end),...
                                                        exitflag3_3,...%was: (1:idxPartition2),...
                                                        [], ... % was: ieprior,...
                                                        [], ... % was: stdprior,...
                                                        Ie2_3(:,end),...
                                                        Directives);
      % and conquer!
      % AICFull = [AICc1_3,AICc2_3];
      polycoefs = cat(2,polycoefs1_3,polycoefs2_3,polycoefs3_3);
      best_order = [best_order1_3,best_order2_3,best_order3_3];
      n_params = [n_params1_3,n_params2_3,n_params3_3];
      ne = [ne1_3,ne2_3,ne3_3];
      neEnd = [neEnd1_3,neEnd2_3,neEnd3_3];
      Ie = [Ie1_3,Ie2_3,Ie3_3];
      exitflag = [exitflag1_3,exitflag2_3,exitflag3_3];
      nSteps = [nSteps1_3,nSteps2_3,nSteps3_3];
      idx_out = [idx_out1_3,idx_out2_3,idx_out3_3];
      %load recurse_idx.mat r_idx
      %r_idx(4) = r_idx(4) + 1;
      %save('recurse_idx.mat','r_idx')
      %save(['recurse-checks-4-',sprintf('-%04i',r_idx)])
    else
      % AICFull = AICfull(1)*ones(size(dt));
      % Veni Vidi Conquestatum! When current region is sufficiently
      % well modeled and we just return back up
      % ne_End = neEnd;
      idx_out = idx_in;
    end
  end
  % savename = sprintf('ut-%03i-%03i.mat',idx_in(1),numel(idx_in));
  % save(savename,'neEnd')
  % disp(['Leaving recurse_AICfit, nt: ',num2str(numel(idx_out))])

end
