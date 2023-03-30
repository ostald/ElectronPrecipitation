function [AICcSec,polycoefs,best_order,n_params,ne,neEnd,Ie,exitflags] = AICcFitParSeq(pp,ppstd,alpha,dt,ne00,A,polycoefs,ieprior,stdprior,Ie_prev,Directives)
% function [AICcSec,ne,neEnds] = AICcFitParSeq(Out,ne00,idxStart,idxSteps)
% AICCFITPARSEQ - 
% Input parameters
%   Out.A, Out.pp, Out.maxorder, Out.polycoefs, Out.ieprior,
%   Out.ppstd, Out.alpha, Out.dt, Out.Ec, Out.dE, Out.ieprior,
%   Out.ErrType, Out.ErrWidth
% Input parameters that will be split
%   Out.pp, Out.polycoefs, Out.ppstd, Out.alpha, Out.dt,
%   Out.ieprior,
% 
% Input parameters that will be unaffected by 
%   Out.A, Out.maxorder, Out.Ec, Out.dE, Out.ErrType, Out.ErrWidth
% 
% Outputs:
% xOut.polycoefsX, xOut.exitflagX, xOut.outputX, xOut.AICcX,
% xOut.best_orderX, xOut.IeX, xOut.neX, xOut.neEndX

  S_type = Directives.Ietype;
  fms_opts = optimset('fminsearch');
  fms_opts.Display = 'off';%'off';%'final';
  fms_opts.MaxFunEvals = 1e4;
  fms_opts.MaxIter = 1e6;
  fms_opts.TolFun = 1e-8;
  fms_opts.TolX = 1e-8;


  A = Directives.A;
  n_meas = sum(isfinite(pp),'all');
  n_tsteps = numel(dt);

  AICcSec = nan(Directives.maxorder,n_tsteps);
  ne = 0*pp;
  neEnd = 0*pp;
  
  % then run the fit for all different numbers of nodes
  for nn = 1:Directives.maxorder
    % Here we have the parameters for the Generalized
    % Ellison-Ramaty distribution 
    % Ie(E) = I_hat*E^gamma1*exp(-abs((E-E0)/dE)^gamma2):
    %             log(I_hat) dE  E_0  gamma_1, gamma_2
    if strcmp(lower(S_type),'ger')
      % default parameters for the GER-type spectra
      X_defaults = [10,        4,  0.5, 1,       1];
      X_minall   = [-1         0.3 0   -2        0.5];
      X_MAXall   = [35        30  30    3        4];
      X0 = X_defaults(1:nn+1);
      X0min = X_minall(1:nn+1);
      X0max = X_MAXall(1:nn+1);
    else
      % defaults for the exponential of polynomials
      if nn==1
        % X0 = repmat([20,-5],2,1); % This seemed very peculiar! BG-20220927
        X0 = repmat([20,-5],1,1); 
      else
        % X0 = [polycoefs(nn-1,1:nn,1),-0.5];
        X0 = polycoefs{nn-1,1}(1,:);
        X0(:,end+1) = -0.5;
      end
      X0min = -Inf(size(X0));
      X0max = Inf(size(X0));
      X0max(:,end) = 0; % the coefficient to the largest power of E has
      
    end
    if Directives.InterpSpec
      X0(2,:) = 1;
      X0min(2,:) = 0.8;
      X0max(2,:) = 1.2;
    end
    % Straight interpolation of parameters doesnt seem to work all
    % that well. Tested on the GER-distribution, next test for
    % polynomials.
    % X0 = repmat(X_defaults(1:nn+1),Directives.InterpSpec+1,1);
    % X0min = repmat(X_minall(1:nn+1),Directives.InterpSpec+1,1);
    % X0max = repmat(X_MAXall(1:nn+1),Directives.InterpSpec+1,1);

    % If that too stumbles then we have to go for something like a
    % 10-20 % variation of parameters during a period.
    % if Directives.InterpSpec
    % X0 = X_defaults(1:nn+1));
    % X0(2,:) = 1;
    % X0min = X_minall(1:nn+1);
    % X0max = X_MAXall(1:nn+1);
    % X0min(2,:) = 0.8;
    % X0max(2,:) = 1.2;
    % end
    
% $$$     if isnan(polycoefs(nn,1,1,1))
% $$$       if nn==1
% $$$         X0 = repmat([20,-5],2,1);
% $$$       else
% $$$         % X0 = [polycoefs(nn-1,1:nn,1),-0.5];
% $$$         X0 = squeeze(polycoefs(nn-1,1:nn,1,:))';
% $$$         X0(:,end+1) = -0.5;
% $$$       end
% $$$     else % nn > 1 
% $$$          % then we take the parameters from the look-ahead fittings 
% $$$       X0 = polycoefs(nn,1:(nn+1),1);
% $$$     end
    
    % do not limit the coefficients, this is handled with a
    % prior now
    
    if ~isempty(ieprior)
      [x,fval,exitflag] = fminsearchbnd( @(p) ElSpec_fitfun(p, ...
                                                        pp, ...
                                                        ppstd, ...
                                                        ne00, ...
                                                        A, ...
                                                        alpha, ...
                                                        dt, ...
                                                        Directives.Ec, ...
                                                        Directives.dE , ...
                                                        'integrate', ...
                                                        ieprior, ...
                                                        stdprior, ...
                                                        n_meas, ...
                                                        Directives.ErrType, ...
                                                        Directives.ErrWidth,...
                                                        S_type), ...
                                         X0, X0min, X0max, fms_opts);
    else
      % the normal fit
      [x,fval,exitflag] = fminsearchbnd( @(p) ElSpec_fitfun(p, ...
                                                        pp, ...
                                                        ppstd, ...
                                                        ne00, ...
                                                        A, ...
                                                        alpha, ...
                                                        dt, ...
                                                        Directives.Ec, ...
                                                        Directives.dE, ...
                                                        'integrate', ...
                                                        [], ...
                                                        [], ...
                                                        n_meas, ...
                                                        Directives.ErrType, ...
                                                        Directives.ErrWidth,...
                                                        S_type), ...
                                         X0, X0min, X0max, fms_opts);
    end
    
    % collect results
    for i_t = 1:n_tsteps
      polycoefs{nn,i_t} = x;
    end
    exitflags(nn,1:n_tsteps) = exitflag;
    AICcSec(nn,1:n_tsteps) = fval;
    
  end

  % pick the best AICc value
  [minAIC,locminAIC] = min(AICcSec(:,1));
  n_params(1:n_tsteps) = locminAIC+1;
  best_order(1:n_tsteps) = locminAIC;
  % In med samma sak som i updateAICc.m
  curr_pars = polycoefs{locminAIC,1};
  if min(size(curr_pars)) < 2
    for i1 = n_tsteps:-1:1
      Ie(:,i1) = model_spectrum(curr_pars, ...
                                Directives.Ec, S_type);
    end
  else
    dbstop
    X0 = curr_pars(1,:);
    Xend = X0.*curr_pars(2,:);
    t_curr = cumsum(dt);
    for it4s = numel(dt):-1:1
      now_pars = ( X0 + ...
                   ( t_curr(it4s) - t_curr(1))/(t_curr(end) - t_curr(1))*... ...
                   ( Xend - X0 ) );
      Ie(:,it4s) = model_spectrum( now_pars, Directives.Ec, S_type );
    end

  end
  ne0 = ne00;
  for i1 = 1:n_tsteps
    
    % the modeled Ne (which should match with the measurements)
    ne(:,i1) = integrate_continuity_equation( ne0, ...
                                              Ie(:,i1), ...
                                              dt(i1), ...
                                              Directives.A, ...
                                              Directives.dE(:), ...
                                              alpha(:,i1), ...
                                              'integrate' ...
                                              );
    % modeled Ne at end of the time step
    neEnd(:,i1) = integrate_continuity_equation( ne0, ... 
                                                 Ie(:,i1), ...
                                                 dt(i1),...
                                                 Directives.A, ...
                                                 Directives.dE(:), ...
                                                 alpha(:,i1), ...
                                                 'endNe' ...
                                                 );
    ne0 = neEnd(:,i1);
  end

end
