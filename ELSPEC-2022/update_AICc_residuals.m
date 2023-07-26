function AICres = update_AICc_residuals(ne,stdne,ne0,A,alpha,dt,X0,E,dE, ...
                           integtype,IePrior,stdPrior , nmeas,eType,eWidth,varargin)
%
% Calculate corrected Akaike information criterion value for model
% value of the  next electron density profile.
%
% AICres = update_AICc_residuals(ne,stdne,ne0,A,alpha,dt,E0,Y0,gamma,E,integtype)
%
% INPUT:
%  ne         measured electron density profile
%  stdne      standard deviations of ne
%  ne0        modeled electron density profile from the previous time step
%  A          Ion production profile matrix
%  alpha      effective recombination rates
%  dt         time step [s]
%  X0         coefficients of the polynomial
%  E          energies of the dense "model grid"
%  dE         widths of the energy bins
%  integtype  type input for integrate_continuity_equation
%             ('endNe','integrate',or 'equilibrium')
%  IePrior    Apriori flux
%  stdPrior   standard deviation of the apriori
%  eType      type of error distribution: 's', 'n', 't', or 'l'.
%  eWidht     width of transition from normal to long-tailed distribution in units of sigma
%
%
% OUTPUT:
%  AICres     array with scale-and-weighted residuals - for optimising
%             corrected Akaike information criterion value and model
%             selsction
%
%
% IV 2017
% BG 2022
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

n_t = numel(dt);

% the high-resolution spectrum
if nargin < 16 || isempty(varargin{1})
  for it4s = size(X0,1):-1:1
    s(it4s,:) = model_spectrum( X0(it4s,:) , E );
  end
else
  s_type = varargin{1};
  t = cumsum(dt);
  if size(X0,1)> 1
    Xend = X0(1,:).*X0(2,:);
  else
    Xend = X0(1,:);
  end
  if numel(dt) > 1
    for it4s = numel(dt):-1:1
      X_curr = X0(1,:) + (t(it4s)-t(1))/(t(end)-t(1))*(Xend - X0(1,:));
      s(it4s,:) = model_spectrum( X_curr, E, s_type );
    end
  else
    s(1,:) = model_spectrum( X0, E, s_type );  
  end
end
% the updated model ne profile
if strcmp(integtype,'equilibrium')
    nemod = integrate_continuity_equation(ne0(:),s(:),dt(1),A,dE(: ...
                                                      ),alpha(:),integtype);
    nemeas = ne(:,1);
    stdmeas = stdne(:,1);
elseif strcmp(integtype,'equilibriumavg')
  ninteg = length(dt);
  nemod = NaN(length(ne0),ninteg);
  neEnd = NaN(length(ne0),ninteg);
  nemod(:,1) = integrate_continuity_equation(ne0(:),s(1,:)',dt(1),A,dE(:),alpha(:,1),'equilibrium');
  if ninteg > 1
    for tt=2:ninteg
      nemod(:,tt) = integrate_continuity_equation(neEnd(:,tt-1),s(1,:)',dt(tt),A,dE(:),alpha(:,min(end,tt)),'equilibrium');
    end
  end
  nemeas = ne;
  stdmeas = stdne;
else
    switch lower(integtype)
      case 'endne'
        integ_type_end = 'endNe';
      case 'integrate'
        integ_type_end = 'endNe';
      case 'linearend'
        integ_type_end = 'linearend';
      case 'linearint'
        integ_type_end = 'linearend';
      otherwise
        error(['Error. \nUnknown integtype %s in update_AICc' ] , type )
    end

    ninteg = length(dt);
    nemod = NaN(length(ne0),ninteg);
    neEnd = NaN(length(ne0),ninteg);
    if ~exist('s','var')
        keyboard
    end
    nemod(:,1) = integrate_continuity_equation(ne0(:),s(1,:)',dt(1),A,dE(:),alpha(:,1),integtype);
    neEnd(:,1) = integrate_continuity_equation(ne0(:),s(1,:)',dt(1),A,dE(:),alpha(:,1),integ_type_end);
    if ninteg > 1
        for tt=2:ninteg
            nemod(:,tt) = integrate_continuity_equation(neEnd(:,tt-1),s(min(end,tt),:)',dt(tt),A,dE(:),alpha(:,min(end,tt)),integtype);
            neEnd(:,tt) = integrate_continuity_equation(neEnd(:,tt-1),s(min(end,tt),:)',dt(tt),A,dE(:),alpha(:,min(end,tt)),integ_type_end);
        end
    end
    nemeas = ne;
    stdmeas = stdne;
end

% the information criterion value
%AIC = AICc( nemeas(:) , stdmeas(:).^2 , nemod(:) , length(X0) , length(nemeas(:)) );
% AIC = AICc( nemeas(:) , stdmeas(:).^2 , nemod(:) , length(X0) , ...
%            nmeas, eType,eWidth );
%AIC = AICc( nemeas(:) , stdmeas(:).^2 , nemod(:) , numel(X0) , ...
%            nmeas, eType,eWidth );
nParams = size(X0,2).*max(1,size(X0,1)/2);
AICres = AICc_residuals( nemeas(:) , stdmeas(:).^2 , nemod(:) , nParams , ...
                      nmeas, eType,eWidth );
% Ougth to be numel(X0) above: BG-20220630?

% slight regularization for all coefficients
% BG: skipping this because below "reason"
% AIC = [AIC; % BG: Dont understans this term? + sum(X0.^2./1e5.*10.^[1:size(X0,2)],'all') - ...
AICres = [AICres; 
          -1e25*diff(s(end-1:end)).*(diff(s(end-1:end))<0);       % These two terms bias
          -1e25*diff(s([end-10,end])).*(diff(s([end-10,end]))<0)];% agains increasing Ie 

%% the spectrum should go to zero at the high energy end
%AIC = AIC + s(end)^2;

%% reject negative densities, these are possible only with the
%% linear model...
%AIC = AIC + any(nemod(:)<1)*1e10;


% try to damp oscillations by means of smoothing...
% allow larger variations when difference between previous
% solution and the new measured Ne is large..
%if ~isempty(Ieprev)
%    %    AIC = AIC + sum(((Ieprev(:)-s(:)).^2)./1e18);
%    %    AIC = AIC + sum(((Ieprev(:)-s(:)).^2)/mean(9e18*(ne0-ne).^2./stdne.^2));
%    %    AIC = AIC + sum(((Ieprev(:)-s(:)).^2)/(1e15*mean(abs(ne0.^4-ne.^4)./stdne.^4)));
%
%end

% try with a prior model (from an equilibrium solution...)
if ~isempty(IePrior)
    AICres = [AICres; 
           2*(s(:) - IePrior(:))./(stdPrior(:)).^2];
end

% tell the solver that flux at 100 eV is probably not huge...
%s100 = model_spectrum(X0,100);
%AIC = AIC + s100.^2/1e18;


end
