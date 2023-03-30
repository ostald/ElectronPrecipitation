function ElSpecOut = ElSpec_iqt_ic(varargin)
% ElSpec Estimate primary energy spectra of precipitating electrons from
% EISCAT incoherent scatter data.
%
%
% The primary energy spectra are modeled  using a parametrix model
% of the form
% s = E.*exp(P_n(E))
% where P_n is an n'th order polynomial.
%
% Corrected Akaike information criterion is used for seleting the
% optimal order of the polynomial P_n.
%
% ElSpecOut = ElSpec(varargin)
%
%
% INPUT:
%
%  name-value pairs, the available arguments are listed below.
%
%  At least  ppdir OR fitdir must be defined.
%
%  For example ElSpec('ppdir',<path-to-pp-files>)
%
%  ppdir        Path to GUISDAP raw densities. See details.
%  fitdir       Path to GUISDAP parameter fit results. See details.
%  experiment   EISCAT experiment name, e.g. 'beata' or 'arc1'
%  radar        Radar, 'uhf', 'vhf', or 'esr'
%  version      Experiment version number.
%  hmin         Minimum altitude in km, must be between 80 and 150 km.
%  hmax         Maximum altitude in km, must be between 80 and 150 km.
%  btime        Analysis start time [year,month,day,hour,minute,second]
%  etime        Analysis end time [year,month,day,hour,minute,second]
%  ionomodel    Ion production model, 'Fang' or 'Sergienko'. See details.
%  recombmodel  Recombination model, 'Rees', 'ReesO2+', 'ReesN2+', 'ReesNO+',
%               'SheehanGrO2+', 'SheehanGrN2+', 'SheehanGrNO+',
%               'delPozo1', or 'delPozo2'. See details.
%  integtype   Type of continuity function integration, 'endNe',
%              'integrate', or 'equilibrium'. See details.
%  egrid        Energy grid [eV]
%  maxorder     range of the number of nodes in the spline
%               expansion.
%  plotres      plot results with ElSpecPlot? 1=plot, 0=no plot,
%               default 1
%  tres         "type" of time resolution 'best' or 'dump' [default]
%  emin         smallest energy to include in FAC and power
%               estimation, default 1 keV
%  ieprior      prior model test....,requires a lot of manual
%                work... use empty array to skip the prior...
%  stdprior
%  ninteg       number of ne-profiles to use in each integration
%               period, default 6
%  nstep        number of ne-slices in each time-step, default 1
%  saveiecov    logical, should the large covariance matrices of
%               the flux estimates be saved? default false.
%  alpha_eff    effective recombination rates calculated from IC
%               integration
%  iri_ic       species densities at times ts calculated with IC integration    
%
% OUTPUT:
%  ElSpecOut    A MATLAB structure with fields:...
%
%
%
% Details:
%
%  The analysis uses two sets of GUISDAP results as input: a set of
%  high-time resolution raw electron densities, and a set of lower
%  resolution 'normal' fit results. These cannot be produced in a
%  single GUISDAP run, because time resolution of the raw densities
%  would then be matched to time resolution of the full fit. The
%  raw densities in ppdir can be produced by means of setting
%
%  analysis_ppshortlags=1
%
%  and
%
%  analysis_pponly=1
%
%  in GUISDAP. The
%  GUISDAP output will be experiment- and
%  radar-specific. Currently, we have a readers for uhf beata and
%  arc1 only.
%
%  Two models for ion production profiles by monoenergetic
%  electrons are implemented. These are Sergienko1993 and
%  Fang2010. See
%  help ion_production_Sergienko1993
%  and
%  help ion_production_Fang2010
%  for details.
%
%  Eight recombination models are implemented. 'Rees' is based on
%  ion abundances from IRI and values from Rees
%  (1989). 'SheehanGr' uses values for ground state ions from
%  Sheehan and St.-Maurice (2004). 'SheehanEx' uses values for
%  vibrationally excited ions from Sheehan and St.-Maurice (2004).
%  'ReesO2+','ReesN2+', and 'ReesNO+' are Rees recombination rates for
%  pure O2+, N2+, adn NO+, correspondingly. 'SheehanGrO2+',
%  'SheehanGrN2+' and 'SheehanGrNO+' are the corresponding Sheehan
%  and St.-Maurice (2004) recombination rates for vibrational
%  ground-states. 'delPozo1' and 'delPozo2' are simple analytic
%  approximations. See
%  help effective_recombination_rate
%  for details.
%
%  Three different methods for integrating the electron continuity
%  equation are available. 'endNe' ignores the actual integration,
%  and returns the electron density profile at end of the time
%  step. 'integrate' returns average of the time-dependent density
%  profile over the time-step dt. 'equilibrium' ignores the initial
%  condition and assumes a photochemical equilibrium. 'equilibrium'
%  is always forced for the very first time step, where the initial
%  condition is inknown. See
%  help integrate_continuity_equation
%
%
%  The software uses tabulated values of IRI parameters and uses
%  the MATLAB Aerospace toolbox implementation of MSIS model. The
%  IRI tables are provided as part of the distribution,
%  and are automatically read by the function 'modelParams'. An R
%  function for creating new tables, createModelParameterFile, is
%  provided within the MATLAB distribution. The R function uses
%  packages IRI2016 and R.matlab. The latter is a standard CRAN
%  pacakge, whereas the package IRI2016 can be requested from Ilkka
%  Virtanen (ilkka.i.virtanen@oulu.fi).
%
% IV 2017, 2018
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

persistent ElSpecLicenseShown


% display the copyright message only once
if isempty(ElSpecLicenseShown)
    ElSpecLicenseShown = true;
    disp('  ')
    disp('Copyright 2015-2022, Ilkka Virtanen, Bjorn Gustavsson and Juha Vierinen')
    disp(['This is free software, licensed under GNU GPL version 2 or later.'])
    disp('  ')
    disp(['This program is distributed in the hope that it will be ' ...
          'useful, '])
    disp(['but WITHOUT ANY WARRANTY; without even the implied ' ...
          'warranty of '])
    disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.')
    disp('See the GNU General Public License for more details')
    pause(1)
end


% power profiles
defaultPPdir = [];
checkPPdir = @(x) (exist(x,'dir')|isempty(x)|strcmp(x,'simu'));

% plasma parameters
defaultFitdir = [];
checkFitdir = @(x) (exist(x,'dir')|isempty(x));

% name of the experiment
defaultExperiment = 'beata';
checkExperiment = @(x) (isstr(x));

% radar, only uhf exps implemented at the moment
defaultRadar = 'uhf';
validRadar = {'uhf','vhf','esr','42m','32m','esr42m','esr32m'};
checkRadar = @(x) any(validatestring(x,validRadar));

% experiment version
defaultVersion = 1;
checkVersion = @(x) (isnumeric(x) & length(x)==1);

% minimum altitude
defaultHmin = 80;
checkHmin = @(x) (x>=80 & x <150);

% maximum altitude
defaultHmax = 150;
checkHmax = @(x) ( x>80 & x<=150);

% ion production model
defaultIonomodel = 'Fang';
validIonomodel = {'Fang','Sergienko'};
checkIonomodel = @(x) any(validatestring(x,validIonomodel));

% recombination model
defaultRecombmodel = 'SheehanGr';
validRecombmodel = {'SheehanGr','SheehanEx','Rees','ReesO2+', ...
                    'ReesN2+','ReesNO+','SheehanGrO2+', ...
                    'SheehanGrN2+','SheehanGrNO+','delPozo1','delPozo2'};
checkRecombmodel = @(x) any(validatestring(x,validRecombmodel));

% type of integration
defaultIntegtype = 'integrate';
validIntegtype = {'integrate','equilibrium','endNe'};
checkIntegtype = @(x) any(validatestring(x,validIntegtype));

% energy grid
defaultE = logspace(1,5,300);
checkE = @(x) (isnumeric(x) & length(x)>1 & all(x>0));

% maximum order of the polynomial model
defaultMaxorder = 5;
checkMaxorder = @(x) (length(x)==1 & all(x>=1));

% plot the results?
defaultPlotres = false;
checkPlotres = @(x) ((isnumeric(x)|islogical(x))&length(x)==1);

% time resolution (best available or EISCAT dump length)
defaultTres = 'dump';
validTres = {'dump','best'};
checkTres = @(x) any(validatestring(x,validTres));

% minimum energy for FAC and power calculation
defaultEmin = 1e3;
checkEmin = @(x) (isnumeric(x) & length(x)==1 & x>0);

% a priori differential number flux
defaultIeprior = [];
checkIeprior = @(x) (ismatrix(x));

% a priori standard deviations
defaultStdprior = [];
checkStdprior = @(x) (ismatrix(x));

% number of time-slices integrated together
defaultNinteg = 6; % 30 s in beata, note that only the first is used
                   % with full weight. the time resolution will be
                   % effectively one slice for high Ne
checkNinteg = @(x) (isnumeric(x) & length(x)==1 & all(x>0));

% number of time-slices to proceed on each time-step
defaultNstep = 1;
checkNstep = @(x) (isnumeric(x) & length(x)==1 & all(x>0));

% save the  large flux covariance matrix?
defaultSaveiecov = false;
checkSaveiecov = @(x) ((isnumeric(x)|islogial(x))&length(x)==1);

% standard deviation for model Ne in the lowest observed altitude when the guisdap fit has failed
defaultBottomStdFail = 1e10;
checkBottomStdFail = @(x) (isnumeric(x)&length(x)==1 & x>0);


% start time
defaultBtime = [];

% end time
defaultEtime = [];

% function for checking the time limits
checkTlim = @(x) (isnumeric(x) & length(x)==6 & all(x>=0));

% default time-direction
% defaultTdir = 1; 
% checkTdir = @(x) (isnumeric(x) & length(x)==1 & ((x==1) || (x==-1)));

% Type of neg-log-likelihood-function
% For normal-distributed errors standard least-squares fitting is
% appropriate, for this choise use 's' or 'n'. For more
% heavy-tailed noise-distributions use a Lorentzian ('l') or tanh ('t')
% weighting function with a width-parameter specifying at how many
% sigma the error-distribution starts to significantly deviate from
% a Gaussian (typically 3-4?)
defaultErrType = 's';
validErrType = {'s','n','t','l'}; % for standard/normal, tanh-weighting,
                                  % log-weighting (Cauchy-distributed noise)
checktErrType = @(x) any(validatestring(x,validErrType));

% sigma-width of transition to long-tailed residfual distribution
% for the robust fitting AICc
defaultErrWidth = 3;
checktErrWidth = @(x) (length(x)==1 & all(x>=1));

% Outlier-settings
defaultOutliers = ones(1,5);
checkOutliers = @(x) (all(isnumeric(x)) & size(x,2)==5 & size(x,1)>=1 & all(x(:)>0));

% Outfilename-settings
defaultOutfilename = '';
checkOutfilename = @(x) isstring(x)| ischar(x);

% time resolution (best available or EISCAT dump length)
% type of electron-spectra parameterization
defaultIetype = 'p';
validIetype = {'p','ger'};
checkIetype = @(x) any(validatestring(x,validIetype));

% Do linear interpolation of electron-spectra
% 
defaultInterpSpec = 0;
checkInterpSpec = @(x) ((numel(x) == 1) & isnumeric(x) & ...
                        ( (x == 0) || (x == 1) ) );

% use only field-aligned data?
defaultFAdev = 3;
checkFAdev = @(x) (isnumeric(x)&length(x)==1);
                    

% custom densities
defaultiri_ic = 0;
checkiri_ic = @(x) 1; %(ismatrix(x));

% custom recombintaion rates
defaultalpha_eff = 0;
checkalpha_eff = @(x) (ismatrix(x) || 0);


if exist('inputParser') %#ok<EXIST> 
  % parse the input
  p = inputParser;
  
  addParameter(p,'btime',defaultBtime,checkTlim);
  addParameter(p,'etime',defaultEtime,checkTlim);
  addParameter(p,'ppdir',defaultPPdir,checkPPdir);
  addParameter(p,'fitdir',defaultFitdir,checkFitdir);
  addParameter(p,'experiment',defaultExperiment,checkExperiment);
  addParameter(p,'radar',defaultRadar,checkRadar);
  addParameter(p,'version',defaultVersion,checkVersion);
  addParameter(p,'hmin',defaultHmin,checkHmin);
  addParameter(p,'hmax',defaultHmax,checkHmax);
  addParameter(p,'ionomodel',defaultIonomodel,checkIonomodel);
  addParameter(p,'recombmodel',defaultRecombmodel,checkRecombmodel);
  addParameter(p,'integtype',defaultIntegtype,checkIntegtype);
  addParameter(p,'egrid',defaultE,checkE);
  addParameter(p,'maxorder',defaultMaxorder,checkMaxorder);
  addParameter(p,'plotres',defaultPlotres,checkPlotres);
  addParameter(p,'tres',defaultTres,checkTres);
  addParameter(p,'emin',defaultEmin,checkEmin);
  addParameter(p,'ieprior',defaultIeprior,checkIeprior);
  addParameter(p,'stdprior',defaultStdprior,checkStdprior);
  addParameter(p,'ninteg',defaultNinteg,checkNinteg);
  addParameter(p,'nstep',defaultNstep,checkNstep);
  addParameter(p,'saveiecov',defaultSaveiecov,checkSaveiecov);
%   addParameter(p,'Tdir',defaultTdir,checkTdir);
  addParameter(p,'fadev',defaultFAdev,checkFAdev);
  addParameter(p,'bottomstdfail',defaultBottomStdFail,checkBottomStdFail);
  addParameter(p,'ErrType',defaultErrType,checktErrType);
  addParameter(p,'ErrWidth',defaultErrWidth,checktErrWidth);
  addParameter(p,'Outliers',defaultOutliers,checkOutliers);
  addParameter(p,'Outfilename',defaultOutfilename,checkOutfilename);
  addParameter(p,'Ietype',defaultIetype,checkIetype);
  addParameter(p,'InterpSpec',defaultInterpSpec,checkInterpSpec);

  addParameter(p,'alpha_eff',defaultalpha_eff,checkalpha_eff);
  addParameter(p,'iri_ic',defaultiri_ic,checkiri_ic);

  parse(p,varargin{:})
  
  %out = struct;
  
  % collect all the inputs to the output list
  out = p.Results;
  out.E = out.egrid;
  
else
  
  def_pars.btime = defaultBtime;
  def_pars.etime = defaultEtime;
  def_pars.ppdir = defaultPPdir;
  def_pars.fitdir = defaultFitdir;
  def_pars.experiment = defaultExperiment;
  def_pars.radar = defaultRadar;
  def_pars.version = defaultVersion;
  def_pars.hmin = defaultHmin;
  def_pars.hmax = defaultHmax;
  def_pars.ionomodel = defaultIonomodel;
  def_pars.recombmodel = defaultRecombmodel;
  def_pars.integtype = defaultIntegtype;
  def_pars.egrid = defaultE;
  def_pars.maxorder = defaultMaxorder;
  def_pars.plotres = defaultPlotres;
  def_pars.tres = defaultTres;
  def_pars.emin = defaultEmin;
  def_pars.ieprior = defaultIeprior;
  def_pars.stdprior = defaultStdprior;
  def_pars.ninteg = defaultNinteg;
  def_pars.nstep = defaultNstep;
  def_pars.saveiecov = defaultSaveiecov;
  
  %def_pars.Tdir = defaultTdir;
  def_pars.ErrType = defaultErrType;
  def_pars.ErrWidth = defaultErrWidth;
  def_pars.Outliers = defaultOutliers;
  def_pars.Outfilename = defaultOutfilename;
  def_pars.Ietype = defaultIetype;
  def_pars.InterpSpec = defaultInterpSpec;
  
  def_pars.iri_ic = defaultiri_ic;
  def_pars.alpha_eff = defaultalpha_eff;


  out = parse_pv_pairs(def_pars,varargin);
  out.E = out.egrid;
  
end
% check that we have either ppdir or fitdir, or both
if isempty(out.ppdir) & isempty(out.fitdir)
    error('Either ppdir or fitdir must be given')
end

% we have input routines only for uhf arc1 and beata power profiles
if ~isempty(out.ppdir)
  wrongexp = false;
  if strcmp(out.radar,'uhf')
    if ~any(strcmp(out.experiment,{'beata','arc1'}))
      wrongexp = true;
    end
  else
    wrongexp = true;
  end
  if wrongexp
    error(['Power profile input routines are available only for ' ...
           'UHF arc1 and UHF beata experiments.'])
  end
end

% check that nstep is not larger than ninteg
% $$$ if out.nstep > out.ninteg
% $$$     error('nstep must not be larger than ninteg.')
% $$$ end

% Then print all the inputs
fprintf('\nElSpec input arguments\n\n')

fprintf('%12s: %s\n','ppdir',out.ppdir);
fprintf('%12s: %s\n','fitdir',out.fitdir);
fprintf('%12s: %s\n','experiment',out.experiment);
fprintf('%12s: %s\n','radar',out.radar);
fprintf('%12s: %i\n','version',out.version);
fprintf('%12s: %5.1f\n','hmin [km]',out.hmin);
fprintf('%12s: %5.1f\n','hmax [km]',out.hmax);
fprintf('%12s: ','btime');
for kk=1:length(out.btime)
    fprintf('%i ',out.btime(kk));
end
fprintf('\n');
fprintf('%12s: ','etime');
for kk=1:length(out.etime)
    fprintf('%i ',out.etime(kk));
end
fprintf('\n');
fprintf('%12s: %s\n','ionomodel',out.ionomodel);
fprintf('%12s: %s\n','recombmodel',out.recombmodel);
fprintf('%12s: %s\n','integtype',out.integtype);
% $$$ fprintf('%12s: ','E [keV]');
% $$$ for kk=1:(length(out.egrid)/10)
% $$$     if kk>1
% $$$         fprintf('%14s','');
% $$$     end
% $$$     fprintf(' %5.2f',out.E((kk-1)*10+1:min(length(out.E),kk*10))/1000);
% $$$     fprintf('\n');
% $$$ end
fprintf('%12s: %i\n','maxorder',out.maxorder);
fprintf('%12s: %i\n','plotres',out.plotres);
fprintf('%12s: %s\n','tres',out.tres);
fprintf('%12s: %10.0f\n','Emin [eV]',out.emin);
fprintf('%12s: %i\n','ninteg',out.ninteg);
% fprintf('%12s: %i\n','nstep',out.nstep);




% read the data and model values
% a hack to allow reading in simulated data
if strcmp(out.ppdir,'simu')
    load(out.fitdir);
    out.h = simudata.h;
    out.ts = simudata.ts;
    out.te = simudata.te;
    out.pp = simudata.pp;
    out.ppstd = simudata.ppstd;
    out.par = simudata.par;
    out.parstd = simudata.parstd;
    out.iri = simudata.iri;
else
    readIRI = false;
    if any(strcmp(out.recombmodel,{'SheehanGr','SheehanEx','Rees'}))
        readIRI = true;
    end
    [out.h,out.ts,out.te,out.pp,out.ppstd,out.par,out.parstd,out.iri,out.f107,out.f107a,out.f107p,out.ap,out.loc,out.azel,out.I] = ...
        readFitData( out.ppdir , out.fitdir , out.hmin , out.hmax , ...
                     out.btime , out.etime , out.experiment , out.radar , ...
                     out.version , out.tres , readIRI, p.Results.fadev , p.Results.bottomstdfail);
    if strcmp(out.recombmodel,'SheehanGrFlipchem')
    out.iri = calculateFlipchemComposition(out.ts,out.h,out.par,out.pp,out.loc,out.iri);
    end
    % warn about the ESR compositions
    if strcmp(p.Results.radar,'esr')
        if readIRI
            disp(['WARNING: With the selected recombination model ElSpec calculates the ion  compositions ' ...
                  'with EISCAT Tromso coordinates, the compositions may be ' ...
                  'incorrect for the ESR radar'])
            % disp(['WARNING: ElSpec calculates the ion and neutral compositions ' ...
            %       'with EISCAT Tromso coordinates, the compositions may be ' ...
            %       'incorrect for the ESR radar'])
        end
    end

    if numel(out.iri_ic) > 1
        %replace density model with IC values
        out.iri = out.iri_ic;
    end

end
% nt = 2^floor(log2(numel(out.ts)));
nt = numel(out.ts); % 128*floor(numel(out.ts)/128); 
% Change above since there is no FFT-type reason to stick to powers
% of 2 for the time. We should be just as fine integrating over
% odd-number of time-steps if we modify the partitioning to be up
% to floor(end/2) or ceil(end/2) and then the rest. /BG 20220926 
out.ts = out.ts(1:nt);
out.te = out.te(1:nt);
out.pp = out.pp(:,1:nt);
out.par = out.par(:,:,1:nt);
out.ppstd = out.ppstd(:,1:nt);
out.parstd = out.parstd(:,:,1:nt);
out.iri = out.iri(:,:,1:nt);

% $$$ idx_Out2 = [91, 123, 421];
% $$$ idx_Out1 = [12];
% $$$ idx_Outdt = [1,2];
% $$$ idx_Outdz = [5,10,2];
% $$$ idx_OutCne = [3 3 10];

if any(isfinite(out.Outliers))
  for i1 = 1:size(out.Outliers,1)
    i1_1 = out.Outliers(i1,1);
    i1_2 = out.Outliers(i1,2);
    i2_1 = out.Outliers(i1,3);
    i2_2 = out.Outliers(i1,4);
    C = out.Outliers(i1,5);
    out.pp(i1_1:i1_2,i2_1:i2_2) = out.pp(i1_1:i1_2,i2_1:i2_2)*C;
  end
end

% disp('************************************')
% disp('TESTING with SIC composition, ElSpec, line 318!')
% disp('************************************')
% %load('/Users/ilkkavir/Documents/EISCAT_LAB/ElectronPrecipitationSpectra/SIC_composition/Ilkalle_composition.mat');
% load('Ilkalle_composition.mat');
% tsic = posixtime(datetime(t,'convertfrom','datenum'));
% %interpolate to the correct height grid
% O2p = zeros(length(out.h),length(tsic));
% Op = zeros(length(out.h),length(tsic));
% NOp = zeros(length(out.h),length(tsic));
% for indd = 1:length(tsic)
%    O2p(:,indd) = interp1(Hgt,O2plus(:,indd),out.h,'linear','extrap');
%    Op(:,indd) = interp1(Hgt,Oplus(:,indd),out.h,'linear','extrap');
%    NOp(:,indd) = interp1(Hgt,NOplus(:,indd),out.h,'linear','extrap');
% end
% %interplation in time
% for indd = 1:length(out.h)
%    out.iri(indd,8,:) = interp1(tsic,NOp(indd,:),out.ts,'linear','extrap');
%    out.iri(indd,9,:) = interp1(tsic,O2p(indd,:),out.ts,'linear','extrap');
%    out.iri(indd,10,:) = interp1(tsic,Op(indd,:),out.ts,'linear','extrap');
% end


% time step sizes
out.dt = diff( out.te );
out.dt = [out.dt(1);out.dt(:)];

% some dimensions and initializations
nt = length(out.ts); % number of time steps
nh = length(out.h);  % number of heights
nE = length(out.E) - 1; % number of energy bins
out.ne = NaN(nh,nt);   % an array for electron density estimates
out.neEnd = NaN(nh,nt);   % an array for electron density estimates
                          % at integration end points
out.neEndStd = NaN(nh,nt);   % an array for electron density estimates
out.neEndCov = NaN(nh,nh,nt);   % an array for electron density estimates

out.Ie = NaN(nE,nt);   % an array for flux estimates
out.AICc = NaN(out.maxorder,nt); % an array for the information
                                  % criterion values

out.IeCov = NaN(nE,nE,nt);   % Covariances of flux estimates
out.IeStd = NaN(nE,nt); % standard deviation of flux estimates
%out.alpha = NaN(nh,nt); % an array for the effective recombination
%                        % rates

out.q = NaN(nh,nt); % an array for ion production rates
out.polycoefs = cell(out.maxorder,nt);
out.polycoefsCovar = NaN(out.maxorder,out.maxorder+1,out.maxorder+1,nt);
out.exitflag = NaN(out.maxorder,nt);
out.output = cell(out.maxorder,nt);
out.best_order = NaN(1,nt);
out.FAC = NaN(1,nt);
out.FACstd = NaN(1,nt);
out.Pe = NaN(1,nt);
out.PeStd = NaN(1,nt);

% options for the MATLAB fit routines
fms_opts = optimset('fminsearch');
fms_opts.Display = 'off';%'off';%'final';
fms_opts.MaxFunEvals=1e4;
fms_opts.MaxIter=1e6;
fms_opts.TolFun=1e-8;
fms_opts.TolX=1e-8;

% name of the output file
% outfilename = ['ElSpec_',datestr(datetime(round(out.ts(1)), ...
%                                            'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'-',datestr(datetime(round(out.te(end)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'_',out.experiment,'_',out.radar,'_',out.ionomodel,'_',out.recombmodel,'_',out.integtype,'_',num2str(out.ninteg),'_',num2str(out.nstep),'_',out.tres,'_',datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'];

if ~isfield(out,'Outfilename') ||  isempty(out.Outfilename)
  outfilename = ['ElSpec_',...
                 datestr(datevec(unixtime2mat(out.ts(1))),'yyyymmddTHHMMSS'),'-',...
                 datestr(datevec(unixtime2mat(out.ts(end))),'yyyymmddTHHMMSS'),'_',...
                 out.experiment,'_',...
                 out.radar,'_',...
                 out.ionomodel,'_',...
                 out.recombmodel,'_',...
                 out.Ietype,'_',...
                 out.integtype,'_',...
                 num2str(out.ninteg),'_'...
                 out.tres,'_',...
                 datestr(now,'yyyymmddTHHMMSS'),'.mat'];
else
  outfilename = out.Outfilename;
end
out.outfilename = outfilename;
fprintf('Will write results in:\n %s\n',outfilename);



% Update the ion production matrix
[A,Ec,dE] = ion_production(out.E, ...
                           out.h*1000, ...
                           out.iri(:,4,1), ...
                           out.iri(:,5,1), ...
                           out.iri(:,6,1), ...
                           out.iri(:,7,1), ...
                           out.iri(:,1,1),...
                           out.ionomodel, out.I);
A(isnan(A)) = 0;


% update the effective recombination coefficient. Assume N2+
% density to be zero
if numel(out.alpha_eff) > 1
    %replace standart model with IC values
    out.alpha = out.alpha_eff; 
else
    for it = 1:numel(out.ts)
      out.alpha(:,it) = effective_recombination_coefficient(out.h, ...
                                                    out.par(:,3,it), ...
                                                    out.iri(:,9,it), ...
                                                    out.iri(:,9,it).*0, ...
                                                    out.iri(:,8,it), ...
                                                    out.recombmodel );
    end
end

% save interval is 100 step, independently from the time resolution
ndtsave = 20;%ceil(mean(120./out.dt));
breaknow = 0;
% iteration over all time steps

% Directives.A
% Directives.maxorder
% Directives.Ec
% Directives.dE
% Directives.ErrType
% Directives.ErrWidth
Directives.A = A;
Directives.maxorder = out.maxorder;
Directives.Ec = Ec; 
Directives.dE = dE;
Directives.ErrType = out.ErrType;
Directives.ErrWidth = out.ErrWidth;
Directives.ninteg = out.ninteg;
Directives.Ietype = out.Ietype;
Directives.InterpSpec = out.InterpSpec;
out.A = A;
out.Ec = Ec;
out.dE = dE;

pp = out.pp;
ppstd = out.ppstd;
alpha = out.alpha;
dt = out.dt;
polycoefs = out.polycoefs;
ieprior = out.ieprior;
stdprior = out.stdprior;
% 1st fit spectrum for ne0 from the first few electron-density-profiles
[AICc1,polycoefs1,best_order1,n_params1,ne1,neEnd1,Ie1] = AICcFitParSeq(pp(:,1:Directives.ninteg),...
                                                  ppstd(:,1:Directives.ninteg),...
                                                  alpha(:,1:Directives.ninteg),...
                                                  30*dt(1:Directives.ninteg),...
                                                  mean(pp(:,1:Directives.ninteg),2,'omitnan'),...
                                                  A,...
                                                  polycoefs(:,1:Directives.ninteg),...
                                                  [],...
                                                  [],...
                                                  [],... % Ie(t) constant
                                                  Directives);
ne0 = sqrt(A*(Ie1(:,1).*out.dE')./alpha(:,1));% neEnd(:,end);
Ie0 = Ie1;
% $$$ [AICc1,polycoefs1,best_order1,n_params1,ne1,neEnd1,Ie1] = AICcFitParSeq(pp(:,1:Directives.ninteg),...
% $$$                                      ppstd(:,1:Directives.ninteg),...
% $$$                                      alpha(:,1:Directives.ninteg),...
% $$$                                      dt(1:Directives.ninteg),...
% $$$                                      ne00,...
% $$$                                      A,...
% $$$                                      polycoefs(:,:,1:Directives.ninteg),...
% $$$                                      [],...
% $$$                                      [],...
% $$$                                      Directives);
% $$$ ne0 = sqrt(A*(Ie1(:,1).*out.dE')./alpha(:,1));% neEnd(:,end);
%### Magical number: 128 - TODO make adjustable
iStart = 1;
iIntervall = 128;
while iStart < numel(dt)
  
%   while iStart+iIntervall-1 > numel(dt) & iIntervall > 1
%     iIntervall = floor(iIntervall/2); % Attempt at handling
%                                       % not-powers of 2 BG-20220926
%   end
%   if iIntervall == 1 || iStart+iIntervall-1 > numel(dt)
%     break
%   end
  iEnd = iStart+iIntervall-1;
  if iEnd > numel(dt) || abs(iEnd-numel(dt)) < 16
      iEnd = numel(dt);
  end
  disp([datestr(now,'HH:MM:SS'),' ',num2str([iStart iIntervall iEnd])])
  iHalfway = floor((iStart+iEnd)/2);
  [A,Ec,dE] = ion_production(out.E, ...
                           out.h*1000, ...
                           out.iri(:,4,iHalfway), ...
                           out.iri(:,5,iHalfway), ...
                           out.iri(:,6,iHalfway), ...
                           out.iri(:,7,iHalfway), ...
                           out.iri(:,1,iHalfway),...
                           out.ionomodel, out.I);

  A(isnan(A)) = 0;


  % 2nd fit spectrum for the entire event
  [cAICc,cpolycoefs,cbest_order,cn_params,cne,cneEnd,cIe,cexitflag] = AICcFitParSeq(pp(:,iStart:iEnd),...
                                                    ppstd(:,iStart:iEnd),...
                                                    alpha(:,iStart:iEnd),...
                                                    dt(iStart:iEnd),...
                                                    ne0,...
                                                    A,...
                                                    polycoefs(:,iStart:iEnd),...
                                                    ieprior,...
                                                    stdprior,...
                                                    Ie0(:,end),...
                                                    out);
  % 3rd enter the recursion
  % disp(['Calling recurse_AICfit, nt: ',num2str(numel(dt(iStart:iEnd)))])
  [cne,cneEnd,cIe,cpolycoefs,cbest_order,cn_params,cexitflag,cnSteps] = recurse_AICfit(cne,...
                                                    cneEnd,...
                                                    pp(:,iStart:iEnd),...
                                                    ppstd(:,iStart:iEnd),...
                                                    alpha(:,iStart:iEnd),...
                                                    dt(iStart:iEnd),...
                                                    ne0,...
                                                    A,...
                                                    cpolycoefs,...
                                                    cbest_order,...
                                                    cn_params,...
                                                    cIe,...
                                                    1:numel(dt(iStart:iEnd)),...
                                                    cexitflag,...
                                                    ieprior,...
                                                    stdprior,...
                                                    Ie0(:,end),...
                                                    out);
  % disp(['returned from recurse_AICfit, nt: ',num2str(numel(cnSteps))])
  ne0 = cneEnd(:,end);
  if iStart == 1
    ne_all         = cne;
    neEnd_all      = cneEnd;
    Ie_all         = cIe;
    polycoefs_all  = cpolycoefs;
    best_order_all = cbest_order;
    n_params_all   = cn_params;
    exitflag_all   = cexitflag;
    nSteps_all     = cnSteps;
  else
    ne_all         = [ne_all,         cne];
    neEnd_all      = [neEnd_all,      cneEnd];
    Ie_all         = [Ie_all,         cIe];
    polycoefs_all  = [polycoefs_all,  cpolycoefs];
    best_order_all = [best_order_all, cbest_order];
    n_params_all   = [n_params_all,   cn_params];
    exitflag_all   = [exitflag_all,   cexitflag];
    nSteps_all     = [nSteps_all,     cnSteps];
  end
  iStart = iStart + iIntervall;
end
% 5th collect output
out.polycoefs = polycoefs_all;
out.best_order = best_order_all;
out.params = n_params_all;
out.ne = ne_all;
out.neEnd = neEnd_all;
out.Ie = Ie_all;
out.exitflag = exitflag_all;
out.nSteps = nSteps_all;
% out.AICc = AICc;
% 6th calculate the covariances
ne0Cov = diag(mean(ppstd(:,1:Directives.ninteg).^2,2)/Directives.ninteg);
% error estimation for the coefficients from linearized theory
% $$$ xCovar = ElSpec_error_estimate_polycoefs( out.pp(:,1:Directives.ninteg), ...
% $$$                                           ppstd(:,1:Directives.ninteg), ...
% $$$                                           ne0, ...
% $$$                                           ne0Cov, ...
% $$$                                           A, ...
% $$$                                           out.alpha(:,1), ...
% $$$                                           out.dt(1:Directives.ninteg), ...
% $$$                                           out.Ec, ...
% $$$                                           out.dE, ...
% $$$                                           'neEnd', ...% TO-BE-FIXED!
% $$$                                           out.ieprior, ...
% $$$                                           out.stdprior, ...
% $$$                                           polycoefs1( best_order1(1), ...
% $$$                                                   1:best_order1(1)+1 , ...
% $$$                                                   1 ), ...
% $$$                                           numel(out.pp(:,1:Directives.ninteg)), ...
% $$$                                           out.ErrType, ...
% $$$                                           out.ErrWidth ...
% $$$                                           );

% $$$ it = 1;
% $$$ while it+out.nSteps(it)<=numel(out.te)
% $$$   
% $$$   xCovar = ElSpec_error_estimate_polycoefs( out.pp(:,it+[0:nSteps(it)-1]), ...
% $$$                                             out.ppstd(:,it+[0:nSteps(it)-1]),...
% $$$                                             ne0, ...   % These buggers
% $$$                                             ne0Cov,... % are tricky!
% $$$                                             A, ...
% $$$                                             out.alpha(:,it+[0:nSteps(it)-1]), ...
% $$$                                             out.dt(it+[0:nSteps(it)-1]), ...
% $$$                                             out.Ec, ...
% $$$                                             out.dE, ...
% $$$                                             out.integtype, ...
% $$$                                             out.ieprior, ...
% $$$                                             out.stdprior, ...
% $$$                                             out.polycoefs( ...
% $$$                                                 out.best_order(it), ...
% $$$                                                 1:(out.best_order(it)+1), ...
% $$$                                                 it ), ...
% $$$                                             numel(out.pp(:,it+[0:nSteps(it)-1])), ...
% $$$                                             out.ErrType, ...
% $$$                                             out.ErrWidth ...
% $$$                                             );
% $$$   
% $$$   % flux covariance matrix
% $$$   for tt = it:(it+out.nSteps(it)-1)
% $$$     out.IeCov(:,:,tt) = ElSpec_error_estimate_Ie( xCovar(1:(out.best_order(it)+1),...
% $$$                                                          1:(out.best_order(it)+1)), ...
% $$$                                                   out.polycoefs( ...
% $$$                                                       out.best_order(it), ...
% $$$                                                       1:(out.best_order(it)+1), ...
% $$$                                                       it), ...
% $$$                                                   out.Ec );
% $$$   end
% $$$   [out.neEndStd(:,(it+out.nSteps(it)-1)) , out.neEndCov(:,:,(it+out.nSteps(it)-1))]= ...
% $$$       ElSpec_error_estimate_ne2( ne0, ...
% $$$                                  out.polycoefs(out.best_order(it), 1:(out.best_order(it)+1),it), ...
% $$$                                  xCovar, ...
% $$$                                  sum(out.dt(it:(it+out.nSteps(it)-1))), ...
% $$$                                  A, ...
% $$$                                  out.Ec, ...
% $$$                                  out.dE(:), ...
% $$$                                  out.alpha(:,(it+out.nSteps(it)-1)), ...
% $$$                                  'endNe' );
% $$$ 
% $$$   it = it + out.nSteps(it);
% $$$ end
% 7 calculate FAC and power-flux and such

% ion production...
for it = numel(out.ts):-1:1
  out.q(:,it) = A*(out.Ie(:,it).*out.dE');
end

% chi-square of residuals
out.chisqr = mean(((out.ne - out.pp)./out.ppstd).^2);

% FAC
% integrate the spectra and multiply with e to get the current
% density
Eind_fac = out.Ec >= out.emin;

for it = numel(out.ts):-1:1
  out.FAC(it) = sum(out.Ie(Eind_fac,it).*out.dE(Eind_fac)')*1.60217662e-19;
% $$$   out.FACstd(it) = sqrt(out.dE(Eind_fac) * out.IeCov(Eind_fac,Eind_fac,it) ...
% $$$                         * out.dE(Eind_fac)')*1.60217662e-19;
end

% Power carried by the precipitating electrons
EdE = out.dE(Eind_fac)'.*out.Ec(Eind_fac)';
for it = 1:numel(out.ts)
  out.Pe(it) = sum(out.Ie(Eind_fac,it).*EdE)*1.60217662e-19;
% $$$   out.PeStd(it) = sqrt( EdE' * out.IeCov(Eind_fac, ...
% $$$                                          Eind_fac,it) ...
% $$$                         * EdE ) * 1.60217662e-19;
end
ElSpecOut = out;
try
  save(outfilename,'ElSpecOut')
catch
  disp(['Failed to save ElSpecOut into file:',outfilename])
end
end
