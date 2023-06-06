function ElSpec_IC_iter(iter, log_dir, ppdir)
%% ElSpec extended with Ion Chemitry
%% Based on Example-script for ElSpec_iqt usage
% This script should be possible to adapt for anyone with basic skills in
% matlab-usage. This test runs the ElSpec_iqt and ElSpec_qt functions
% testing the outlier-resilience with Pearson-type 6 statistics for the
% electron-density-estimates

%% 1, Setting up the matlab-path
% Simply modify the path below to point to the directory where ELSPEC-2022
% is installed. That should be it.
%addpath /bigdata/Campaigns/ELSPEC-2022 -end
%addpath 'ELSPEC-2022' -end

%% 2 Setup of parameters controlling ELSPEC

% Energy grid - between 10 ad 100 keV in 400 logarithmic/exponential steps
egrid = logspace(1,5,400);

% Data directories
% The paths to the directories with the ionospheric parameters and the
% power-profiles.
%fitdirUM = '/media/bgu001/5f5e8978-a828-4fd4-aabf-2032a3fb895b/Data/EISCAT/Analysed/2006-12-12_arc1_4@uhf';
%fitdir = '/mnt/data/bjorn/EISCAT/Analysed/2006-12-12_arc1_4@uhf';
fitdir = '../Data/Eiscat/fit';
%ppdirUM = '/media/bgu001/5f5e8978-a828-4fd4-aabf-2032a3fb895b/Data/EISCAT/tmp-ionlines/2006-12-12_arc1_4@uhf-pp';
%ppdir = '/mnt/data/bjorn/EISCAT/Analysed/2006-12-12_arc1_4@uhf-pp';
%ppdir = '../Data/Eiscat/pp';
% Flag for specifying which EISCAT-experiment it is
experiment = 'arc1';
% Altitude-limits.
hmax = 150;
hmin = 95;
% Time-limits
btime = [2006, 12, 12, 19, 30, 0];
etime = [2006, 12, 12, 19, 35, 0];
% Selection of which ionisation-profile method to use
ionomodel = 'Sergienko';
%recombmodel = 'SheehanGrFlipchem';
% and which type of continuity-integration-method to use
integtype = 'integrate';
% Time-resolution to use
tres = 'best';
% Maximum order of polynomial to try. The electron-flux is calculated as an
% the product of the energy with an exponential of an n-th order polynomial:
% Ie(E) = E*exp(p_n(E))
maxorder = 5;
% ninteg set to 5, here that only means that this is the number of
% time-steps to integrate for the initial condition
ninteg = 20;
% For additional settings see the help of ElSpec, ElSpec_iqt and ElSpec_qt

%% Outlier-specification
% The arrays here specify the outliers, and are of the format 
%              iz1 iz2 it1 it2 C
    % OUTLIERS{1} = [ 12  13  91  92 10;
    %                12 21 123 124 10;
    %                12 21 421 422 10];
    % OUTLIERS{2} = [12 21 91 92 10;
    %                12 21 123 124 30;
    %                12 21 421 422 30];
    % OUTLIERS{3} = [12 13 91 92 30;
    %                12 51 123 124 10;
    %                12 51 421 422 3];
    % OUTLIERS{4} = [12 21 91 92 5;
    %                12 51 123 124 3;
    %                12 51 421 422 2];
    % OUTLIERS{5} = [12 21 91 92 2;
    %                12 51 123 124 2;
    %                12 51 421 422 2];
% and the electron-densities at altitudes between z(iz1) and z(iz2) at
% times between t(it1) and t(it2) will be multiplied with C, the last
% element on each row. Here 3 outliers are inserted in the power-profiles.

% Pearson-type 6 statistics assumed
ErrType = 'l'; % L for Lorentzian.
% for i1 = numel(OUTLIERS):-1:1
%   Outliers = OUTLIERS{i1};
%   Outname = sprintf('ElSpec-Outliers-L-5-%02i.mat',i1);
%   disp(Outname)
%   disp(Outliers)
%   %                                       'Tdir',Tdir,...  (time direction)
%   ElSpecQT_Outliers_L5{i1} = ElSpec_qt('fitdir',fitdir,...
%                                        'ppdir',ppdir,...
%                                        'experiment',experiment,...
%                                        'hmax',hmax,'hmin',hmin,...
%                                        'btime',btime,'etime',etime,...
%                                        'ionomodel',ionomodel,...
%                                        'integtype',integtype,...
%                                        'egrid',egrid,...
%                                        'tres',tres,...
%                                        'ErrType',ErrType,...
%                                        'MaxOrder',maxorder,...
%                                        'ninteg',ninteg,...
%                                        'Outliers',Outliers,...
%                                        'Outfilename',Outname);
%   ElSpecPlot(ElSpecQT_Outliers_L5{i1})
%   [fnm1,fnm2,fnm3] = fileparts(ElSpecQT_Outliers_L5{i1}.Outfilename); 
%   print('-depsc','-vector',[fnm2,'-QT-l-5']);
%   dstr = sprintf('Done with loop S i1: %i at %s',i1,datestr(now,'HH:MM:SS'));
%   disp(dstr)
%   close(gcf)
%   
% end
% 
% for i1 = numel(OUTLIERS):-1:1
%   Outliers = OUTLIERS{i1};
%   Outname = sprintf('ElSpec-iqtOutliers-L-5-%02i.mat',i1);
%   disp(Outname)
%   disp(Outliers)
%   ElSpecQT_iqtOutliers_L5{i1} = ElSpec_iqt('fitdir',fitdir,...
%                                        'ppdir',ppdir,...
%                                        'experiment',experiment,...
%                                        'hmax',hmax,'hmin',hmin,...
%                                        'btime',btime,'etime',etime,...
%                                        'ionomodel',ionomodel,...
%                                        'integtype',integtype,...
%                                        'egrid',egrid,...
%                                        'tres',tres,...
%                                        'ErrType',ErrType,...
%                                        'MaxOrder',maxorder,...
%                                        'ninteg',ninteg,...
%                                         'Outliers',Outliers,...
%                                        'Outfilename',Outname);
%   ElSpecPlot(ElSpecQT_iqtOutliers_L5{i1})
%   [fnm1,fnm2,fnm3] = fileparts(ElSpecQT_iqtOutliers_L5{i1}.Outfilename) ;
%   print('-depsc','-vector',[fnm2,'-IQT-l-5']);
%   dstr = sprintf('Done with loop S i1: %i at %s',i1,datestr(now,'HH:MM:SS'));
%   disp(dstr)
%   close(gcf)
%   
% end

%try soemthing with same partitioning:
% if iter > 1
%     m2 = iter - 2;
%     nStepsm2 = load(fullfile('..', log_dir, ["/ElSpec-iqt_IC_" + m2 + ".mat"])).ElSpecOut.nSteps;
%     m1 = iter - 1;
%     nStepsm1 = load(fullfile('..', log_dir, ["/ElSpec-iqt_IC_" + m1 + ".mat"])).ElSpecOut.nSteps;
%     if nStepsm2 == nStepsm1
%         nstep = nStepsm2;
%         disp('nstep copied')
%     else
%         nstep = 1;
%     end
% else
%     nstep = 1;
% end

if iter > 0
    j = iter - 1;
    icdir = fullfile('..',log_dir,["IC_" + j + ".mat"]);
    %icdir = "../" + log_dir + "IC_" + j + ".mat";
    icdata = load(icdir);
    iri_ic = icdata.elspec_iri_sorted;
    alpha_eff = icdata.eff_rr;
    ne_init = icdata.ne_init;

    elspec_m1 = fullfile('..',log_dir,["ElSpec-iqt_IC_" + j + ".mat"]);
    nsteps_old = load(elspec_m1).ElSpecOut.nSteps;
    %ninteg = nsteps_old(1); did not result in converging behaviour
else
    iri_ic = 0;
    alpha_eff = 0;
    ne_init = 0;
end

%Outname = sprintf('../' + log_dir + 'ElSpec-iqt_IC_0');
Outname = fullfile("..", log_dir, ["ElSpec-iqt_IC_" + iter]);
disp(Outname)


ElSpecQT_iqtOutliers_L5 = ElSpec_iqt_ic('fitdir',fitdir,...
                                       'ppdir',ppdir,...
                                       'experiment',experiment,...
                                       'hmax',hmax,'hmin',hmin,...
                                       'btime',btime,'etime',etime,...
                                       'ionomodel',ionomodel,...
                                       'integtype',integtype,...
                                       'egrid',egrid,...
                                       'tres',tres,...
                                       'ErrType',ErrType,...
                                       'MaxOrder',maxorder,...
                                       'ninteg',ninteg,...
                                       'Outfilename',Outname,...
                                       'iri_ic', iri_ic, ...
                                       'alpha_eff', alpha_eff, ...
                                       'iteration', iter, ...
                                       'ne_init', ne_init);

%with recombmodel:
% ElSpecQT_iqtOutliers_L5 = ElSpec_iqt_ic('fitdir',fitdir,...
%                                        'ppdir',ppdir,...
%                                        'experiment',experiment,...
%                                        'hmax',hmax,'hmin',hmin,...
%                                        'btime',btime,'etime',etime,...
%                                        'ionomodel',ionomodel,...
%                                        'recombmodel', recombmodel, ...
%                                        'integtype',integtype,...
%                                        'egrid',egrid,...
%                                        'tres',tres,...
%                                        'ErrType',ErrType,...
%                                        'MaxOrder',maxorder,...
%                                        'ninteg',ninteg,...
%                                        'Outfilename',Outname,...
%                                        'iri_ic', iri_ic, ...
%                                        'alpha_eff', alpha_eff, ...
%                                        'iteration', iter, ...
%                                        'ne_init', ne_init);
ElSpecPlot(ElSpecQT_iqtOutliers_L5, ieelim = [10, 14], faclim = [0, 10], plim = [0, 60]);
[fnm1,fnm2,fnm3] = fileparts(ElSpecQT_iqtOutliers_L5.Outfilename) ;
print('-depsc','-vector',[Outname]);
print('-dpdf', '-painters', [Outname]);
dstr = sprintf('Done with loop S i1: 0 at %s',datestr(now,'HH:MM:SS'));
disp(dstr)
close(gcf)
  

%% Two modified plotting-functions
% In addition to ElSpecPlot, ElSpecPlot2 and ElSpecPlotSmall 2 new
% plotting-functions are added:
%fig1 = figure;
%ElSpecPlotIeNePpRes(ElSpecQT_iqtOutliers_L5,fig1.Number);
%fig2 = figure;
%ElSpecPlotRes(ElSpecQT_iqtOutliers_L5,fig2.Number);
%ElSpecPlotIeNePpRes(ElSpecQT_iqtOutliers_L5);
%ph = ElSpecPlotTres(ElSpecQT_iqtOutliers_L5);
%set(ph,'color','r','linestyle','--','marker','.');
end
