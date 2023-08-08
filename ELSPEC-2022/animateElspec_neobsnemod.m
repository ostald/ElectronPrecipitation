function M = animateElspec_neobsnemod(ElSpecRes,iT2plot,ops)
% plotconf_neobsnemod - observed and modeled electron-densities +/- sigma
%   plotconf_neobsnemod plots the observed and ELSPEC-modeled electron
%   densities with the 1 and 2 sigma confidence-intervalls estimated from
%   the observations. Optionally one can plot the estimated
%   ionization-profile as well. The confidence-intervalls are plotted in a
%   JackKnife-style around the modeled electron-density, which is suitable
%   considering that the measurements are each one sample out of some
%   random-distributed distribution with some expected value, for which we
%   have the modeled electron-density and some spread - that we have a
%   standard-deviation for. Since we have close-to-normal-distributed
%   random distributed measurements it is appropriate to look at the
%   modeled electron-density as an estimate of the expected values of the
%   electron densities with some normal-distributed probability for the
%   measurements with a width of the standard deviation.
% 
% Calling:
%  ph = plotconf_neobsnemod(ElSpecRes,iT2plot,ops)
%  ops = plotconf_neobsnemod
% Input:
%  ElSpecRes - struct with ElSpec-results as returned from the ElSpec,
%              ElSpec_iqt and ElSpec_iqtcl functions
%  iT2plot   - index in time-dimension to plot the altitude-proofiles for,
%              scalar integer.
%  ops       - struct with options controlling the plot, struct with
%              default options are retured when plotconf_neobsnemod is
%              called without input arguments.
% Output:
%  ph - handles to lines and patches generated, the patch-handles are last
%       in the array.
%  ops - struct with default options, returned when calling the function
%        without input-arguments.

%% Default options:
dops.ax = [5e9 1e12 95 150];
dops.aftsz = 12;
dops.lftsz = 12;
dops.tftsz = 15;
dops.pmsz = 12;
dops.facealpha = 0.4;
if nargin == 0
  % return the def-ops when there are no input arguments
  M = dops;
  return
end
% If there is a user-supplied options-struct then replace the
% default-fields with the user-options 
if nargin > 2 && isstruct(ops)
  dops = merge_structs(dops,ops);
end
if nargin < 2 || isempty(iT2plot)
  iT2plot = 1:numel(ElSpecRes.ts);
end

M(numel(iT2plot)) = getframe(gcf);
for it = 1:numel(iT2plot)
  iT = iT2plot(it);
  clf
  plotconf_neobsnemod2(ElSpecRes,iT,dops);
  M(it) = getframe(gcf);
end