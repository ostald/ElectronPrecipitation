function ph = plotconf_neobsnemod2(ElSpecRes,iT2plot,ops)
% plotconf_neobsnemod2 - observed and modeled electron-densities +/- sigma
%   plotconf_neobsnemod plots the observed and ELSPEC-modeled electron
%   densities with the 1 and 2 sigma confidence-intervalls estimated from
%   the observations, together with panels of ionization-profiles and
%   the ionization-produciton-matrix and the primary electron-spectra. The
%   confidence-intervalls are plotted in a JackKnife-style around the
%   modeled electron-density, which is suitable considering that the
%   measurements are each one sample out of some random-distributed
%   distribution with some expected value, for which we have the modeled
%   electron-density and some spread - that we have a standard-deviation
%   for. Since we have close-to-normal-distributed random distributed
%   measurements it is appropriate to look at the modeled electron-density
%   as an estimate of the expected values of the electron densities with
%   some normal-distributed probability for the measurements with a width
%   of the standard deviation. 
% 
% Calling:
%  ph = plotconf_neobsnemod2(ElSpecRes,iT2plot,ops)
%  ops = plotconf_neobsnemod2
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
  ph = dops;
  return
end
% If there is a user-supplied options-struct then replace the
% default-fields with the user-options 
if nargin > 2 && isstruct(ops)
  dops = merge_structs(dops,ops);
end
axftsz = dops.aftsz;
lftsz = dops.lftsz;
tftsz = dops.tftsz;
pl_m_sz = dops.pmsz;
ax = dops.ax;
facealpha = dops.facealpha;

%% Extract the ElSpec-results
h = ElSpecRes.h;
ts = ElSpecRes.ts;
ne_mod = ElSpecRes.ne;
ne_obs = ElSpecRes.pp;
ne_obsstd = ElSpecRes.ppstd;
q_mod = ElSpecRes.q;
Ie = ElSpecRes.Ie;
A = ElSpecRes.A;
A(~isfinite(A(:))) = 0;
E = ElSpecRes.E;
%% Get to plotting
iT = iT2plot(1); % safe-guard against array-input for iT2plot
% plotting layout
sphIe = subplot(3,5,1:2);
sph_A = subplot(3,5,[6 7 11 12]);
sph_q = subplot(3,5,[8 13]);
sphne = subplot(3,5,[9 10 14 15]);


%% Electron-spectral panel
axes(sphIe) %  = subplot(3,5,1:2);
loglog(E/1e3,Ie([1:end,end],iT))
axis([0.08 50 1e6 1e12])
ylabel('I_e (/m^2/s/eV)','fontsize',lftsz)
title('Electron flux','fontsize',lftsz)
grid on
set(gca,'YMinorGrid','off','XMinorGrid','on','XTickLabel','')

%% Ionization-profile-matrix panel
axes(sph_A) %  = subplot(3,5,[6 7 11 12]);
pcolor(E/1e3,h,log10(A(:,[1:end,end]))),shading flat
set(gca,'XScale','log')
xlabel('Energy (keV)','fontsize',lftsz)
ylabel('height (km)','fontsize',lftsz)
title('log10(ionization profile matrix)','fontsize',tftsz)
caxis([-6 0]+max(caxis)) % clim is too new. Fuck off with the whining!
set(gca,'fontsize',axftsz,'tickdir','both')
axis([0.08 50 95 150])

%% Ionization-profile panel
axes(sph_q) % = subplot(3,5,[8 13]);
phq = semilogx(q_mod(:,iT),h,'c-');
title('Ionization-profile','fontsize',tftsz)
xlabel('q_e (/m^3/s)','fontsize',lftsz)
set(gca,'YTickLabel','','fontsize',axftsz,'tickdir','both','box','off')
axis([1e6 1e12 95 150])
set(gca,'XTick',[1e6 1e9 1e12])
grid on

%% Electron-density panel
axes(sphne) % = subplot(3,5,[9 10 14 15]);
ph2d = semilogx(ne_obs(:,iT),h,'b');
hold on
[hJPP1,~] = JackKnifePlotY(ne_mod(:,iT),...
                           h(:),...
                           max(1,ne_mod(:,iT) - ne_obsstd(:,iT)*2),...
                           ne_mod(:,iT) + ne_obsstd(:,iT)*2);
set(hJPP1,'FaceColor',[1 0.7 0.7],'facealpha',facealpha)
[hJPP2,hJPL2] = JackKnifePlotY(ne_mod(:,iT),...
                               h(:),...
                               max(1,ne_mod(:,iT) - ne_obsstd(:,iT)),...
                               ne_mod(:,iT) + ne_obsstd(:,iT)*1);
set(hJPP2,'FaceColor',[1 0.5 0.5],'facealpha',facealpha)
set(hJPL2,'color','r','linewidth',2)
set(gca,'YTickLabel','','tickdir','both','box','off')
phnep = semilogx(ne_mod(:,max(1,iT-1)),h,'m-');
phne = semilogx(ne_mod(:,iT),h,'r-','linewidth',2);
phneO = semilogx(ne_obs(:,iT),h,'b.-','markersize',pl_m_sz);
legend([phneO,phne,phnep,phq],...
       '$n_e(t)$',...
       '$\tilde{n}_e(t)$',...
       '$\tilde{n}_e(t-dt)$',...
       'location','northwest',...
       'interpreter','latex',...
       'fontsize',axftsz)
set(gca,'fontsize',axftsz,'YTickLabel','')
title('electron density','fontsize',tftsz)
axis(ax)
% ylabel('height (km)','fontsize',lftsz)
xlabel('electron density (/m^3)','fontsize',lftsz)
grid on
slh = suplabel(char(datetime(ts(iT),'ConvertFrom','posixtime','Format','yyyy-MM-dd hh:mm:ss.SS')),'t');
set(slh,'fontsize',lftsz)
if nargout == 1
  ph = [phneO,phne,phnep,hJPL2,hJPP2,hJPP1];
end

delete(ph2d)
end