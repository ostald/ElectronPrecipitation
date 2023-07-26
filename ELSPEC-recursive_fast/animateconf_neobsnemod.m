function M = animateconf_neobsnemod(ElSpecRes,ops)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dops.ax = [5e9 1e12 95 150];
dops.aftsz = 12;
dops.lftsz = 12;
dops.tftsz = 15;
dops.plot_q = 1;
dops.pmsz = 12;
if nargin == 0
  M = dops;
  return
end
if nargin > 1 && isstruct(ops)
  dops = merge_structs(dops,ops);
end
axftsz = dops.aftsz;
lftsz = dops.lftsz;
tftsz = dops.tftsz;
plot_q = dops.plot_q;
pl_m_sz = dops.pmsz;
ax = dops.ax;

h = ElSpecRes.h;
ts = ElSpecRes.ts;
ne_mod = ElSpecRes.ne;
ne_obs = ElSpecRes.pp;
ne_obsstd = ElSpecRes.ppstd;
% q_mod = ElSpecRes.A*ElSpecRes.Ie;
q_mod = ElSpecRes.q;
M(numel(ts)) = getframe(gcf);
for iT = 1:23 % numel(ts)
  clf
  semilogx(ne_obs(:,iT),h,'b')
  hold on
  [hJPP,~] = JackKnifePlotY(ne_mod(:,iT),...
                            h(:),...
                            max(1,ne_mod(:,iT) - ne_obsstd(:,iT)*2),...
                            ne_mod(:,iT) + ne_obsstd(:,iT)*2);
  set(hJPP,'FaceColor',[1 0.7 0.7],'facealpha',0.5)
  [hJPP,hJPL] = JackKnifePlotY(ne_mod(:,iT),...
                               h(:),...
                               max(1,ne_mod(:,iT) - ne_obsstd(:,iT)),...
                               ne_mod(:,iT) + ne_obsstd(:,iT)*1);
  set(hJPP,'FaceColor',[1 0.5 0.5],'facealpha',0.5)
  set(hJPL,'color','r','linewidth',2)
  if plot_q
    phq = semilogx(q_mod(:,iT),h,'c-');
  else
    phq = [];
  end
  phnep = semilogx(ne_mod(:,max(1,iT-1)),h,'m-');
  phne = semilogx(ne_mod(:,iT),h,'r-','linewidth',2);
  phneO = semilogx(ne_obs(:,iT),h,'b.-','markersize',pl_m_sz);
  % legend([phneO,phne,phnep,phq],'n_e(t)','\tilde{n}_e(t)','\tilde{n}_e(t-dt)','q_e(t)','location','northwest','interpreter','latex')
  legend([phneO,phne,phnep,phq],...
         '$n_e(t)$',...
         '$\tilde{n}_e(t)$',...
         '$\tilde{n}_e(t-dt)$',...
         '$q_e(t)$',...
         'location','northwest',...
         'interpreter','latex',...
         'fontsize',axftsz)
  set(gca,'fontsize',axftsz)
  title(char(datetime(ts(iT),'ConvertFrom','posixtime','Format','yyyyMMdd HH:mm:ss')),'fontsize',tftsz)
  axis(ax)
  ylabel('height (km)','fontsize',lftsz)
  xlabel('electron density (/m^3)','fontsize',lftsz)
  %drawnow
  grid on
  M(iT) = getframe(gcf);
end

end