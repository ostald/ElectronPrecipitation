clf

if 1
  sph211 = subplot(2,1,1);
  pcolor(rem(ts/3600,24),Ec/1e3,log10(ElSpecQT_iqtOutliers_L5{i1}.IeStd./(ElSpecQT_iqtOutliers_L5{i1}.Ie+1))),shading flat
  set(gca,'YScale','log','fontsize',12,'TickDir','both')     
  set(sph211,'Position',get(sph211,'Position')+[-0.03 0 0 0])
  axis([ax5t(1:2) 0.5 50])
  timetick
  ylabel('E (eV)')
  caxis([-2 1 ])
  cblh1 = colorbar_labeled('\sigma_{I_e}/I_e','log','fontsize',12);
end

if 1
  sph212 = subplot(2,1,2);
  pcolor(rem(ts/3600,24),h,(real(ElSpecQT_iqtOutliers_L5{i1}.neEndStd)./(ElSpecQT_iqtOutliers_L5{i1}.ne))),shading flat
  caxis([0 .15 ])
  set(gca,'fontsize',12,'TickDir','both')     
  set(sph212,'Position',get(sph212,'Position')+[-0.03 0 0 0])
  axis([ax5t(1:2) 95 150])
  timetick
  ylabel('height (km)')
  xlabel('time (UT)')
  cblh1 = colorbar_labeled('\sigma_{n_e}/n_e','linear','fontsize',12);
end

