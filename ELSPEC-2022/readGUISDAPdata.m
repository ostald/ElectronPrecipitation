function [h,ts,te,pp,ppstd,par,parstd,loc,azel,I] = readGUISDAPdata( ppdir , ...
                                                  fitdir , hmin , ...
                                                  hmax , tmin , tmax ...
                                                  , exp , radar , ...
                                                      version , tres , FAdev  ...
                                                      )
%
% Read GUISDAP raw densities (power profiles) and fit results from
% files in ppdir and fitdir.
%
% [h,ts,te,pp,ppstd,par,parstd] =
%    readGUISDAPdata( ppdir , fitdir , hmin , hmax , tmin , tmax , ...
%    exp , radar , version , tres , FAdev)
%
%
% INPUT:
%  ppdir   path to GUISDAP power profile data
%  fitdir  path to GUISDAP fit results
%  hmin    minimum height [km]
%  hmax    maximum height [km]
%  tmin    begin time, array [year, month, day, hour , minute ,
%          second], use empty array to read all data
%  tmax    end time, array [year, month, day, hour , minute ,
%          second], use empty array to read all data
%  exp     EISCAT experiment name
%  radar   EISCAT radar name ['U','V']
%  version EISCAT experiment version number [1,2,3,...]
%  tres    "type" of time resolution 'best' or 'dump'
%  FAdev   maximum beam direction deviation from field-aligned  [deg]
%  azel    azimuth and elevation of the radar beam
%  I       magnetic inclination angle (deg)
%
%
% OUTPUT:
%  h       heights [km]
%  ts      integration period start times [unixtime]
%  te      integration period end times [unixtime]
%  pp      raw densities, length(h)xlength(ts) matrix [m^-3]
%  ppstd   standard deviations of the raw densities [m^-3]
%  par     lenth(h)x4xlength(ts) array of plasma parameters
%  parstd  standard deviations of the plasma parameters
%
% The four plasma parameters are Ne [m^-3], Ti [K], Te [k], and Vi
% [ms^-1]. Failed iterations and results with chi-squared larger
% than 10 are replaced with NaNs.
%
% IV 2016 - 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% a special case for empty ppdir
if isempty(ppdir)
    [hpar,ts,te,par,parstd,loc,azel,I] = readGUISDAPpar( fitdir , FAdev );
    if isempty(hpar)
        h = [];
        pp = [];
        ppstd = [];
        return
    end
    hind = hpar(:,1)>=hmin & hpar(:,1)<=hmax;
    if isempty(tmin)
        t1 = -Inf;
    else
        t1 = date2unixtime(tmin);
    end
    if isempty(tmax)
        t2 = Inf;
    else
        t2 = date2unixtime(tmax);
    end
    it = find((ts >= t1) & (te <= t2));
    ts = ts(it);
    te = te(it);
    h = hpar(hind,1)';
    par = par(hind,:,it);
    parstd = parstd(hind,:,it);
    pp = squeeze(par(:,1,:));
    ppstd = squeeze(parstd(:,1,:));


    % the guisdap error estimates are not good enough for us, try
    % to calculate a sample average
    [d1 d2 d3] = size(parstd);
    for i1 = 1:d1
        for i3 = 1:d3
            tmpne = par(max(1,i1-1):min(d1,i1+1),1,max(1,i3-1): ...
                        min(d3,i3+1));
            ppstd(i1,i3) = std(tmpne(:),'omitnan');
        end
    end
    
    return
end

% read power profiles
[hpp,tspp,tepp,pp1,ppstd1,locpp,azelpp,Ipp] = readGUISDAPpp( ppdir , exp , radar , ...
                                            version , tres , FAdev );

% read fit results
[hpar,tspar,tepar,par1,parstd1,locpar,azelpar,Ipar] = readGUISDAPpar( fitdir , FAdev );

if isempty(hpp) | isempty(hpar)
    h = [];
    ts = [];
    te = [];
    pp = [];
    ppstd = [];
    par = [];
    parstd = [];
    loc = [];
    azel = [];
    I = [];
    return
end

% pick the location
loc = locpp;
if isempty(loc)
    loc = locpar;
end
if isempty(loc)
    switch lower(radar)
      case 'uhf'
        loc = [69.583 19.21 0.03];
      case 'vhf'
        loc = [69.583 19.21 0.03];
      case 'esr'
        loc = [78.1530   16.0290    0.4380];
      case '42m'
        loc = [78.1530   16.0290    0.4380];
      case '32m'
        loc = [78.1530   16.0290    0.4380];
      otherwise
        error(['radar location not found from data files and unknown radar site' radar])
    end
end
% azimuth and elevation
azel = azelpp;
% magnetic inclination
I = Ipp;

% now we have data in different grids 

% first pick the most common height step from hpp, this will be our
% height resolution
hstps = diff( hpp );
hstps = hstps(:);
hstps(isnan(hstps)) = [];
hstp = mode( hstps );

% then pick the most common first height, this is where we will
% start from
hbegs = hpp(1,:);
hbegs(isnan(hbegs)) = [];
hbeg = mode( hbegs );

% then the most common largest height
hmaxs = max( hpp );
hmaxs(isnan(hmaxs)) = [];
hend = mode(hmaxs);

% then finally the heights
h = hbeg:hstp:hend;
h = h(h>=hmin);
h = h(h<=hmax);
nh = length(h);

% use the times from power profiles
ts = tspp;
te = tepp;
if isempty(tmin)
    t1 = -Inf;
else
    t1 = date2unixtime(tmin);
end
if isempty(tmax)
    t2 = Inf;
else
    t2 = date2unixtime(tmax);
end
it = find((ts >= t1) & (te <= t2));
ts = ts(it);
te = te(it);
nt = length(te);

% interpolate everything to the height grid h and pick only
% time-slices it
pp = NaN(nh,nt);
ppstd = NaN(nh,nt);
for k=1:nt
    % remove NaN values from pp before interpolation
    ii = ~isnan(hpp(:,it(k)));
    pp(:,k) = interp1(hpp(ii,it(k)),pp1(ii,it(k)),h);
    ppstd(:,k) = interp1(hpp(ii,it(k)),ppstd1(ii,it(k)),h);
end

% the fit results, the fit results can be optionally skipped with
% empty fitdir.
par = NaN(nh,4,nt);
parstd = NaN(nh,4,nt);

if ~isempty(fitdir)
    % ntpar = length(tepar);
    % par2 = NaN(nh,4,ntpar);
    % parstd2 = NaN(nh,4,ntpar);
    % %for every timestep, interpolate in height:
    % for k=1:ntpar 
    %     ii = ~isnan(hpar(:,k));
    %     % interpolate, some tricks at edges (i.e. extra repetition in the
    %     % front and the end (?? there might be a bug => fixed))
    %     p_curr = par1(ii,:,k);
    %     p_currstd = parstd1(ii,:,k);
    %     par2(:,:,k) = interp1([0;hpar(ii,k);1e6],p_curr([1,1:end,end],:),h,'linear','extrap'); 
    %     parstd2(:,:,k) = interp1([0;hpar(ii,k);1e6],p_currstd([1,1:end,end],:),h,'linear','extrap');
    %     % before:
    %     %par2(:,:,k) =
    %     %interp1([0;hpar(ii,k);1e6],[par1(ii(1),:,k);par1(ii,:,k);par1(ii(end),:,k)],h,'linear','extrap');         %is that correct? i think par1(ii(end),:,k) may be wrong, bc indexing by array does not work if you specify one element in the array: ii(end)
    %     %parstd2(:,:,k) = interp1([0;hpar(ii,k);1e6],[parstd1(ii(1),:,k);parstd1(ii,:,k);parstd1(ii(end),:,k)],h,'linear','extrap');
    % end
    % %then, for every height, interpolate in time:
    % for k=1:nh
    %     for l=1:4
    %         %replaced te, tepar by tmid, tmidpar to use central points
    %         par(k,l,:) = interp1([0;tmidpar;1e20],[squeeze(par2(k,l,1));squeeze(par2(k,l,:));squeeze(par2(k,l,end))],tmid,'linear','extrap');
    %         parstd(k,l,:) = interp1([0;tmidpar;1e20],[squeeze(parstd2(k,l,1));squeeze(parstd2(k,l,:));squeeze(parstd2(k,l,end))],tmid,'linear','extrap');
    %     end
    % end

    
    %replace interpolation with 2d interpolation for a smooth
    %altitude profile

    % first, interpolate in time to fill all NANs with finite values
    % go wide with clipping, so that 2d interpolation later has no need to
    % extrapolate
    tmidpar = (tspar + tepar)/2;
    i_tpar_min = find(tmidpar>=t1, 1)-1;
    i_tpar_max = find(tmidpar<=t2, 1, 'last')+1;
    tpar_ = tmidpar(i_tpar_min:i_tpar_max);
    ntpar_ = length(tpar_);
    
    nhpar = size(hpar, 1);
    
    par3 = NaN(nhpar,4,ntpar_);
    parstd3 = NaN(nhpar,4,ntpar_);
    for k=1:nhpar
        for l=1:4
            ii = ~isnan(squeeze(par1(k, l, :)));
            p_curr = squeeze(par1(k,l,ii));
            p_currstd = squeeze(parstd1(k,l,ii));
            par3(k,l,:) = interp1([0;tmidpar(ii);1e20],[p_curr(1); p_curr; p_curr(end)],tpar_,'pchip','extrap');
            parstd3(k,l,:) = interp1([0;tmidpar(ii);1e20],[p_currstd(1); p_currstd; p_currstd(end)],tpar_,'pchip','extrap');
        end
    end

    %Then interpolate to a uniform, coarse hpar_ ≈ hpar (might not be
    %necessary, since relative deviations in hpar are <1%)
    hpar_ = mean(hpar, 2);
    i_hpar_min = find(hpar_>=hmin, 1)-1;
    i_hpar_max = find(hpar_<=hmax, 1, 'last')+1;
    hpar_ = hpar_(i_hpar_min:i_hpar_max);
    nhpar_ = length(hpar_);

    par2 = NaN(nhpar_,4,ntpar_);
    parstd2 = NaN(nhpar_,4,ntpar_);
    for k=1:ntpar_ 
        ii = ~isnan(hpar(:,k));
        % interpolate, some tricks at edges (i.e. extra repetition in the
        % front and the end (?? there might be a bug => fixed with p_curr)
        p_curr = par3(ii,:,k);
        p_currstd = parstd3(ii,:,k);
        par2(:,:,k) = interp1([0;hpar(ii,k);1e6],p_curr([1,1:end,end],:),hpar_,'pchip','extrap'); 
        parstd2(:,:,k) = interp1([0;hpar(ii,k);1e6],p_currstd([1,1:end,end],:),hpar_,'pchip','extrap');
        
    end

    %finally, do 2d interpolation to upscale to desired resolution
    tmid = (ts+te)/2;
    [X, Y] = meshgrid(tpar_, hpar_);
    [Xq, Yq] = meshgrid(tmid, h);
    for l=1:4
        par(:, l, :) = interp2(X, Y, squeeze(par2(:, l, :)), Xq, Yq, 'makima');
        parstd(:, l, :) = interp2(X, Y, squeeze(parstd2(:, l, :)), Xq, Yq, 'makima');
    end
end

end
