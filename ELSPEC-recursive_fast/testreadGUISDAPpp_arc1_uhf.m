function [h,ts,te,pp,ppstd,loc,azel] = testreadGUISDAPpp_arc1_uhf( ff , version ...
                                                      , tres )
%
% [h,ts,te,pp,ppstd,loc] = readGUISDAPpp_arc1_uhf( ff , version )
%
% Read raw electron densities from GUISDAP UHF arc1 result files.
% Sufficient power profiles are available ONLY IF the analysis was
% run with parameters analysis_ppshortlags=1 and analysis_pponly=1
%
% This function is based on guisdap_pp_reader by Björn Gustavsson.
%
% INPUT:
%   ff       a char vector of file names
%   version  exp version number, ignored at the moment
%   tres     "type" of time resolution 'best' or 'dump'
%
% OUTPUT:
%  h     heights [km]
%  ts    integration start times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  te    integration end times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  pp    a matrix of electron density values [m^-3]
%  ppstd standard deviations of the electron densities [m^-3]
%  loc   [latitude (deg north) longitude (deg east) height (km)]
%        of the radar site
%  azel  beam azimuth and elevation angles
%
% IV 2016

switch version
  case 1
    nprofstep = 3;
  case 2
    nprofstep = 2;
  otherwise
    error('unknown version number')
end

%nprofstep = 1; % BG-20221004: I dont understand the nprofstep-variable?
% load the first file to get the matrix sizes
load(ff{1});
i_end = [find(diff(r_pprange)<0);numel(r_pprange)];% Should give all ends also the last one? BG-20221004
i_start = [0; i_end(1:end-1)]+1;
%i_end = i_end(1:nprofstep:end); % there are N profiles, we want only the
                        % first one...
nh_all = unique(i_end-i_start,'stable')+1;
nf = length(ff);

loc = r_XMITloc;
for iR = numel(nh_all):-1:1
  nh = nh_all(iR);
  ppmat{iR} = NaN(nh,length(i_end)/numel(nh_all));

  % note to self: this is the only place where both options are
  % checked. After this if-statement, we asssume that tres is 'best'
  % if it is not 'dump'. Remember to change the other if-statements
  % if other options are implemented.
  if lower(tres)=='dump'
    pp{iR} = NaN(nh,nf);
    ppstd{iR} = NaN(nh,nf);
    ts = NaN(nf,1);
    te = NaN(nf,1);
    h{iR} = NaN(nh,nf);
    azel = NaN(nf,2);
  elseif lower(tres)=='best'
    nppdump = (length(i_end)/numel(nh_all));
    npp = nf * nppdump;
    pp{iR} = NaN(nh,npp);
    ppstd{iR} = NaN(nh,npp);
    ts = NaN(npp,1);
    te = NaN(npp,1);
    h{iR} = NaN(nh,npp);
    azel = NaN(npp,2);
  else
    error(['Unknown time resolution type',tres])
  end
end
idx_profs = ones(size(nh_all));
for k = 1:length(ff)

    load(ff{k})

    % there are several profiles in the file, find their ends
    i_end = [find(diff(r_pprange)<0);numel(r_pprange)]; % was (BG-20221004): find(diff(r_pprange)<0);
    i_start = [1;i_end(1:end-1)+1];
%     % start indices of profiles

    % start and end times of integration, converted into unix time
    if lower(tres) == 'dump'
      ts(k) = date2unixtime(r_time(1,:));
      te(k) = date2unixtime(r_time(2,:));
    else
      ts( (k-1)*nppdump+1 ) = date2unixtime(r_time(1,:));
      te( k*nppdump ) = date2unixtime(r_time(2,:));
      dtpp = (te( k*nppdump ) - ts( (k-1)*nppdump+1 ) ) / nppdump;
      for kpp = 2:nppdump
          ts( (k-1)*nppdump+kpp ) = ts( (k-1)*nppdump+1 ) + (kpp-1)*dtpp;
          te( k*nppdump-kpp+1 ) = te( k*nppdump ) - (kpp-1)*dtpp;
      end
    end
    % From here iterate forward over the individual power-profiles
    % and assign the heights, pp, ppstd to the respective elements of the
    % range-groups. This will not exactly work for 'dump' resolution in
    % time, but shouldn't be bad for 'best'
    if lower(tres)=='best'
      for iDX = 1:numel(i_end)
        curr_idx = i_start(iDX):i_end(iDX);
        curr_set = find(numel(curr_idx)==nh_all);
        h{curr_set}(:,idx_profs(curr_set)) = r_pprange(curr_idx);
        pp{curr_set}( : ,idx_profs(curr_set)) = r_pp(curr_idx);
        ppstd{curr_set}( :,idx_profs(curr_set)) = r_pperr(curr_idx);
        idx_profs(curr_set) = idx_profs(curr_set)+1;
      end

    % heights and beam directions
    %if lower(tres) == 'dump'
    %else
      %h( : , (k-1)*nppdump+1 ) = r_pprange(i2(1):i3(1))*sin(r_el*pi/ ...
      %                                                  180);
      azel((k-1)*nppdump+1,:) = [r_az r_el];
      for kpp=2:nppdump
          %h( : , (k-1)*nppdump+kpp ) = h( : , (k-1)*nppdump+1 );
          azel((k-1)*nppdump+kpp,:) = [r_az r_el];
      end
    else % then lower(tres)=='dump'
      i2 = [1;i_end(nprofstep:nprofstep:end-1)+1];
      % end indices of profiles
      i3 = [i_end(1:nprofstep:end)];%length(r_pp)];
      % read the power profiles
      for i4 = 1:length(i2),
        ppmat(:,i4) = r_pp(i2(i4):i3(i4));
      end
      h(:,k)  = r_pprange(i2(1):i3(1))*sin(r_el*pi/180);
      azel(k,:) = [r_az r_el];
      % integrate the profiles together
      pp(:,k) = mean(ppmat,2);
      ppstd(:,k) = std(ppmat,0,2)./sqrt(size(ppmat,2));
    %else
    %    % keep the profiles separate
    %    pp( : , ((k-1)*nppdump+1) : (k*nppdump) ) = ppmat;
    %    ppstd( : , (k-1)*nppdump+1 ) = std(ppmat,0,2);
    %    for kpp = 2:nppdump
    %        ppstd( : , (k-1)*nppdump+kpp ) = ppstd( : , (k-1)*nppdump+1 );
    %    end
    end

end

if lower(tres)=='best'
  % Then we need to merge the multiple power-profiles. Here we use standard
  % variance-weighted averages for all ranges where we have multiple
  % samples of the power-profile. That is:
  % 
  % PP = sum(pp/ppstd.^2)/sum(1/ppstd.^2)
  % 
  % and
  %
  % ppstd = sqrt(1/sum(1/ppstd.^2))
  PP = zeros(max(nh_all),npp);
  one_over_var = zeros(max(nh_all),npp);
  for i_set = 1:numel(nh_all)
    for i_z = 1:nh_all(i_set)       
      PP(i_z,:) = PP(i_z,:) + pp{i_set}(i_z,:)./ppstd{i_set}(i_z,:).^2;
      one_over_var(i_z,:) = one_over_var(i_z,:) + 1./ppstd{i_set}(i_z,:).^2;
    end
  end
  pp = PP./one_over_var;
  ppstd = sqrt(1./one_over_var);
end
% smoother variances...
if lower(tres)=='best'
  ppstd2 = ppstd; 
  for ppind=1:npp
    ppstd2(:,ppind) = std(pp(:,max(1,ppind-5):min(npp,ppind+5))');
  end
  ppstd = ppstd*2/3 + ppstd2*1/3;
end

end
