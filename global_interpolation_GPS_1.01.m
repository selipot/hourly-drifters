%Shane Elipot, University of Miami, RSMAS, selipot@rsmas.miami.edu, January 2017
% global interpolation GPS: code to apply the LOWESS filter + linear
% interpolation on the global GPS drifter dataset
% see Elipot et al, http://dx.doi.org/10.1002/2016JC011716
% as an example run as batch job in Matlab with the parallel toolbox and 10
% workers
% run with jj = batch('global_interpolation_GPS_update','Pool',10,'CaptureDiary',true);

% this code may not run as it is and is put on this repository for openness on our methods. 

% this script requires to have the Jlab toolbox for Matlab that includes several routines to deal with cell variables and geographical data

% to be applied on the new GPS drifters and the ones continuing from
% previous iteration

% load info for the list of IDs and respective drogue off dates
% see deployment files and directory files
load('drifterCELL.mat','Droff','IDunique');

% load GPS drifter data 
% Here the drifter data were arranged in cell arrays of column vector
load('drifterCELL.mat','id','num','lon','lat','qual');
% num is a cell array of times as returned by Matlab datenum.m
% lon, lat are cell arrays of longitude and latitude
% id is a cell array of id numbers, repeated for each observations of a single trajectory
% qual is a cell array of GPS drifter quality indices, either all 0 or >9

% if all quality indices are zero this means they are all ok;
% otherwise, they are not

qual9 = cell(length(qual),1);
for k = 1:length(qual)
    if all(qual{k} == 0)
        qual09{k} = true(size(qual{k}));
    else
        qual09{k} =  qual{k}>=9 ;
    end
end

% reduce, i.e keep only all zeros or quality 9 index only
[lon2,lat2,num2,qual2,id2] = cellindex(lon,lat,num,qual,id,qual09);

% make sure the longitudes are between 0 and 360
for k = 1:length(lon2)
    q = find(lon2{k}<0);
    if ~isempty(q)
        lon2{k}(q) = lon2{k}(q)+360;
    end
end

M = length(id2); %
tr = 1/24; % interpolation interval = 1 hour
b = 2;% plus or minus b points to define the bandwith
p = 1;% fit linear trend function
    
numgs = cell(M,1);
ids = numgs;
longs = numgs;
latgs = numgs;
uvgs = numgs;
gap = numgs;
drg = numgs;
latstd = numgs;
lonstd = numgs;
uvgstd = numgs;

tic;
parfor k = 1:M; %number of drifter
    
    try
        
        disp(['k = ' num2str(k) ' working on drifter #' num2str(id2{k}(1)) ' the time is ' datestr(now)]);
        
        % this preprocessing look for periods of time when the drifter is
        % static; take the middle time
        
        [long,latg,tg] = gps_reduce(lon2{k},lat2{k},num2{k});
        
        %disp([num2str(length(num2{k})-length(tg)) ' GPS points removed out of ' num2str(length(num2{k})) ' in drifter #' num2str(id2{k}(1))]);
        
        % Main interpolation routine ; last argument is a limit on size of gaps
        [longlatgs,cilonlat,~] = LatLonLocalWess(tg,unwrap(long*pi/180)*180/pi,latg,p,b,+Inf);
        % sometimes the variance or uncertainty estimation fails
        % GPS error is assumed isotropic, resulting in isotropic velocity
        % error but not isotropic lon/lat error
        
        % then linearly interpolate to hourly time steps
        numgs{k} = [ceil(tg(1)*24)/24:tr:floor(tg(end)*24)/24]';
         
        % make an id cell
        ids{k} = id2{k}(1)*ones(size(numgs{k}));
        
        [latgs{k},latstd{k}] = piecelinvar(tg,imag(longlatgs(:,1)),numgs{k},imag(cilonlat(:,1)));
        
        % longitude is treated carefully; unwrap result for linear
        % interpolation then recast between 0 and 360
        [dum1,lonstd{k}] = piecelinvar(tg,unwrap(real(longlatgs(:,1))*pi/180)*180/pi,numgs{k},real(cilonlat(:,1)));
        longs{k} = mod(dum1,360);
        
        [u,ciu] = piecelinvar(tg,real(longlatgs(:,2)),numgs{k},real(cilonlat(:,2)));
        [v,civ] = piecelinvar(tg,imag(longlatgs(:,2)),numgs{k},imag(cilonlat(:,2)));
        uvgs{k} = u+1i*v;
        uvgstd{k} = ciu+1i*civ;
        
        % drogue status
        qid = find(IDunique==id2{k}(1));
        % if isnan(Droff) this means that the drogue is still on at the time of latest update
        if isnan(Droff(qid))
            drg{k} = ones(size(numgs{k}));
        else
            drg{k} = double(numgs{k}<=Droff(qid));
        end
        
        % calculate size of gaps
        gap{k} = NaN*ones(size(numgs{k}));
        for m = 1:length(numgs{k})
            qf = find(tg>=numgs{k}(m));
            qb = find(tg<=numgs{k}(m));
            if ~isempty(qf) && ~isempty(qb)
                gap{k}(m) = tg(qf(1))-tg(qb(end));
            end
        end
     
        %fprintf(['GPS drifter #' num2str(id2(k)) ' done']);
        fprintf('segment k=%d for ID=%d done\n',[k (double(id2{k}(1)))]);

    catch
        
        fprintf('segment k=%d for ID=%d failed\n',[(k) (double(id2{k}(1)))]);
        
    end
    
end
  
toc;

drifterGPS.id = ids;
drifterGPS.num = numgs;
drifterGPS.lon = longs;
drifterGPS.lonstd = lonstd;
drifterGPS.lat = latgs;
drifterGPS.latstd = latstd;
drifterGPS.uv = uvgs;
drifterGPS.uvstd = uvgstd;
drifterGPS.gap = gap;
drifterGPS.drogue = drg;

save('drifterGPS_update.mat','drifterGPS');

% then I conduct a number of diagnostics to remove failed interpolation and interpolation over gaps larger than 12 h etc.

return
