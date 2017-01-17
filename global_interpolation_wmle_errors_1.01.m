%Shane Elipot, University of Miami, RSMAS, selipot@rsmas.miami.edu, January 2017
% code to apply the WMLE interpolation to Argos drifters
% and calculate errors by bootstrapping
% see Elipot et al, http://dx.doi.org/10.1002/2016JC011716

% this code may not run as it is and is put on this repository for openness on our methods. 

% this script requires to have the Jlab toolbox for Matlab that includes several routines to deal with cell variables
% and geographical data

% run in Matlab with parallel computing toolbox as an example as a batch
% job
% jj = batch('global_interpolation_wmle_errors_update','Pool',10,'CaptureDiary',true);

% this code assume that you are loading the data as Matlab cell arrays 
% and these data are for Argos-tracked drifters
load('drifterCELL.mat'],'IDunique'); % array of AOML unique IDs
load('drifterCELL.mat'],'num'); % cell array of times as returned by datenum.m
load('drifterCELL.mat'],'lon','lat'); % cell arrays of geographical corrdinates
load('drifterCELL.mat'],'qual'); % cell array of Argos quality index
load('drifterCELL.mat'],'Exp'); % cell array of experiment number
load('drifterCELL.mat'],'Droff'); % date of drogue loss, NaN if still on

% dirout
dirout = 'individual_1.01/';

M = size(IDunique,1);

tr = 1/24; % interpolating interval in days

TUTC = cell(M,1);
LON = cell(M,1);
LAT = cell(M,1);
UV = cell(M,1);
H = cell(M,1);
ELONLAT = cell(M,1);
FVAL = cell(M,1);
DRG = cell(M,1);

tic;

parfor m = 1:M
    
    try
        fprintf('m=%d ID=%d\n',[(m) (double(IDunique(m)))]);
        
        [lona,lata,ta,qa] = argos_reduce(unwrap(lon{m}*pi/180)*180/pi,lat{m},num{m},qual{m});
        
        % 3 points are necessary
        if length(ta)>=3
            
            % replace unknown or 0 quality with quality -1
            qa(qa~=1 & qa~=2 & qa~=3) = -1;
            
            % define interpolation times
            TUTC{m} = [ceil(ta(1)*24)/24:tr:floor(ta(end)*24)/24]';
                   
            % calculate size of gaps first and estimate only for gaps < 12h
            ga = NaN*TUTC{m};
            for k = 1:length(ga)
                qf = find(ta>=TUTC{m}(k));
                qb = find(ta<=TUTC{m}(k));
                ga(k) = ta(qf(1))-ta(qb(end));
            end
            
            % initialize output
            lonlat = Inf*(1+1i)*ones(length(TUTC{m}),2);
            elonlat = lonlat;
            hout = Inf*(ones(length(TUTC{m}),2*2));
            fval = Inf*ones(length(TUTC{m}),1);
            qq = ga <= 0.5;
            % the unwrap should not be necessary but let's keep it anyway
            [dum1,dum2,dum3,dum4] = LatLonLikelihoodHessian(ta,unwrap(lona*pi/180)*180/pi,lata,TUTC{m}(qq),1,2,qa,2,'boot');
            %[lonlat,hout,elonlat,fval]
            lonlat(qq,:) = dum1;
            hout(qq,:) = dum2;
            elonlat(qq,:) = dum3;
            fval(qq) = dum4;
            
            % position
            % reestablish longitudes between 0 and 360
            % the mod operation unfortunatelly here convert the Inf to NaN
            LON{m} = mod(real(lonlat(:,1)),360);
            LAT{m} = imag(lonlat(:,1));
            % errors in lat lon and their derivatives
            ELONLAT{m} = elonlat;
            % velocity
            UV{m} = lonlat(:,2);
            
            NN = sum(~isnan(hout.')).';
            grms = sqrt(nanmean(hout'.^2))';
            H{m} = [ga grms];
            % likelihood value
            FVAL{m} = fval;
            
            % drogue status
            if isnan(Droff(m))
                DRG{m} = ones(size(TUTC{m}));
            else
                DRG{m} = double(TUTC{m}<=Droff(m));
            end
           
            % factor for calculating confidence intervals
            cifac = 3.1824./sqrt(NN);
            % save a text version
            fileout = [dirout 'd' num2str(IDunique(m))];
            fid = fopen(fileout,'w');
            dmy = datevec(TUTC{m});
            xout = [double(IDunique(m))*ones(size(lonlat,1),1) dmy(:,1:4) real(lonlat(:,1)) ...
                imag(lonlat(:,1)) round(10^5*cifac.*real(elonlat(:,1))) round(10^5*cifac.*imag(elonlat(:,1))) ...
                real(lonlat(:,2)) imag(lonlat(:,2))  real(cifac.*elonlat(:,2)) imag(cifac.*elonlat(:,2)) ga grms DRG{m}];
            %fprintf(fid,{ff},xout.');
            fprintf(fid,'%8d %4d %3d %3d %3d %9.5f %9.5f %8d %8d %8.4f %8.4f %8.4f %8.4f %6.4f %6.4f %3d\n',xout.');
            fclose(fid);
            
        end
    catch
        fprintf('m=%d ID=%d failed\n',[(m) (double(IDunique(m)))]);
        
    end
    
end
toc;

% A number of post-processing steps ensue to remove interpolating gaps that are too large and organize the data.

