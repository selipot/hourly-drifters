function [lonlat,varargout] = LatLonLocalWess(ti,loni,lati,p,b,varargin)
% LATLONLOCALWESS is the LOcally WEighted Scatter plot Smoothing estimator
% with a variable bandwidth , See Fans and Gijbels page 24, adapted to 
% observations of time, longitude and latitude.
% [LON,LAT] = LATLONLOCALWESS(TI,LONI,LATI,P,B) implements the local 
% polynomial estimator of order P with a variable bandwith of B points. 
% LON and LAT are column vectors of the same length as TI, LONI and LATI.
% The algorithm iterates 3 times.
    
% Shane Elipot, 2014, version 0

% sort the input?
[ti,I] = sort(ti);
loni = loni(I);
lati = lati(I);

ti = ti(:);
loni = loni(:);
lati = lati(:);

% initialize the output to the input by default
%lon = NaN*ones(size(ti));
%lat = lon;
%lon = loni;
%lat = lati;

% how many optional output arguments?
nout = nargout - 1;

% initialize the output
lonlat = NaN*(1+1i)*ones(size(ti(:),1),p+1);
r = 0*loni;

if nout == 1
    varlonlat = NaN*lonlat;
elseif nout == 2
    varlonlat = NaN*lonlat;
    cd = NaN*ones(size(ti(:),1),2); % return the condition number of the matrix inversion   
end

[n1,~] = size(ti);

if isempty(varargin)
    hcutoff = +Inf;
else
    hcutoff = varargin{1};
end

d = ones(size(ti));

N = 3;% number of iterations

pp = p+2;

for m = 1:N
    for k = 1:length(ti)
        c = 0;
        
        b2 = b;
        
       % initialize q to no data point
       % initialize condition of matrix 
        
        g = 0;
        q = [];
        
        while (length(q) < (pp)) || (g < 1) % while ... 
            
            h = 0.5*(ti(min(n1,k+b2))-ti(max(1,k-b2)))+0*10/(24*60);
            w = kernelHH(ti-ti(k),h);
            w = w.*d;
            
            % need to make sure the problem is not underdetermined: if it is, no estimation possible, or lower the order?
            q = find(w~=0);
            
            %disp([h k b2 length(q) c]);

            w = w(q);
            ti2 = ti(q);
            loni2 = loni(q);
            lati2 = lati(q);
            
            P = length(q);
            
            [xi2,yi2] = latlon2xy(lati2,loni2,lati(k),loni(k));
            
            X = zeros(length(ti2),p+1);
            z = ti2-ti(k);
            z = z(:);
            for j = 0:p
                X(:,j+1) = z.^j;
            end
            W = diag(w,0);
            Y = [xi2(:)+1i*yi2(:)];
            R = transpose(X)*W*X;
            cR = cond(R);
            %disp(['cond(R)=' num2str(cR)]);
            if cR < 1/eps
                
                g = 1;
                
                A = R\transpose(X)*W;
                beta = A*Y;
                
                if m == N && nout >= 1
                    % if only 2 points are used, error should be close to
                    % zero, not Infinity when divided by size(z,1)-2 = 0
                    %e1 = real(X*beta - Y)'*real(X*beta - Y)/(max(size(z,1)-2,2));
                    %e2 = imag(X*beta - Y)'*imag(X*beta - Y)/(max(size(z,1)-2,2));
                    % use a GPS accuracy of 20/10^5 of degrees of latitude = 22 m,
                    % as diagnosed from significant digits
                    e1 = 22*10^-3;
                    e2 = 22*10^-3;
                    % there might be an error here: see eq. 2.112 of
                    % Wunsch's book; should be
                    %C1 = R\X'*W*(e1*eye(length(q)))*W*X/R;
                    %C2 = R\X'*W*(e2*eye(length(q)))*W*X/R;                    
                    C1 = R\X'*(e1*eye(length(q)))*W*X/R;
                    C2 = R\X'*(e2*eye(length(q)))*W*X/R;
%                     if C1>mean(xi2.^2)
%                         C1 = NaN;
%                         C2 = NaN;
%                     end
                    % uncertainty in meters?
                    % varlonlat(k,1) = 1000*(C1(1,1).^0.5+1i*C2(1,1).^0.5);
                    % if k == 1879; disp(['step' num2str(k)]);disp([e1 e2]);disp([X]);disp(W);end
                end
                
                [lat,lon] = xy2latlon(real(beta(1)),imag(beta(1)),lati(k),loni(k));
                
                if m == N && nout >= 1
                    % factor for 95% confidence interval
                    fac = tinv(1-0.025,P-1);
                    [latu,lonu] = xy2latlon(real(beta(1))+fac*C1(1,1).^0.5,imag(beta(1))+fac*C2(1,1).^0.5,lati(k),loni(k));
                    % error in longitude latitude following a
                    % great circle?
                    varlonlat(k,1) = abs(lonu-lon) + 1i*abs(latu-lat);
                    %if k == 1879; disp([lonu latu]);end
                end
                
                if m<N
                    r(k) = spheredist(lat,lon,lati(k),loni(k));% the residual is xy - 0+1i*0 the origin of tangent plane
                end
                
                lonlat(k,1) = lon+1i*lat;
                
                if nout == 2
                    cd(k,1) = cR;
                    cd(k,2) = length(q);
                end
                
                if p > 0
                    % calculate derivative of position, in m/s^p
                    lonlat(k,2:p+1) = 1000*(86400.^-(1:p)).*factorial(1:p).*beta(2:p+1).';
                    % calculate variance of estimates
                    if  m == N && nout >= 1
                        foo1 = diag(C1).^0.5;
                        foo2 = diag(C2).^0.5;
                        % uncertainties in meters per second
                        varlonlat(k,2:p+1) = fac*(1000*(86400.^-(1:p)).*factorial(1:p).*(foo1(2:p+1)+1i*foo2(2:p+1)).');
                    end
                end
                
            end
            
            co = (length(q) < pp) + 2*(g<1);
            
            switch co
                
                case 1
                    %disp(['co=' num2str(co) ';k=' num2str(k) ';n=' num2str(m) '; increasing b to ' num2str(b2) ' because of undetermined system q=' num2str(length(q))]);
                    b2 = b2+1;
                    c = c+1;
                case 2
                    %disp(['co=' num2str(co) ';k=' num2str(k) ';n=' num2str(m) '; increasing b to ' num2str(b2) ' because of cond(R)=' num2str(cR)]);
                    b2 = b2+1;
                    c = c+1;
                    
                case 3
                    %disp(['co=' num2str(co) ';k=' num2str(k) ';n=' num2str(m) '; increasing b to ' num2str(b2) ' because of cond(R)=' num2str(cR) ' and undetermined system q='  num2str(length(q))]);
                    b2 = b2+1;
                    c = c+1;
            end
            
            % limit on possible infinite while statement
            if c > 10

                break
                
            end
            
        end
        
       
        
    end
    
    if m<N
        M = nanmedian(abs(r));
        d = kernelB(r/(6*M));
        % there is an issue of NaN here? how to deal with it?
        qnan = isnan(d);
        d(qnan) = 1;
    end
   
end

if nout == 1
    varargout{1} = varlonlat;
elseif nout == 2
    varargout{1} = varlonlat;
    varargout{2} = cd;
end

%return

figure
plot(loni,lati,'.k');
hold on
plot(360+lonlat(:,1),'b');
xx = 360+[1;1]*[real(lonlat(:,1))]';
yy = [imag(lonlat(:,1)+varlonlat(:,1)) imag(lonlat(:,1)-varlonlat(:,1))].';
line(xx,yy);
yy = [1;1]*[imag(lonlat(:,1))]';
xx = 360+[real(lonlat(:,1)+varlonlat(:,1)) real(lonlat(:,1)-varlonlat(:,1))].';
line(xx,yy);

function K = kernelHH(t,h)
% the kernel is zero outside of the normalized bandwith 1 by construction
    K = (1-abs(t/h).^3);
    qn = K<0;
    K(qn) = 0;
    K = (1/h)*(70/81)*K.^3;
    
function K = kernelB(t)
% the kernel is zero outside of the normalized bandwith 1 by construction
    K = (1-abs(t).^2).^2;
    qn = abs(t)>=1;
    K(qn) = 0;
    


