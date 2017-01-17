function [lonlat,hout,varargout] = LatLonLikelihoodHessian(ti,loni,lati,t,p,b,cin,d,varargin)
%  [LONLAT,HOUT,VARARGOUT] = LATLONLIKELIHOODHESSIAN(TI,LONI,LATI,T,P,B,CIN,D,VARARGIN)
% P is the order of polynomial
% B is the number of points to use or before and after the times in t
% CIN is the error class
% D is a number indicating the error distribution pdf model
% 0: Gaussian PDF, 2: t-location-scale PDF
% -1 is a possible option for D for  returning the unweighted least square solution  
% option is for boostraping or calculating Hessian: 'boot' or 'hess'
% Shane Elipot, University of Miami, RSMAS, selipot@rsmas.miami.edu, January 2017

% radius of the Earth
R = 6371000;

% sort the input; should already be sorted
[ti,I] = sort(ti);
loni = loni(I);
lati = lati(I);
ti = ti(:);
loni = loni(:); 
lati = lati(:);

ni = size(ti,1);

% initialize the output
lonlat = NaN*ones(length(t),p+1)*(1+1i);

% output parameter h
% returns the temporal distance to all points used

hout =  NaN*ones(length(t),2*b);

% how many optional input argument
nin = nargin - 8;
if nin == 1
    option = varargin{1};
end

% how many optional output arguments?
nout = nargout - 2;
if nout == 1
    switch d
        case -1
            ye = NaN*ones(length(t),p+1); % output error
        case {0,2}
            ye = NaN*ones(length(t),2*(p+1)); % output error
    end
    
    elseif nout == 2
    switch d
        case -1
            ye = NaN*ones(length(t),p+1); % output error
        case {0,2}
            ye = NaN*ones(length(t),(p+1)); % output error
    end
    cd = NaN*ones(length(t),1); % something about the minimization?
end

% form of the errors as a input? 
if length(cin) == 1
    ci = cin*ones(size(ti));
elseif (length(cin) == ni)
    % sort and make column vector
    ci = cin(I);
    ci = ci(:);
end

% options for optimization
options1 = optimoptions('fminunc','Display','off','Algorithm','trust-region','GradObj','on');%
% options2 = optimoptions('fminunc','Display','off','Algorithm','quasi-newton','GradObj','off');
% convert class to unbiased Gaussian errors for first guess
%  Numbers estimated from comparison of Argos and GPS data
% real part is for zonal, imaginary part for meridional
% resolution of 1 m implies 1./(111.12*1000)=8.99*10^-6 latitude resolution; round up
% to 9*10^-6 or 10^-5;

ei = NaN*ci*(1+1i);
ei(ci==1) = (641+1i*473)*10^-5;
ei(ci==2) = (390+1i*341)*10^-5;
ei(ci==3) = (294+1i*280)*10^-5;
ei(ci==-1) = (445+1i*364)*10^-5;

switch d
    
    case {-1,0}
        % joint normal distribution
          % longitude
        % class 1
        mux(ci==1,1) = 10*10^-5;
        sigmax(ci==1,1) = 641*10^-5;
        % class 2
        mux(ci==2,1) = -9*10^-5;
        sigmax(ci==2,1) = 390*10^-5;
        % class 3
        mux(ci==3,1) = -5*10^-5;
        sigmax(ci==3,1) = 294*10^-5;
        
        % latitude
        % class 1
        muy(ci==1,1) = 6*10^-5;
        sigmay(ci==1,1) = 473*10^-5;
        % class 2
        muy(ci==2,1) = 55*10^-5;
        sigmay(ci==2,1) = 341*10^-5;
        % class 3
        muy(ci==3,1) = 62*10^-5;
        sigmay(ci==3,1) = 280*10^-5;
        
    case 2
        
        % joint t-distribution
        % longitude
        % class 1
        mux(ci==1,1) = 9.28747e-05;
        sigmax(ci==1,1) = 0.00448337;
        nux(ci==1,1) = 3.75179;
        % class 2
        mux(ci==2,1) = -8.61962e-05;
        sigmax(ci==2,1) = 0.00297111;
        nux(ci==2,1) = 4.70047;
        % class 3
        mux(ci==3,1) = -4.52899e-05;
        sigmax(ci==3,1) = 0.00235822;
        nux(ci==3,1) = 5.54998;
        % unknown class
        mux(ci==-1,1) = -3.81203e-05;
        sigmax(ci==-1,1) = 0.00294147;
        nux(ci==-1,1) = 3.38047;
        
        % latitude
        % class 1
        muy(ci==1,1) = 0.000144521;
        sigmay(ci==1,1) = 0.00357287;
        nuy(ci==1,1) = 4.44095;
        % class 2
        muy(ci==2,1) = 0.000540441;
        sigmay(ci==2,1) = 0.00268321;
        nuy(ci==2,1) = 5.09186;
        % class 3
        muy(ci==3,1) = 0.000603532;
        sigmay(ci==3,1) = 0.00221004;
        nuy(ci==3,1) = 5.11416;
        % unknown class
        muy(ci==-1,1) = 0.000483519;
        sigmay(ci==-1,1) = 0.0026613;
        nuy(ci==-1,1) = 4.07786;
end

% rate of rotation of the Earth
omega = 7.292e-5;

for k = 1:length(t)
    
    % find the b closest points in time
    dt = ti-t(k);

    if 1

        %disp(['Selecting ' num2str(b) ' points on each side']);
%         % set number on each side
%         qp = find(ti>=t(k));
%         qn = find(ti<=t(k));
%         % deal with end points
%         if length(qn) == 1 && length(qp) >= 3            
%             qnp = [qn(1) qp(1):1:min(qp(1)+b-1,ni)];
%         elseif length(qp) == 1 && length(qn) >=3
%             qnp = [max(1,qn(end)-b+1):1:qn(end) qp(end)];
%         else
%             qnp = unique([max(1,qn(end)-b+1):1:qn(end) qp(1):1:min(qp(1)+b-1,ni)]);
%         end
%         q = qnp;
%         
          qf = find(dt>=0);
          qb = find(dt<=0);

          if length(qf)==1 && length(qb)>=b+1
              qnp = sort(unique([qb(end-b:end) ; qf(1)]));
          elseif length(qb)==1 && length(qf)>=b+1
              qnp = sort(unique([qb(1) ; qf(1:1+b)]));
          else
              qnp = sort(unique([qb(end-(b-1):end) ; qf(1:1+b-1)]));
          end
          q = qnp;
          
    else
        %disp(['Selecting the closest ' num2str(b) ' points']);
       [~,q] = sort(abs(dt),'ascend');
        if all(dt(q(1:b))>0)
            qn = find(dt<=0);
            q = [q(1:max(1,b-1)) ; qn(end)];
        elseif all(dt(q(1:b))<0)
            qp = find(dt>=0);
            q = [q(1:max(1,b-1)) ; qp(1)];
        else
            q = q(1:b);
        end
    end
    
    % data to use
    lati2 = lati(q);
    loni2 = loni(q);
    ti2 = ti(q);
    ei2 = ei(q);
    ci2 = ci(q);
    dt2 = dt(q); 
    %h2 = max(ti2)-t(k);
    %h1 = t(k)-min(ti2);
    foo = NaN*ones(1,2*b);
    foo(1:length(ti2)) = t(k)-ti2';
    hout(k,:) = foo;
    
    switch d
        
      case {-1,0}
            mu = mux(q)+1i*muy(q);
            sigma = sigmax(q)+1i*sigmay(q);
        case 2
            mu = mux(q)+1i*muy(q);
            sigma = sigmax(q)+1i*sigmay(q);
            nu = nux(q)+1i*nuy(q);
    end
    
    % need to make sure the problem is not underdetermined:
    % if it is, no estimation possible, or lower the order maximum 
    % order of polynomial achievable
    % the order of the polynomial plus one cannot be more than the number
    % of data and no less than 0
    pmax = max(0,min(length(q)-1,p));
    
     % build matrices for first guess, least square
    if length(q)>=pmax+1
        
        X = zeros(length(ti2),pmax+1);
        z = ti2-t(k);
        z = z(:);
        for j = 0:pmax
            X(:,j+1) = z.^j; 
        end
        Y = loni2(:)+1i*lati2(:);

        % first guess is unweighted least square
        w00 = ones(size(ei2));
        W00 = diag(w00,0);
        A = transpose(X)*W00*X;
        beta0 = A\transpose(X)*W00*Y;
        

        if 1 % inverse distance weighting
            
            s0 = 20/(24*60);
            s = s0+abs(dt2(:));
            c = s.^-1;
            w1 = c*length(c)/sum(c);
            w2 = w1;
           
        else %polynomial weighting
            
            h2 = max(ti2,[],1);
            h1 = min(ti2,[],1);
            hh = h2(end)-h1(1);
            c = weight(abs(dt2(:)),hh,1);
            w1 = c*length(c)/sum(c);
            w2 = w1;
            
        end
        
        switch d
          
          case -1
          
            beta = beta0;
            W = ones(size(w1(:)))+1i*ones(size(w2(:)));
            [fval,~] = likelihood0([real(beta);imag(beta)],Y,X,W,mu,sigma);
            
            if exist('ye','var')
                % calculate the variance for the least square solution?
                RR1 = diag(real(ei2)).^2;
                RR2 = diag(imag(ei2)).^2;
                Cxx1 = A*transpose(X)*RR1*X/A;
                Cxx2 = A*transpose(X)*RR2*X/A;
                betaV = diag(Cxx1).^0.5+1i*diag(Cxx2).^0.5;
            end
            
          case 0
            
            % weighted likelihood with Gaussian errors
            W = w1(:)+1i*w2(:);
            %unconstrained
            [beta,fval,~,~,~,hessian] = fminunc(@(theta) likelihood0(theta,Y,X,W,mu,sigma),[real(beta0);imag(beta0)],options1);
            %[val,grd] = likelihood0(beta,Y,X,W,mu,sigma);
            %[beta,fval,~,~,~,hessian] = fminunc(@(theta) likelihood0(theta,Y,X,W,mu,sigma),[real(beta0);imag(beta0)],options2);
            %[val,grd] = likelihood0(beta,Y,X,W,mu,sigma);
            %H = HessN(X,W,sigma);
            
             if strcmp(option,'boot')
                 % bootstrap
                 betaboot = cell(length(ti2)-1,1);
                 JJ = nchoosek(1:length(ti2),length(ti2)-1);
                 for jj = 1:size(JJ,1)
                     dum1 = real(W(JJ(jj,:)));
                     dum1 = dum1*length(dum1)/sum(dum1);
                     dum2 = imag(W(JJ(jj,:)));
                     dum2 = dum2*length(dum2)/sum(dum2);
                     WJ = dum1+1i*dum2;
                     [betaboot{jj},~,~,~,~,~] = fminunc(@(theta) likelihood0(theta,Y(JJ(jj,:)),X(JJ(jj,:),:),...
                         WJ,mu(JJ(jj,:)),sigma(JJ(jj,:))),[real(beta0);imag(beta0)],options1);
                 end
                 H = 1./var(cell2mat(betaboot'),0,2)';
             else
                 H = diag(hessian);
             end
            
            %disp(H);disp(diag(hessian));
            foo = (H(:).^-1/(1)).^0.5;
            %foo2 = diag(hessian).^-0.5;
            %disp(foo);disp(foo2);
            betaV(1,:) = foo(1:pmax+1,:)+1i*foo(pmax+2:2*(pmax+1),:);
            %betaV(2,:) = foo2(1:pmax+1,:)+1i*foo2(pmax+2:2*(pmax+1),:);

            beta = beta(1:size(X,2))+1i*beta(size(X,2)+1:2*size(X,2));

          case 2
          
            W = w1(:)+1i*w2(:);
             % option 1 grad is provided
             if 1
                 [beta,fval,~,~,~,hessian] = fminunc(@(theta) likelihood2(theta,Y,X,W,mu,sigma,nu),[real(beta0);imag(beta0)],options1);
             else                  % for global interpolation and calculating uncertainties only
                 beta = [real(beta0);imag(beta0)];
                 fval = 0;
             end
             
             if strcmp(option,'boot')
                 % bootstrap
                 betaboot = cell(length(ti2)-1,1);
                 JJ = nchoosek(1:length(ti2),length(ti2)-1);
                 for jj = 1:size(JJ,1)
                     dum1 = real(W(JJ(jj,:)));
                     dum1 = dum1*length(dum1)/sum(dum1);
                     dum2 = imag(W(JJ(jj,:)));
                     dum2 = dum2*length(dum2)/sum(dum2);
                     WJ = dum1+1i*dum2;
                     [betaboot{jj},~,~,~,~,~] = fminunc(@(theta) likelihood2(theta,Y(JJ(jj,:)),X(JJ(jj,:),:),...
                         WJ,mu(JJ(jj,:)),sigma(JJ(jj,:)),nu(JJ(jj,:))),[real(beta0);imag(beta0)],options1);
                  end
                 H = 1./var(cell2mat(betaboot'),0,2)';
             else
                 H = diag(hessian);
             end
             
             %disp('option2')
             %[beta,fval,~,~,~,hessian] = fminunc(@(theta) likelihood2(theta,Y,X,W,mu,sigma,nu),[real(beta0);imag(beta0)],options2);
             
             %disp('fminsearch');
             %[beta2,fval2] = fminsearch(@(theta) likelihood2(theta,Y,X,W,mu,sigma,nu),[real(beta0);imag(beta0)]);

            % alternate way: calculate exact Hessian; only for order 1 of
            % polynomial model
            %H = HessT(beta,Y,X,W,mu,sigma,nu);

            foo = (H(:).^-1/(1)).^0.5;
            %foo2 = diag(hessian).^-0.5;
            
            betaV(1,:) = foo(1:pmax+1,:)+1i*foo(pmax+2:2*(pmax+1),:);
            %betaV(2,:) = foo2(1:pmax+1,:)+1i*foo2(pmax+2:2*(pmax+1),:);

            beta = beta(1:size(X,2))+1i*beta(size(X,2)+1:2*size(X,2));

        end
        
        % return results
        lonlat(k,:) = beta; 
        cd(k,1) = fval;
        
        % convert derivative of longitude and latitudes for speed and
        % acceleration
        if pmax+1 > 1
            lonlat(k,2) = (24*60*60)^(-1)*(R*cosd(imag(beta(1)))*(2*pi/360)*real(beta(2)) ...
                + 1i*R*(2*pi/360)*imag(beta(2)));
        elseif pmax+1 > 2
            lonlat(k,2) = (24*60*60)^(-1)*(R*cosd(imag(beta(1)))*(2*pi/360)*real(beta(2)) ...
                + 1i*R*(2*pi/360)*imag(beta(2)));
            lonlat(k,3) = (24*60*60)^(-2)*(R*cosd(imag(beta(1)))*(2*pi/360)*real(beta(3)) ...
                - R*((2*pi/360)*imag(beta(2)))^2*sind(imag(beta(1)))+ 1i*R*(2*pi/360)*imag(beta(3)));
        end
        
         if exist('ye','var')
             
             ye(k,:) = betaV(:).';
         
         end
        
    end
    
end

if nout >= 1
    
    varargout{1} = ye;

end

if nout >= 2

    varargout{2} = cd;

end

    function [val,grd] = likelihood0(theta,Y,X,W,mu,sigma)
        
        % likelihood that also returns the gradient
        
        Zr = bsxfun(@times,real(Y) - real(mu) - X*theta(1:size(X,2)),1./real(sigma));
        Zi = bsxfun(@times,imag(Y) - imag(mu) - X*theta(size(X,2)+1:2*size(X,2)),1./imag(sigma));
        val = sum(0.5*(real(W).*Zr.^2+ imag(W).*Zi.^2));
        
        if nargout > 1
            % assumes it is order p = 1
            grd(1,1) =  -sum(real(W).*bsxfun(@times,real(Y) - real(mu) - X*theta(1:size(X,2)),1./(real(sigma)).^2));
            grd(2,1) = -sum(real(W).*bsxfun(@times,real(Y) - real(mu) - X*theta(1:size(X,2)),1./(real(sigma)).^2).*X(:,2));
            grd(3,1) =  -sum(imag(W).*bsxfun(@times,imag(Y) - imag(mu) - X*theta(size(X,2)+1:2*size(X,2)),1./(imag(sigma)).^2));
            grd(4,1) = -sum(imag(W).*bsxfun(@times,imag(Y) - imag(mu) - X*theta(size(X,2)+1:2*size(X,2)),1./(imag(sigma)).^2).*X(:,2));
        end
        
    function [val,grd] = likelihood2(theta,Y,X,W,mu,sigma,nu)

        Zr = real(Y) - real(mu) - X*theta(1:size(X,2));
        Zi = imag(Y) - imag(mu) - X*theta(size(X,2)+1:2*size(X,2));
        valx = 0.5*(real(nu)+1).*real(W).*log(1+((real(nu).*real(sigma).^2).^-1).*Zr.^2);
        valy = 0.5*(imag(nu)+1).*imag(W).*log(1+((imag(nu).*imag(sigma).^2).^-1).*Zi.^2);
        val = sum(valx+valy,1);

        if nargout > 1
            % assume it is order 1
            Zr = real(Y) - real(mu) - X*theta(1:size(X,2));
            n1 = (real(nu)+1).*real(W).*Zr.*(real(nu).*real(sigma).^2).^-1;
            d1 = 1+((real(nu).*real(sigma).^2).^-1).*Zr.^2;

            Zi = imag(Y) - imag(mu) - X*theta(size(X,2)+1:2*size(X,2));
            n2 = (imag(nu)+1).*imag(W).*Zi.*(imag(nu).*imag(sigma).^2).^-1;
            d2 = 1+((imag(nu).*imag(sigma).^2).^-1).*Zi.^2;
            
            grd = [-sum(n1./d1);...
                -sum(X(:,2).*n1./d1);...
                -sum(n2./d2);...
               -sum(X(:,2).*n2./d2)];
        
        end
                
    function H = HessN(X,W,sigma)
    % works only for order 1
    H = NaN*ones(4,1);
    term1 = real(sigma).^-2.*real(W);
    term2 = imag(sigma).^-2.*imag(W);
    H(1,1) = sum(term1,1);
    H(3,1) = sum(term2,1);
    term3 =  term1.*X(:,2).^2;
    term4 =  term2.*X(:,2).^2;
    H(2,1) = sum(term3,1);
    H(4,1) = sum(term4,1);
        
    function H = HessT(beta,Y,X,W,mu,sigma,nu)
        % works for order1 only
        
        H = NaN*ones(4,1);
        
        Zr = real(Y) - real(mu) - X*beta(1:size(X,2));
        
        Zi = imag(Y) - imag(mu) - X*beta(size(X,2)+1:2*size(X,2));
        
        term1 = ((real(nu)+1).*real(W)).*(real(nu).*real(sigma).^2+Zr.^2).^-1 -...
            2*(real(nu)+1).*real(W).*Zr.^2.*(real(nu).^2.*real(sigma).^4.*(1+(real(nu).*real(sigma).^2).^-1.*Zr.^2).^2).^-1;
        
        term2 = ((imag(nu)+1).*imag(W)).*(imag(nu).*imag(sigma).^2+Zi.^2).^-1 -...
            2*(imag(nu)+1).*imag(W).*Zi.^2.*(imag(nu).^2.*imag(sigma).^4.*(1+(imag(nu).*imag(sigma).^2).^-1.*Zi.^2).^2).^-1;
        
        term3 = ((real(nu)+1).*real(W).*X(:,2).^2).*(real(nu).*real(sigma).^2+Zr.^2).^-1 - ...
            2*((real(nu)+1).*real(W).*Zr.^2.*X(:,2).^2).*(real(nu).^2.*real(sigma).^4.*(1+(Zr.^2.*(real(nu).*real(sigma).^2).^-1)).^2).^-1;
        
        term4 = ((imag(nu)+1).*imag(W).*X(:,2).^2).*(imag(nu).*imag(sigma).^2+Zi.^2).^-1 - ...
            2*((imag(nu)+1).*imag(W).*Zi.^2.*X(:,2).^2).*(imag(nu).^2.*imag(sigma).^4.*(1+(Zi.^2.*(imag(nu).*imag(sigma).^2).^-1)).^2).^-1;
        
        H(1,1) = sum(term1,1);
        H(2,1) = sum(term3,1);
        H(3,1) = sum(term2,1);
        H(4,1) = sum(term4,1);

function W = weight(t,h,g)

    %g = 1;
    W = kernelB(t/h,g)/h;

function K = kernelB(t,g) % beta family kernel for g = 0,1,2, ...
    
% the kernel is zero outside of the normalized bandwith 1 by construction
    K = (1-t.^2);
    qn = K<0;
    K(qn) = 0;
    K = 0.75*K.^g;
    
        

