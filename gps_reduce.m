function [lona,lata,ta] = gps_reduce(loni,lati,ti)
% code to reduce GPS data before interpolating
%Shane Elipot, University of Miami, RSMAS, selipot@rsmas.miami.edu, January 2017

% first remove exactly repeated data
X = unique([loni lati ti],'rows','stable');
loni = X(:,1);
lati = X(:,2);
tai = X(:,3);

% then remove repeated times and average positions

% preallocate
ta_tmp = Inf*tai;
lona_tmp = tai;
lata_tmp = tai;

n = 0;
while max(size(tai))>1
   
   dt = tai - tai(1);
   q = find(dt == 0);
   n = n+1;
   ta_tmp(n,1) = tai(1);
   lona_tmp(n,1) = mean(loni(q));
   lata_tmp(n,1) = mean(lati(q));
    
   tai(q) = [];
   loni(q) = [];
   lati(q) = [];

end

% if there is one last point alone
if length(tai) == 1
    
    n = n+1;
    ta_tmp(n,1) = tai(1);
    lona_tmp(n,1) = loni(1);
    lata_tmp(n,1) = lati(1);
    
    tai(1) = [];
    loni(1) = [];
    lati(1) = [];
end


% shrink
q = ta_tmp==Inf;
ta_tmp(q) = [];
lona_tmp(q) = [];
lata_tmp(q) = [];

% now remove static points

% preallocate
ta = Inf*ta_tmp;
lona = ta;
lata = ta;

n = 0;
while max(size(ta_tmp))>1
    
    dxa = spheredist(lata_tmp,lona_tmp);% dxa(1)=0 by construction
    q = find(dxa==0);% q is always at least of size 1
    n = n + 1;
    ta(n,1) = 0.5*(ta_tmp(q(1))+ta_tmp(q(end)));
    lona(n,1) = lona_tmp(q(1));
    lata(n,1) = lata_tmp(q(1));
    
    ta_tmp(q) = [];
    lona_tmp(q) = [];
    lata_tmp(q) = [];

end

% if there is one last point alone
if length(ta_tmp) == 1
    
    n = n + 1;
    ta(n,1) = ta_tmp(1);
    lona(n,1) = lona_tmp(1);
    lata(n,1) = lata_tmp(1);
    
    ta_tmp(1) = [];
    lona_tmp(1) = [];
    lata_tmp(1) = [];

end

% shrink
q = ta==Inf;
ta(q) = [];
lona(q) = [];
lata(q) = [];

return

figure;
subplot(2,1,1);
hold on;plot(ti,loni,'.');
plot(ta,lona,'o');
axis tight
subplot(2,1,2);
hold on;plot(ti,lati,'.');
plot(ta,lata,'o');
axis tight
