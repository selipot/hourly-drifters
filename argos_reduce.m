function [lona,lata,ta,qa] = argos_reduce(loni,lati,ti,qi);
% Shane Elipot, University of Miami, RSMAS, selipot@rsmas.miami.edu, January 2017
% code to reduce the quality-controlled Argos data in 20 min windows before interpolation
    
        lona_tmp = loni;
        lata_tmp = lati;
        ta_tmp =  ti;
        qa_tmp =  qi;
        
        n = 0;
        while ~isempty(ta_tmp)
            
            q = find(abs(ta_tmp(1)-ta_tmp)<=20/(24*60));

            n = n+1;
            
            qamax = max(qa_tmp(q));
            I = qa_tmp(q) == qamax;
            
            ta(n,1) = mean(ta_tmp(q(I)));
            lona(n,1) = mean(lona_tmp(q(I)));
            lata(n,1) = mean(lata_tmp(q(I)));
            qa(n,1) = qamax;
            
            ta_tmp(q) = [];
            lona_tmp(q) = [];
            lata_tmp(q) = [];
            qa_tmp(q) = [];
            
        end
        
        dxa = diff(spheredist(lata,lona));
        q = find(dxa == 0);
        
        if ~isempty(q)
            
            foo = qa([q q+1]);
            [~,I] = min(foo,[],2);
            dum = [q q+1];
            qout = NaN*ones(length(q),1);
            for k = 1:length(q)
                qout(k) = dum(k,I(k));
            end
            
            ta(qout) = [];
            lona(qout) = [];
            lata(qout) = [];
            qa(qout) = [];
            
        end
        