function n = getBehaviorSubset(e,b,start,stop)

n.beh = e.beh(start:stop);
n.stim = e.stim(start:stop);

n.n_tot = length(n.beh);
n.n_H = 0;
n.n_M = 0;
n.n_FA = 0;
n.n_CR = 0;
n.n_EL = 0;
for i=1:n.n_tot
    if cell2mat(n.beh(i)) == 'H'
        n.n_H = n.n_H+1;
    elseif cell2mat(n.beh(i)) == 'M'
        n.n_M = n.n_M+1;
    elseif cell2mat(n.beh(i)) == 'FA'
        n.n_FA = n.n_FA+1;    
    elseif cell2mat(n.beh(i)) == 'CR'
        n.n_CR = n.n_CR+1;
    elseif cell2mat(n.beh(i)) == 'EL'
        n.n_EL = n.n_EL+1;
    end
end

if n.n_H+n.n_M==0
    n.n_go = 1;
else
    n.n_go = n.n_H+n.n_M;
end
if n.n_FA+n.n_CR==0
    n.n_nogo = 1;
else
    n.n_nogo = n.n_FA+n.n_CR;
end
if n.n_H+n.n_FA==0
    n.n_lick = 1;
else
    n.n_lick = n.n_H+n.n_FA;
end
if n.n_M+n.n_CR==0
    n.n_nolick = 1;
else
    n.n_nolick = n.n_M+n.n_CR;
end
n.n_correct = n.n_H+n.n_CR;

n.r_H = n.n_H/n.n_go;  
n.r_M = n.n_M/n.n_go;  
n.r_FA = n.n_FA/n.n_nogo;  
n.r_CR = n.n_CR/n.n_nogo;
n.r_EL = n.n_EL/n.n_tot;
n.r_correct = n.n_correct/n.n_tot;
[n.dp,n.c] = dprime_simple(n.r_H,n.r_FA);
