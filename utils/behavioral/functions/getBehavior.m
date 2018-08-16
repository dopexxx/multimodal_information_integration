function [b] = getBehavior(e,binSize)

if nargin == 1
  binSize = 30;
end

% behavior count

b.n_tot = e.n_trials;
b.n_H = 0;
b.n_M = 0;
b.n_FA = 0;
b.n_CR = 0;
b.n_EL = 0;
for i=1:b.n_tot
    if cell2mat(e.beh(i)) == 'H'
        b.n_H = b.n_H+1;
    elseif cell2mat(e.beh(i)) == 'M'
        b.n_M = b.n_M+1;
    elseif cell2mat(e.beh(i)) == 'FA'
        b.n_FA = b.n_FA+1;    
    elseif cell2mat(e.beh(i)) == 'CR'
        b.n_CR = b.n_CR+1;
    elseif cell2mat(e.beh(i)) == 'EL'
        b.n_EL = b.n_EL+1;
    end
end

% global rates

if b.n_H+b.n_M==0
    b.n_go = 1;
else
    b.n_go = b.n_H+b.n_M;
end
if b.n_FA+b.n_CR==0
    b.n_nogo = 1;
else
    b.n_nogo = b.n_FA+b.n_CR;
end
if b.n_H+b.n_FA==0
    b.n_lick = 1;
else
    b.n_lick = b.n_H+b.n_FA;
end
if b.n_M+b.n_CR==0
    b.n_nolick = 1;
else
    b.n_nolick = b.n_M+b.n_CR;
end
b.n_correct = b.n_H+b.n_CR;

b.r_H = b.n_H/b.n_go;  
b.r_M = b.n_M/b.n_go;  
b.r_FA = b.n_FA/b.n_nogo;  
b.r_CR = b.n_CR/b.n_nogo;
b.r_EL = b.n_EL/b.n_tot;
b.r_correct = b.n_correct/b.n_tot;
[b.dp,b.c] = dprime_simple(b.r_H,b.r_FA);

% local rates

n = 1;
b.lr_H = zeros(length(e.beh)-(binSize-1),1);
b.lr_M = zeros(length(e.beh)-(binSize-1),1);
b.lr_FA = zeros(length(e.beh)-(binSize-1),1);
b.lr_CR = zeros(length(e.beh)-(binSize-1),1);
b.lr_EL = zeros(length(e.beh)-(binSize-1),1);
b.l_dp = zeros(length(e.beh)-(binSize-1),1);
b.l_c = zeros(length(e.beh)-(binSize-1),1);
ln_H = 0;
ln_M = 0;
ln_FA = 0;
ln_CR = 0;
ln_EL = 0;

for i=1:length(e.beh)
    if i>=binSize
        for j=0:binSize-1
            if cell2mat(e.beh(i-j)) == 'H'
                ln_H = ln_H+1;
            elseif cell2mat(e.beh(i-j)) == 'M'
                ln_M = ln_M+1;    
            elseif cell2mat(e.beh(i-j)) == 'FA'
                ln_FA = ln_FA+1;    
            elseif cell2mat(e.beh(i-j)) == 'CR'
                ln_CR = ln_CR+1;    
            elseif cell2mat(e.beh(i-j)) == 'EL'
                ln_EL = ln_EL+1;   
            end
        end
        
        if ln_H+ln_M==0
            ln_go = 1;
        else
            ln_go = ln_H+ln_M;
        end
        if ln_FA+ln_CR==0
            ln_nogo = 1;
        else
            ln_nogo = ln_FA+ln_CR;
        end
        
        b.lr_H(n) = ln_H/ln_go;
        b.lr_M(n) = ln_M/ln_go;
        b.lr_FA(n) = ln_FA/ln_nogo;
        b.lr_CR(n) = ln_CR/ln_nogo;
        b.lr_EL(n) = ln_EL/binSize;
        [b.l_dp(n),b.l_c(n)] = dprime_simple(b.lr_H(n),b.lr_FA(n));
        ln_H = 0;
        ln_M = 0;    
        ln_FA = 0;    
        ln_CR = 0;
        ln_EL = 0;
        n = n+1;
    end    
end

% discrete local rates

b.dlr_H = zeros(fix(length(b.lr_H)/binSize),1);
b.dlr_M = zeros(fix(length(b.lr_M)/binSize),1);
b.dlr_FA = zeros(fix(length(b.lr_FA)/binSize),1);
b.dlr_CR = zeros(fix(length(b.lr_CR)/binSize),1);
b.dlr_EL = zeros(fix(length(b.lr_EL)/binSize),1);
b.dl_dp = zeros(fix(length(b.l_dp)/binSize),1);
b.dl_c = zeros(fix(length(b.l_c)/binSize),1);
n = 0;
for i=1:fix(length(b.lr_H)/binSize)
    b.dlr_H(i) = b.lr_H(1+n*binSize);
    b.dlr_M(i) = b.lr_M(1+n*binSize);
    b.dlr_FA(i) = b.lr_FA(1+n*binSize);
    b.dlr_CR(i) = b.lr_CR(1+n*binSize);
    b.dlr_EL(i) = b.lr_EL(1+n*binSize);
    b.dl_dp(i) = b.l_dp(1+n*binSize);
    b.dl_c(i) = b.l_c(1+n*binSize);
    n = n+1;
end

end











