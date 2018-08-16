function bs = getBehaviorByStimulus(e)

bs.n_N = 0;
bs.n_NH = 0;
bs.n_NM = 0;
bs.n_NFA = 0;
bs.n_NCR = 0;
bs.n_NEL = 0;
bs.n_V = 0;
bs.n_VH = 0;
bs.n_VM = 0;
bs.n_VFA = 0;
bs.n_VCR = 0;
bs.n_VEL = 0;
bs.n_S = 0;
bs.n_SH = 0;
bs.n_SM = 0;
bs.n_SFA = 0;
bs.n_SCR = 0;
bs.n_SEL = 0;
bs.n_VS = 0;
bs.n_VSH = 0;
bs.n_VSM = 0;
bs.n_VSFA = 0;
bs.n_VSCR = 0;
bs.n_VSEL = 0;

for trial = 1:e.n_trials
    if cell2mat(e.stim(trial)) == 'N'
        bs.n_N = bs.n_N+1;
        if cell2mat(e.beh(trial)) == 'H'
            bs.n_NH = bs.n_NH+1;
        elseif cell2mat(e.beh(trial)) == 'M'
            bs.n_NM = bs.n_NM+1;
        elseif cell2mat(e.beh(trial)) == 'FA'
            bs.n_NFA = bs.n_NFA+1;        
        elseif cell2mat(e.beh(trial)) == 'CR'
            bs.n_NCR = bs.n_NCR+1;    
        elseif cell2mat(e.beh(trial)) == 'EL'
            bs.n_NEL = bs.n_NEL+1;
        end
    elseif cell2mat(e.stim(trial)) == 'V'
        bs.n_V = bs.n_V+1;
        if cell2mat(e.beh(trial)) == 'H'
            bs.n_VH = bs.n_VH+1;
        elseif cell2mat(e.beh(trial)) == 'M'
            bs.n_VM = bs.n_VM+1;
        elseif cell2mat(e.beh(trial)) == 'FA'
            bs.n_VFA = bs.n_VFA+1;        
        elseif cell2mat(e.beh(trial)) == 'CR'
            bs.n_VCR = bs.n_VCR+1;    
        elseif cell2mat(e.beh(trial)) == 'EL'
            bs.n_VEL = bs.n_VEL+1;
        end
    elseif cell2mat(e.stim(trial)) == 'S'
        bs.n_S = bs.n_S+1;
        if cell2mat(e.beh(trial)) == 'H'
            bs.n_SH = bs.n_SH+1;
        elseif cell2mat(e.beh(trial)) == 'M'
            bs.n_SM = bs.n_SM+1;
        elseif cell2mat(e.beh(trial)) == 'FA'
            bs.n_SFA = bs.n_SFA+1;        
        elseif cell2mat(e.beh(trial)) == 'CR'
            bs.n_SCR = bs.n_SCR+1;    
        elseif cell2mat(e.beh(trial)) == 'EL'
            bs.n_SEL = bs.n_SEL+1;
        end     
    elseif cell2mat(e.stim(trial)) == 'V+S'
        bs.n_VS = bs.n_VS+1;
        if cell2mat(e.beh(trial)) == 'H'
            bs.n_VSH = bs.n_VSH+1;
        elseif cell2mat(e.beh(trial)) == 'M'
            bs.n_VSM = bs.n_VSM+1;
        elseif cell2mat(e.beh(trial)) == 'FA'
            bs.n_VSFA = bs.n_VSFA+1;        
        elseif cell2mat(e.beh(trial)) == 'CR'
            bs.n_VSCR = bs.n_VSCR+1;    
        elseif cell2mat(e.beh(trial)) == 'EL'
            bs.n_VSEL = bs.n_VSEL+1;
        end  
    end
end

% correct trials for each condition

bs.n_Nc = bs.n_NH+bs.n_NCR;
bs.n_Vc = bs.n_VH+bs.n_VCR;
bs.n_Sc = bs.n_SH+bs.n_SCR;
bs.n_VSc = bs.n_VSH+bs.n_VSCR;

bs.r_Nc = bs.n_Nc/bs.n_N;
bs.r_Vc = bs.n_Vc/bs.n_V;
bs.r_Sc = bs.n_Sc/bs.n_S;
bs.r_VSc = bs.n_VSc/bs.n_VS;

% early licks for each condition

bs.r_NEL = bs.n_NEL/bs.n_N;
bs.r_VEL = bs.n_VEL/bs.n_V;
bs.r_SEL = bs.n_SEL/bs.n_S;
bs.r_VSEL = bs.n_VSEL/bs.n_VS;

end