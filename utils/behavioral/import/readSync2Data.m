function [sync2] = readSync2Data(file)

sync_file = lvm_import(file);
col1 = sync_file.Segment1.data(:,1);
col2 = sync_file.Segment1.data(:,2);

k=1;
for l = 2:length(col2)    
    if(col2(l-1)<1 & col2(l)>1)
        sync2(k) = col1(l); % for each trial: first frame in which TTL is present
        k = k+1;
    end    
end

sync2 = sync2';
end

