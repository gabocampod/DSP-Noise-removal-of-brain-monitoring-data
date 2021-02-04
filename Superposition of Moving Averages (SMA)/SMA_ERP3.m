function [ReconSignal,AvgArtefact,eeg_in] = SMA_ERP1(INPUT,CurrentFreq)
%% Load data, the filename INPUT is the EEG+tACS DATA, for pre filtered data only

Channels=1;
Data=INPUT';

%% enter stimulation and EEG detailsdetails

T=length(Data(:,1));
Fs=500;
TIME=T/Fs;
% F=5;
F=CurrentFreq;
t=1/Fs:1/Fs:TIME;


%% Moving average windowing details
if F==40
NoSegments=(T/((Fs/F)*2));
AdjSegments=floor(0.05*NoSegments);   
    
else
NoSegments=(T/((Fs/F)*1));
AdjSegments=floor(0.05*NoSegments);
end

if mod(AdjSegments,2)==1
    AdjSegments=AdjSegments-1;
end
SegmentCount=AdjSegments+1;
timestep=(T/NoSegments);

%%

for i=1:1:Channels
    t1=1;
    t2=timestep;
    
    for j=1:1:NoSegments
        Segments(:,j)=Data(t1:t2,i);
        t1=t1+timestep;
        t2=t2+timestep;
    end
    
    for k=1:1:NoSegments
       if k<=(AdjSegments/2)
            S1=1;
            S2=AdjSegments;
            CurrentSegments(:,1:SegmentCount)=Segments(:,1:SegmentCount);
       else if k>=(NoSegments-(AdjSegments/2))
            CurrentSegments(:,1:SegmentCount)=Segments(:,(NoSegments-SegmentCount+1):NoSegments);
           else
            CurrentSegments(:,1:SegmentCount)=Segments(:,(k-AdjSegments/2):(k+AdjSegments/2));
    
           end
       end
       
        AvgSegment(:,k)=mean(CurrentSegments,2); 
        ScaledArtefact(:,k)=AvgSegment(:,k);
    end
    
    AvgArtefact(:,i)=reshape(ScaledArtefact,30000,1);
    TEMPSIG(:,i)=Data(:,i)-AvgArtefact(:,i);
    Signal(:,i)=TEMPSIG(:,i); 

end

eeg_in = Data(:,1);


ReconSignal=Signal;
end
