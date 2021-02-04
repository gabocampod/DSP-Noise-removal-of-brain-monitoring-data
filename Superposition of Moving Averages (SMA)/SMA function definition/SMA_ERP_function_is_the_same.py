import numpy as np

#def SMA_REAL_TIME(INPUT,currentfreq,SamplingFreq):

INPUT = np.arange(20000,40000,1)
currentfreq = 5
SamplingFreq =  500

data = np.zeros(shape=(1, INPUT.size))  # i.e #data is horizontal [0 0 0 0 ]

data = INPUT

datareshape = np.reshape(data, (data.size, 1))  # now data is vertical [0 0 0 ]´


T = len(datareshape)  # length of DATA
F = currentfreq       # Current frequency (could be asked to user)
Fs = SamplingFreq     # sampling frequency (could be asked to user)

TIME = T/Fs
my_new_t = np.arange(1/Fs,TIME,1/Fs)

if F == 40:
    Nosegments = (T // ((Fs / F) * 2))  # in case freq is 40   #// is int divider
    AdjSegments = np.floor(0.05*Nosegments)     #ADJ is an integer still
else:
    Nosegments = (T // ((Fs // F) * 1))  # for 5h and 10hz
    AdjSegments = np.floor(0.05 * Nosegments)

if np.mod(AdjSegments,2)==1:            #if adj segments is odd then:
    AdjSegments=AdjSegments-1;          #make adjsegments even

SegmentCount = int(AdjSegments + 1)     #still an int
timestep = (T // Nosegments)
# Assuming 1 channel


# used variables
segments = np.zeros(shape=(timestep, Nosegments))  # i.e timestep x Nosegments array
CurrentSegments = np.zeros(shape=(timestep, SegmentCount)) #should still be true
ScaledArtefact = np.zeros(shape=(timestep, Nosegments))
AvgArtefact = np.zeros(shape=(T, 1))
Signal = np.zeros(shape=(1, T))
# start of code

t1 = 0
t2 = timestep

for j in range(Nosegments):             #Has not changed
    segments[:, j] = datareshape[t1:t2, 0]
    t1 = t1 + timestep
    t2 = t2 + timestep

for k in range(1,Nosegments+1,1):   #I added that + 1
    if k<=(AdjSegments/2):           #remember ADJ segments is even so this always gives int
        S1 = 0
        S2 = AdjSegments
        CurrentSegments[:, 0:SegmentCount] = segments[:, 0:SegmentCount]
    else:
        if k>=(Nosegments-(AdjSegments/2)):
              CurrentSegments[:,0:SegmentCount] = segments[:,(Nosegments-SegmentCount):Nosegments]
        else:
            CurrentSegments[:, 0:SegmentCount] = segments[:,(int(k-AdjSegments/2)-1):int((k+AdjSegments/2))]

        ScaledArtefact[:, k-1] = CurrentSegments.mean(axis=1)



AvgArtefact = np.reshape(ScaledArtefact, (T, 1), order='F')  #DEFINE AVGARTEFACT

AvgArtefact = np.reshape(AvgArtefact, (1, T))  #TURN INTO HORIZONTAL VECTOR

Signal = data - AvgArtefact

#print(Signal)


print(Nosegments)
print(timestep)
print(AdjSegments)
print(SegmentCount)
print(segments.shape)
print(CurrentSegments.shape)  ##ok
print(ScaledArtefact.shape)
print(AvgArtefact.shape)
print(Signal.shape)

print(Signal)
