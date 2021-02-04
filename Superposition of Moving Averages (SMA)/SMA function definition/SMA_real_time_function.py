import numpy as np
import scipy.io as sio

FULLDATA_contents = sio.loadmat('DATEtime.mat')  # FULLDATA_contents is an array of structures
                                                 #DATAtime contains four structures ALPHA1 ERP1  SHAMALPHA1 AND SHAMERP1

# ALOCATE CORRECT STRUCTURES WITHIN DATA
ALPHA_struct = FULLDATA_contents['ALPHA1']  # ALPHA_STRUCTS is an array of structures with 1 element
ERP_struct = FULLDATA_contents['ERP1']
SHAMALPHA_struct = FULLDATA_contents['SHAMALPHA1']
SHAMERP_struct = FULLDATA_contents['SHAMERP1']

                                            #Each of this is a 1x1 struct with two fields A1 and TIME

# CONVERT TO NPARRAY VARIABLES
ALPHA_VALandTIME = ALPHA_struct[0, 0]  # ALPHA_VALandTIME is an array of two elements
ERP_VALandTIME = ERP_struct[0, 0]
SHAMALPHA_VALandTIME = SHAMALPHA_struct[0, 0]
SHAMERP_VALandTIME = SHAMERP_struct[0, 0]

                                        #now we are saying alphaValandTIME is element 1 of alpha struct in other words we copy all of
                                        #one into the other but the new one is easier to call

# EXAMPLES OF CALLS TO DATA STRUCTURES
# ALPHA_VALandTIME['A1']
# ALPHA_VALandTIME['TIME']

##print(ALPHA_VALandTIME['A1'].shape)
##print(ALPHA_VALandTIME['A1'])
##print(type(ALPHA_VALandTIME['A1']))


### ### ### ### ###
# SMA IMPLEMENTATION

data = np.zeros(shape=(1,ALPHA_VALandTIME['A1'].size))   #i.e #data is horizontal [0 0 0 0 ]

data = ALPHA_VALandTIME['A1']

datareshape = np.reshape(data,(data.size,1)) #now data is vertical [0 0 0 ]´

T = len(datareshape)  #lenght of DATA
F = 5  # Current frequency (could be asked to user)
Fs = 500  # sampling frequency (could be asked to user)


if F == 40:
    Nosegments = (T // ((Fs / F) * 2))  # in case freq is 40   #// is int divider
else:
    Nosegments = (T // ((Fs // F) * 1))  # for 5h and 10hz

AdjSegments = 10
SegmentCount = AdjSegments+1
timestep = (T // Nosegments)
Start = int(0.4 * Fs)
StartSeg = int(((T - Start) // timestep)+1)


# Assuming 1 channel

#used variables
segments = np.zeros(shape=(timestep,Nosegments))   #i.e timestep x Nosegments array
CurrentSegments = np.zeros(shape=(timestep, SegmentCount))
ScaledArtefact = np.zeros(shape=(timestep, (Nosegments-StartSeg)+1))
AvgArtefact = np.zeros(shape=(Start,1))
Signal = np.zeros(shape=(1,Start))
#start of code


t1 = 0
t2 = timestep



for j in range(Nosegments):

    segments[:, j] = datareshape[t1:t2,0]
    t1 = t1+timestep
    t2 = t2 + timestep

##print(segments)
a=0
for k in range (StartSeg,Nosegments+1,1):
    S1 = StartSeg - AdjSegments
    S2 = StartSeg
    CurrentSegments[:,0:SegmentCount]=segments[:,S1-1:S2]
    ScaledArtefact[:, a] = CurrentSegments.mean(axis=1)
    a = a+1

AvgArtefact = np.reshape( ScaledArtefact,(Start,1), order = 'F')
data = np.reshape(data,(1,T))
print(data.shape)
AvgArtefact = np.reshape(AvgArtefact,(1,Start))
print(AvgArtefact.shape)

Signal = data[0,(T-Start):T] - AvgArtefact



print(Nosegments)
print(timestep)
print(Start)
print(StartSeg)
print(segments.shape)
print(CurrentSegments.shape) ##ok
print(ScaledArtefact.shape)
print(AvgArtefact.shape)
print(Signal.shape)

print(Signal)