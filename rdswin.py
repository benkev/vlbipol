import numpy as np
import matplotlib.pyplot as plt

plt.ion()


finp = open("DIFX_59947_066600.s0000.b0000", "br")

while True:
    #
    # Read the SWIN file header
    #
    sync = np.fromfile(finp, dtype=np.uint32, count=1)
    if len(sync) == 0:
        break
    sync = sync[0]
    hver = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    blnum = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    mjd = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    secs = np.fromfile(finp, dtype=np.float64, count=1)[0]
    confix = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    srcix = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    freix = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    pol = finp.read(2) # Read polarization, like 'XY'
    pulsar_bin = np.fromfile(finp, dtype=np.uint32, count=1)[0]
    dweight = np.fromfile(finp, dtype=np.float64, count=1)[0]
    uvw = np.fromfile(finp, dtype=np.float64, count=3)
    vis = np.fromfile(finp, dtype=np.complex64, count=128)
    

    print('''
    20230103_v23003_v001_swin/v23003_01.difx/DIFX_59947_066600.s0000.b0000 
    Bytes   Type    Contains               Value
    1-4     Int     SYNC WORD              0X%8x 
    5-8     Int     BINARY HEADER VERSION  %u 
    9-12    Int     BASELINE NUM           %u 
    13-16   Int     MJD                    %u 
    17-24   Double  SECONDS                %f 
    25-28   Int     CONFIG INDEX           %u
    29-32   Int     SOURCE INDEX           %u
    33-36   Int     FREQ INDEX             %u
    37-38   Char[2] POLARISATION PAIR      %s
    39-42   Int     PULSAR BIN             %u
    43-50   Double  DATA WEIGHT            %f
    51-58   Double  U (METRES)             %f
    59-66   Double  V (METRES)             %f
    67-74   Double  W (METRES)             %f
    ''' % (sync, hver, blnum, mjd, secs, confix, srcix, freix, pol, pulsar_bin,
           dweight, *tuple(uvw)))

    # plt.figure()
    # plt.plot(vis.real[:], vis.imag[:], 'r.'); plt.grid(1)         
    # plt.show()
    # break;

finp.close()
