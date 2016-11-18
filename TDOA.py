import numpy as np

def calcR(hydrophonePositionArray,pingerLocation):
    # R = | hydrophonePositionArray - pingerLocation |

    #hax
    #values inside sq
    r = [0,0,0,0,0]
    R = [0,0,0,0,0]

    for i in range(5):

        r[i] = np.square(hydrophonePositionArray[i,0]+pingerLocation[0]) + \
               np.square(hydrophonePositionArray[i,1]+pingerLocation[1]) + \
               np.square(hydrophonePositionArray[i,2]+pingerLocation[2])
        
        R[i] = np.sqrt(r[i])

    return np.array([R[0],R[1],R[2],R[3],R[4]])

def tau(hydrophonePositionArray,T,soundVelocity):
    # returns TAUm = wave speed * time difference
    # between reciever m and o
    TAU = np.zeros(4);
    for m in range(1,5):
        TAU[m-1] = soundVelocity*T[m] - soundVelocity*T[0];
    return TAU

def constCoef(hydrophonePositionArray,TAU):
    # Am = (2*Xm/TAUm) - (2*X1/TAU1)
    # bm = (2*Ym/TAUm) - (2*Y1/TAU1)
    # cm = (2*Zm/TAUm) - (2*Z1/TAU1)
    # dm = TAUm - TAU1 - (Xm^2 + Ym^2 +Zm^2)/TAUm + (X1^2 + Y1^2 +Z1^2)/TAU1
    coef = np.zeros((3,3));
    D    = np.zeros((3,1));

    subA = 2*hydrophonePositionArray[1,0]/TAU[0];
    subB = 2*hydrophonePositionArray[1,1]/TAU[0];
    subC = 2*hydrophonePositionArray[1,2]/TAU[0];
    subD = (np.square(hydrophonePositionArray[1,0]) + np.square(hydrophonePositionArray[1,1]) + np.square(hydrophonePositionArray[1,2]))/TAU[0]
    
    for m in range(2,5):
        Am = (2*hydrophonePositionArray[m,0]/TAU[m-1]) - subA
        Bm = (2*hydrophonePositionArray[m,1]/TAU[m-1]) - subB
        Cm = (2*hydrophonePositionArray[m,2]/TAU[m-1]) - subC
        tmp = np.square(hydrophonePositionArray[m,0]) + np.square(hydrophonePositionArray[m,1]) + np.square(hydrophonePositionArray[m,2])
        Dm = TAU[m-1] - TAU[0] - (tmp/TAU[m-1]) + subD
        coef[m-2,0] = Am;
        coef[m-2,1] = Bm;
        coef[m-2,2] = Cm;
        D[m-2] = Dm

    return (coef,D)


def calcE(coef,D):
    # [ A2  B2  C2 ] * [ X ]= [ -D2 ]
    # [ A3  B3  C3 ] * [ Y ]= [ -D3 ]
    # [ A4  B4  C4 ] * [ Z ]= [ -D4 ]
    # return (x,y,z) the estimated location of the audio source
    approxE = np.linalg.lstsq(coef,D)
    return approxE[0]


# given emmitter and reciever position
# calculate T (time signal arives at each reciever)
def problem1(hydrophonePositionArray,pingerLocation,soundVelocity):
    R = calcR(hydrophonePositionArray,pingerLocation);
    T = np.divide(R,soundVelocity);
    print ("R",R)
    print ("T",T)
    return T

# given reciever position and signal arival times
# calculate the location of emmitter
def problem2(hydrophonePositionArray,T,soundVelocity):
    # calculate tau
    # solve for coeficients
    # solve system of equations
    TAU = tau(hydrophonePositionArray,T,soundVelocity);
    coef,D = constCoef(hydrophonePositionArray,TAU);
    pingerLocation = calcE(coef,D)
    print ("TAU",TAU)
    print ("coef",coef);
    print ("D",D)
    return pingerLocation




def main():
    # parameters:

    #   hydrophonePositionArray - array of hydrophone position vectors [x,y,z]?
    #   T                       - array of times for when each hydrophone RX'ed the signal
    #   soundVelocity           - wave speed
    #   pingerLocation          - location of the emmitter
    #   TAU                     - time shift given velocity of wave
    #   TAUm = vTm - vTo

    hydrophone0 = [0,0,0]       #center hydrophone
    hydrophone1 = [1,0,0]
    hydrophone2 = [0,1,0]
    hydrophone3 = [-1,0,0]
    hydrophone4 = [0,-1,0]

    hydrophoneArrangement = [hydrophone0,hydrophone1,hydrophone2,hydrophone3,hydrophone4]

    hydrophonePositionArray = np.array(hydrophoneArrangement)

    pingerLocation = np.array([-44,6,5]);

    #speed that sound travels underwater in meters/second @ 20 degrees Celsius
    soundVelocity = 1481;

    print "\n PARAMETERS \n"

    print ("hydrophonePositionArray",hydrophonePositionArray)
    print ("soundVelocity",soundVelocity)
    print ("pingerLocation",pingerLocation)
    print('''


PROBLEM 1: given hydrophonePositionArray and pingerLocation calculate T
    ''')
    T = problem1(hydrophonePositionArray,pingerLocation,soundVelocity);
    print "T",T

    print('''


PROBLEM 2: given hydrophonePositionArray and T calculate pingerLocation
    ''')
    pingerLocation = problem2(hydrophonePositionArray,T,soundVelocity);
    print ("pingerLocation",pingerLocation);
if __name__ == "__main__":
    main()
