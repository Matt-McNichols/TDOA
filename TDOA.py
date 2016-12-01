import numpy as np

def calcR(P,E):
    # R = | P - E |
    ro = np.square(P[0,0]+E[0]) + np.square(P[0,1]+E[1]) + np.square(P[0,2]+E[2])
    r1 = np.square(P[1,0]+E[0]) + np.square(P[1,1]+E[1]) + np.square(P[1,2]+E[2])
    r2 = np.square(P[2,0]+E[0]) + np.square(P[2,1]+E[1]) + np.square(P[2,2]+E[2])
    r3 = np.square(P[3,0]+E[0]) + np.square(P[3,1]+E[1]) + np.square(P[3,2]+E[2])
    r4 = np.square(P[4,0]+E[0]) + np.square(P[4,1]+E[1]) + np.square(P[4,2]+E[2])
    R0 = np.sqrt(ro);  R1 = np.sqrt(r1);
    R2 = np.sqrt(r2);  R3 = np.sqrt(r3);
    R4 = np.sqrt(r4);
    return np.array([R0,R1,R2,R3,R4])

def tau(P,T,v):
    # returns TAUm = wave speed * time difference
    # between reciever m and o
    TAU = np.zeros(4);
    for m in range(1,5):
        TAU[m-1] = v*T[m] - v*T[0];
    return TAU

def constCoef(P,TAU):
    # Am = (2*Xm/TAUm) - (2*X1/TAU1)
    # bm = (2*Ym/TAUm) - (2*Y1/TAU1)
    # cm = (2*Zm/TAUm) - (2*Z1/TAU1)
    # dm = TAUm - TAU1 - (Xm^2 + Ym^2 +Zm^2)/TAUm + (X1^2 + Y1^2 +Z1^2)/TAU1
    coef = np.zeros((3,3));
    D = np.zeros((3,1));
    subA = 2*P[1,0]/TAU[0];
    subB = 2*P[1,1]/TAU[0];
    subC = 2*P[1,2]/TAU[0];
    subD = (np.square(P[1,0]) + np.square(P[1,1]) + np.square(P[1,2]))/TAU[0]
    for m in range(2,5):
        Am = (2*P[m,0]/TAU[m-1]) - subA
        Bm = (2*P[m,1]/TAU[m-1]) - subB
        Cm = (2*P[m,2]/TAU[m-1]) - subC
        tmp = np.square(P[m,0]) + np.square(P[m,1]) + np.square(P[m,2])
        Dm = TAU[m-1] - TAU[0] - (tmp/TAU[m-1]) + subD
        coef[m-2,0] = Am;
        coef[m-2,1] = Bm;
        coef[m-2,2] = Cm;
        D[m-2] = Dm
    return (coef,D)


def calcE(coef,D):
    # [ A2  B2  C2 ] * [ X ]= [ D2 ]
    # [ A3  B3  C3 ] * [ Y ]= [ D3 ]
    # [ A4  B4  C4 ] * [ Z ]= [ D4 ]
    # return (x,y,z) the estimated location of the audio source
    approxE = np.linalg.lstsq(coef,D)
    print approxE;
    return approxE[0]


# given emmitter and reciever position
# calculate T (time signal arives at each reciever)
def problem1(P,E,v):
    R = calcR(P,E);
    T = np.divide(R,v);
    print ("R",R)
    print ("T",T)
    return T

# given reciever position and signal arival times
# calculate the location of emmitter
def problem2(P,T,v):
    # calculate tau
    # solve for coeficients
    # solve system of equations
    TAU = tau(P,T,v);
    coef,D = constCoef(P,TAU);
    E = calcE(coef,D)
    print ("TAU",TAU)
    print ("coef",coef);
    print ("D",D)
    return E




def main():
    # parameters:
    #   P - array of hydrophone position vectors
    #   T - array of times for when each hydrophone RX'ed the signal
    #   v - wave speed
    #   E - location of the emmitter
    #   TAU - time shift given velocity of wave
    #         TAUm = vTm - vTo
    P = np.array([[0,0,1],
                 [1,0,0],
                 [0,1,0],
                 [-1,0,0],
                 [0,-1,0]]);
    E = np.array([-4,6,5]);
    v = 1481; # (m/s) speed of sound in water at 20 deg C

    print('''
PARAMETERS
    ''')
    print ("P",P)
    print ("v",v)
    print ("E",E)
    print('''


PROBLEM 1: given P and E calculate T
    ''')
    T = problem1(P,E,v);
    print "T",T

    print('''


PROBLEM 2: given P and T calculate E
    ''')
    E = problem2(P,T,v);
    print ("E",E);
if __name__ == "__main__":
    main()
