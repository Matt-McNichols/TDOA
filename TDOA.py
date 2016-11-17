import numpy as np
'''
def tau(P,T,m,v):
    # returns TAUm = wave speed * time difference
    # between reciever m and o


def ConstCoef(P,TAU,m):
    # Am = (2*Xm/TAUm) - (2*X1/TAU1)
    # bm = (2*Ym/TAUm) - (2*Y1/TAU1)
    # cm = (2*Zm/TAUm) - (2*Z1/TAU1)
    # dm = TAUm - TAU1 - (Xm^2 + Ym^2 +Zm^2)/TAUm + (X1^2 + Y1^2 +Z1^2)/TAU1

def systemEQ():
    # [ A2  B2  C2 ] * [ X ]= [ -D2 ]
    # [ A3  B3  C3 ] * [ Y ]= [ -D3 ]
    # [ A4  B4  C4 ] * [ Z ]= [ -D4 ]
    # return (x,y,z) the estimated location of the audio source

def calcT(P,E):
    return T;

def calcE(P,T):
    return E;



'''
def main():
# parameters:
#   P - array of hydrophone position vectors
#   T - array of times for when each hydrophone RX'ed the signal
#   v - wave speed
#   E - location of the emmitter
#   TAU - time shift given velocity of wave
#         TAUm = vTm - vTo
    P = np.array([[0,0,0],
                 [1,0,0],
                 [0,1,0],
                 [-1,0,0],
                 [0,-1,0]]);
    E = np.array([4,1,0]);
    T = np.zeros(5);
    v = 1481; # (m/s) speed of sound in water at 20 deg C

    print("P: ",P)
    print("E: ",E)



if __name__ == "__main__":
    main()
