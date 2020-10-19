import cv2
import math
import numpy as np

court = cv2.imread("tennis_court.png")
#cv2.imshow("Tennis court: ",court)
#cv2.waitKey(0)

def draw_trajectory(image,Sx,Sy):

    px_per_m = 46.31
    gnd_y = 396

    for i in range(Sx.shape[0]):

        X_m = Sx[i]
        X_px = int(X_m*px_per_m)

        Y_m = Sy[i]
        Y_px = int(Y_m*px_per_m)
        Y_px = gnd_y - Y_px


        cv2.circle(image,(X_px,Y_px),3,(0,255,0),-1)

    cv2.imshow("Trajectory:",image)
    cv2.imwrite('trajectory.png',image)
    cv2.waitKey(0)

points = 200


V = 45
rps = 200
angle_degrees =7 #degrees


theta =  (math.pi/180)*angle_degrees
Vx_in = V * math.cos(theta)
Vy_in = V * math.sin(theta)

T = 2
t = np.linspace(0, T, points)


Vx = np.zeros(t.shape[0])
Vy = np.zeros(t.shape[0])
Sx = np.zeros(t.shape[0])
Sy = np.zeros(t.shape[0])
dragx = np.zeros(t.shape[0])
dragy = np.zeros(t.shape[0])
ground = np.zeros(t.shape[0])

Vx[0] = Vx_in
Vy[0] = Vy_in
Sx[0] = 0
Sy[0] = 1

C = 0.65

density_air = 1.21
m = 0.059
g = 9.8
A = 0.0075
clash = 1
radius = 0.0342
g = 9.8



V_spin = radius*rps

for i in range (1,t.shape[0]):

    cross_area = math.pi * (radius**2)
    delta_t = t[i] - t[i-1]

    theta = math.atan(Vy[i-1]/Vx[i-1])


    V_net = (Vx[i-1]**2 + Vy[i-1]**2)**0.5

    '''Cd = 0.55 + 1/ (22.5 + 4.2 * (V_net / V_spin)**2.5)**0.4
    Fd = Cd *cross_area * density_air * (V_net**2) /2'''
    Fd = 0.0014 * V_net**2


    Cl = 1 / (2 + (V_net/V_spin))
    Fl = Cl * cross_area * density_air * (V_net**2)/2

    Fx = -Fd*math.cos(theta) - Fl*math.sin(theta)
    Fy =  Fl*math.cos(theta) - Fd*math.sin(theta) - m*g

    Vx[i] = Vx[i-1] + (1/m)* (Fx) *delta_t
    Vy[i] = Vy[i-1] + (1/m)* (Fy) *delta_t

    Sx[i] = Sx[i - 1] + Vx[i - 1] * delta_t + 0.5 * (Fx/m) * delta_t ** 2
    Sy[i] = Sy[i - 1] + Vy[i - 1] * delta_t + 0.5 * (Fy/m) * delta_t ** 2

    if (Sy[i] <= 0):
        #hitting ground....
        print("Bounce at: T= ",t[i]," secs")
        print("Distance: ",Sx[i])
        Vy[i] = -Vy[i]
        Vx[i] = Vx[i]

    '''if (abs(Sx[i] - 11.885) < 0.1) and (Sy[i] <= 0.914):


        Vx[i] = -Vx[i]*0.3
        Vy[i] = Vy[i]*0.1'''






##Sx = Sx[0:clash]
#Sy = Sy[0:clash]


draw_trajectory(court,Sx,Sy)

