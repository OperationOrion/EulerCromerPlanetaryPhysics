"""

"""

import numpy as np
import matplotlib.pyplot as plt

#UNITS
#given to us in class--- im paying attention to your "NOTE' for unit conversion
DAY_to_SEC  = 86400.0                   # (s)     number of seconds in a day
YEAR_to_DAY =   365.00000000102244      # (days)  number of days in a year
DAY_to_YEAR =     0.0027397260273972603 # (years) number of years in a day
SEC_to_YEAR =     3.1709791983764586e-8 # (year)  number of years in a second
YEAR_to_SEC =   1.0/SEC_to_YEAR         # (s)
AU_to_km    =   149597870.700           # (km)
km_to_AU    =   1.0/AU_to_km            # (AU)
AU_to_m     =   AU_to_km * 1000         # (m)
m_to_AU     =   1.0/AU_to_m             # (AU)

#im going to solve this like a physics major
#and write my speeds to meters per second

#intial positions-------------------------------------------------------------
x0= 5.840e-1*AU_to_m
y0= 7.980e-1*AU_to_m
z0=-3.975e-5*AU_to_m
#in AU

#intial velocities------------------------------------------------------------
vx0=-1.417e-2*AU_to_m/DAY_to_SEC
vy0=1.009e-2*AU_to_m/DAY_to_SEC
vz0=2.368e-7*AU_to_m/DAY_to_SEC
#in AU/Day

#integral limits--------------------------------------------------------------
t01=0.0*YEAR_to_SEC
tn1=1.0*YEAR_to_SEC
#set two
t02=0.0*YEAR_to_SEC
tn2=3.0*YEAR_to_SEC
#increments given in final (in years)----------------------------------------
dt=[0.1,0.01,0.001,0.00001]

#universal constants
G= 6.674e-11 
m= 1.989e30 #kgs of sun

#Our function for Forward Euler!!!-------------------------------------------
def fc_solve_2body_fe(G,m,t0,tn,dt,x0,y0,z0,vx0,vy0,vz0):
    
    GM= G * m
    
    #How many times have I written what's below this semester?
    n=int((tn-t0)/dt)+1
    
    #intializing arrays for each position/velocity
    t=np.linspace(t0,tn,num=n)
    x=np.zeros(n,dtype=float)
    y=np.zeros(n,dtype=float)
    z=np.zeros(n,dtype=float)
    vx=np.zeros(n,dtype=float)
    vy=np.zeros(n,dtype=float)
    vz=np.zeros(n,dtype=float)
    
    #Given initial conditions for each postion/velocity
    x[0]=x0
    y[0]=y0
    z[0]=z0
    vx[0]=vx0
    vy[0]=vy0
    vz[0]=vz0
    
    #FE loop given in classroom code week 14 
    for k in range(n-1):
        #as you can see all we did was add z to make it 3D!!!
        r32k=np.sqrt(x[k]**2 + y[k]**2 + z[k])**(3.0)
        # for x pos/velocity
        vx[k+1]=vx[k]-dt*GM*x[k] / r32k
        x[k+1]=x[k]+dt*vx[k]
        # for y pos/velocity
        vy[k+1]=vy[k]-dt*GM*y[k] / r32k
        y[k+1]=y[k]+dt*vy[k]
        # for z pos/velocity
        vz[k+1]=vz[k]-dt*GM*z[k] / r32k
        z[k+1]=z[k]+dt*vz[k]
    
    return t,x,y,z,vx,vy,vz

#Our Function for Euler-Cromer!!!---------------------------------------------
def fc_solve_2body_ec(G,m,t0,tn,dt,x0,y0,z0,vx0,vy0,vz0):
    
    GM = G * m
    
    #number points from t0--->tn
    n = int((tn-t0)/dt)+1
    
    #intializing arrays for each position/velocity
    t  = np.linspace(t0,tn,num=n)
    x  = np.zeros(n,dtype=float)
    y  = np.zeros(n,dtype=float) 
    z  = np.zeros(n,dtype=float)    
    vx = np.zeros(n,dtype=float)
    vy = np.zeros(n,dtype=float)
    vz = np.zeros(n,dtype=float)

    #Given initial conditions for each position/velocity
    x[0]  = x0
    y[0]  = y0
    z[0] = z0
    vx[0] = vx0
    vy[0] = vy0
    vz[0] = vz0
    
    for k in range(n-1):
        #as you can see all we did was add z to make it 3D!!!
        r32k=np.sqrt(x[k]**2 + y[k]**2 + z[k])**(3.0)
        #for x
        vx[k+1] = vx[k]-dt*GM*x[k] / r32k
        x[k+1]  = x[k]+dt*vx[k+1]
        #for y
        vy[k+1] = vy[k]-dt*GM*y[k] / r32k
        y[k+1]  = y[k]+dt*vy[k+1]
        #for z
        vz[k+1] = vz[k]-dt*GM*z[k] / r32k
        z[k+1]  = z[k]+dt*vz[k+1]
    #NOTE THIS IS SLIGHTLY DIFFERENT FROM FE (look at [k+1] @end of velocities)
    return t,x,y,z,vx,vy,vz

#creating 2 plot functions for both methods, Forward Euler and Euler Cromer
#to save time iterating each graph of dt,integral, PvsT, PvsP

#-----------------------------------------------------------------------------
#GRAPHS FOR FORWARD EULER
#-----------------------------------------------------------------------------

def fc_plotting_fe(FE,dt,t,x,y,z,vx,vy,vz):

    #ONE/THREE YEAR GRAPHS----------------------------------------------------
    #x vs t
    plt.figure()
    plt.grid()
    plt.plot(t,x,"-k",label="X-cord vs time dt={},{}".format(dt,FE))
    plt.xlabel("t(s)")
    plt.ylabel("X(m)")
    plt.legend()
    plt.title('Forward Euler Method')
    #y vs t
    plt.figure()
    plt.grid()
    plt.plot(t,y,"-b",label="Y-cord vs time dt={},{}".format(dt,FE))
    plt.xlabel("t(s)")
    plt.ylabel("Y(m)")
    plt.legend()
    plt.title('Forward Euler Method')
    #z vs t
    plt.figure()
    plt.grid()
    plt.plot(t,z,"-g",label="Z-cord vs time dt={},{}".format(dt,FE))
    plt.xlabel("t(s)")
    plt.ylabel("Z(m)")
    plt.legend()
    plt.title('Forward Euler Method')
    #z vs x
    plt.figure()
    plt.grid()
    plt.plot(x,z,"-p",label="Z-cord vs X-cord dt={},{}".format(dt,FE))
    plt.xlabel("x(m)")
    plt.ylabel("z(m)")
    plt.legend()
    plt.title('Forward Euler Method')
    #z vs y
    plt.figure()
    plt.grid()
    plt.plot(y,z,"-y",label="Z-cord vs Y-cord dt={},{}".format(dt,FE))
    plt.xlabel("y(m)")
    plt.ylabel("z(m)")
    plt.legend()
    plt.title('Forward Euler Method')
    #y vs x
    plt.figure()
    plt.grid()
    plt.plot(x,y,"-m",label="Y-cord vs X-cord dt={},{}".format(dt,FE))
    plt.xlabel("x(m)")
    plt.ylabel("y(m)")
    plt.legend()
    plt.title('Forward Euler Method')
    
 #-----------------------------------------------------------------------------
 #GRAPHS FOR EULER CROMER
 #-----------------------------------------------------------------------------
 
def fc_plotting_ec(EC,dt,t,x,y,z,vx,vy,vz):
 #ONE/THREE YEAR GRAPHS-------------------------------------------------------
     #x vs t
    plt.figure()
    plt.grid()
    plt.plot(t,x,"-k",label="X-cord vs time dt={},{}".format(dt,EC))
    plt.xlabel("t(s)")
    plt.ylabel("X(m)")
    plt.legend()
    plt.title('Euler Cromer Method')
     #y vs t
    plt.figure()
    plt.grid()
    plt.plot(t,y,"-b",label="Y-cord vs time dt={},{}".format(dt,EC))
    plt.xlabel("t(s)")
    plt.ylabel("Y(m)")
    plt.legend()
    plt.title('Euler Cromer Method')
     #z vs t
    plt.figure()
    plt.grid()
    plt.plot(t,z,"-g",label="Z-cord vs time dt={},{}".format(dt,EC))
    plt.xlabel("t(s)")
    plt.ylabel("Z(m)")
    plt.legend()
    plt.title('Euler Cromer Method')
     #z vs y
    plt.figure()
    plt.grid()
    plt.plot(x,z,"-p",label="Z-cord vs X-cord dt={},{}".format(dt,EC))
    plt.xlabel("x(m)")
    plt.ylabel("z(m)")
    plt.legend()
    plt.title('Euler Cromer Method')
     #z vs y
    plt.figure()
    plt.grid()
    plt.plot(y,z,"-y",label="Z-cord vs Y-cord dt={},{}".format(dt,EC))
    plt.xlabel("y(m)")
    plt.ylabel("z(m)")
    plt.legend()
    plt.title('Euler Cromer Method')
     #y vs x
    plt.figure()
    plt.grid()
    plt.plot(x,y,"-m",label="Y-cord vs X-cord dt={},{}".format(dt,EC)) 
    plt.xlabel("x(m)")
    plt.ylabel("y(m)")
    plt.legend()
    plt.title('Euler Cromer Method')
    
#loop for every dt, from 0-1 or 0-3 for each formula type and constants
for n in range(len(dt)):
    
    dtn=dt[n]*YEAR_to_SEC
    
#FORWARD EULER DATA-----------------------------------------------------------
    #FE 1 YEAR---------------------------------------------------------------
    FE='F.E. 1 YEAR'
    t,x,y,z,vx,vy,vz=fc_solve_2body_fe(G,m,t01,tn1,dtn,x0,y0,z0,vx0,vy0,vz0) 
    fc_plotting_fe(FE,dt[n],t,x,y,z,vx,vy,vz)
    
    #3 YEAR------------(using tn2 now)---------------------------------------
    FE= 'F.E. 3 YEAR'
    t,x,y,z,vx,vy,vz=fc_solve_2body_fe(G,m,t01,tn2,dtn,x0,y0,z0,vx0,vy0,vz0) 
    fc_plotting_fe(FE,dt[n],t,x,y,z,vx,vy,vz)

#EULER CROMER DATA-----------------------------------------------------------
    #EC 1 YEAR---------------------------------------------------------------
    EC= 'E.C. 1 YEAR'
    t,x,y,z,vx,vy,vz=fc_solve_2body_ec(G,m,t01,tn1,dtn,x0,y0,z0,vx0,vy0,vz0)   
    fc_plotting_ec(EC,dt[n],t,x,y,z,vx,vy,vz)
    
    #EC 3 YEAR----------(using tn2 now)--------------------------------------
    EC= 'E.C. 3 YEAR'
    t,x,y,z,vx,vy,vz=fc_solve_2body_ec(G,m,t01,tn2,dtn,x0,y0,z0,vx0,vy0,vz0)   
    fc_plotting_ec(EC,dt[n],t,x,y,z,vx,vy,vz)

#That's all she wrote...   -M 5/11/2022

