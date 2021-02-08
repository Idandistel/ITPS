import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
from Honey_Laminate import Honey_Laminate

Nz=1 



def Honeycomb (a,t,b,fiber_width,phi,fiber_clearance,n_times_x,n_times_y):
    # exporting 2d honeycomb cell with total length of 2a(1+cos60) over 2asin(60) pixels , mirroring 2 times for ease of programming
    #first creating the upper left cell
    def Honeycell(t,a,Nz):
        theta=60
        Nx=np.ceil(a*(1+np.cos(np.deg2rad(theta)))+b*np.sin(np.deg2rad(theta)))
        Ny=np.ceil(a*(np.sin(np.deg2rad(theta))))+b/2
        ratio=np.tan(np.deg2rad(theta))
        miniphase=np.zeros([int(Nx),int(Ny),Nz])
        miniphase[:int(a/2),int(b/2-t):int(b/2+t),:]=1 # BFS
    

        #diagonal part
        #bottom diagonal
        for dy in range(int((Ny-b/2))):
            
            dx=int((dy)/ratio)
            outerskin=int(a/2+dx-t+b*np.sin(np.deg2rad(theta)))
            miniphase[outerskin:(outerskin+2*t),dy,:]=1

        miniphase[int(outerskin):,int(Ny-t-b/2):int((Ny-b/2+t)),:]=1 # TFS
        #top diagonal part       
        for dy in range (int(Ny-b/2)):
            dx=int((dy)/ratio)
            outerskin=int(a/2+dx-t)
            miniphase[outerskin:(outerskin+2*t),dy+int(b/2),:]=1

   
        return miniphase

    
    def mirror(A,n_axis):
        Aflip=np.zeros(np.size(A))
        output=A
        if n_axis == 0 :
            Aflip =np.flipud(A)
        
        else:
            Aflip =np.fliplr(A)

        output=np.concatenate((output,Aflip),axis=n_axis)
        return output

    def linear_pattern(A,n_axis,n_times):
        output=A
        for i in range(n_times-1):
          output=np.concatenate((A,output),axis=n_axis)    
    
        return output




    theta=60
    #creating base quarter cell
    c=Honeycell(t,a,Nz)
    #plotting the result
    cfinal=c
    fig = plt.figure(1)
    plt.imshow(cfinal[:,:,0], cmap='viridis', interpolation='nearest')

    #mirror over y axis to the right
    cmid=mirror(c,1)
    #mirror over x downwards
    cfinal=mirror(np.flipud(cmid),0)
 
    #plotting the result
    fig = plt.figure(4)
    plt.imshow(cfinal[:,:,0], cmap='viridis', interpolation='nearest')
  
    #linear pattern the base cell n_times_x and n_times_y
    
    cfinal=linear_pattern(cfinal,0,n_times_x)
    cfinal=linear_pattern(cfinal,1,n_times_y)


    #adding fibers
    #creating coordinate vector
    coordinates=np.zeros([6,2])

    # #creating first cell
    # cfinal=Honey_Laminate(cfinal,coordinates,fiber_width,fiber_angle,fiber_clearance)
    for i in range (2):
        for j in range (4):
            #print('cell index ('+str(i)+str(j)+')')
            coordinates[0,:]=[int(a/2-a*(1+np.cos(np.deg2rad(theta))) )                             , int(a*np.sin(np.deg2rad(theta))-a*np.sin(np.deg2rad(theta))) ]
            coordinates[1,:]=[int(a/2+a*np.cos(np.deg2rad(theta)) -a*(1+np.cos(np.deg2rad(theta))))  , int(2*a*np.sin(np.deg2rad(theta))-a*np.sin(np.deg2rad(theta)))]
            coordinates[2,:]=[int(3*a/2 +a*np.cos(np.deg2rad(theta)) -a*(1+np.cos(np.deg2rad(theta))))  , int(2*a*np.sin(np.deg2rad(theta))-a*np.sin(np.deg2rad(theta)))]
            coordinates[3,:]=[int(3*a/2 +2*a*np.cos(np.deg2rad(theta))-a*(1+np.cos(np.deg2rad(theta))) )  ,  int(a*np.sin(np.deg2rad(theta))-a*np.sin(np.deg2rad(theta)))]
            coordinates[4,:]=[int(3*a/2 +a*np.cos(np.deg2rad(theta))  -a*(1+np.cos(np.deg2rad(theta)))) ,     int(-a*np.sin(np.deg2rad(theta)) )                           ]
            coordinates[5,:]=[int(a/2 +a*np.cos(np.deg2rad(theta)) -a*(1+np.cos(np.deg2rad(theta))))  ,       int(-a*np.sin(np.deg2rad(theta)) )                           ]  
            coordinates[:,0] +=j*a*(1+np.cos(np.deg2rad(theta)))+j*b*np.sin(np.deg2rad(theta))
            coordinates[:,1] +=i*(a*(2*np.sin(np.deg2rad(theta)))+b)+(j%2)*(a*np.sin(np.deg2rad(theta))+b/2)
            cfinal=Honey_Laminate(cfinal,coordinates,fiber_width,phi[i,j],fiber_clearance)

    Nx=np.ceil(a*(1+np.cos(np.deg2rad(theta)))+b*np.sin(np.deg2rad(theta)))
    Ny=np.ceil(a*(np.sin(np.deg2rad(theta))))+b/2
    #print(2*Nx*n_times_x,2*Ny*n_times_y,1)
    #print(cfinal.shape)

 


    fig = plt.figure(3)
    plt.imshow(cfinal[:,:,0], cmap='viridis', interpolation='nearest')
    plt.show()


    plt.show()
    return cfinal[:,:,:]
