import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np



def Honey_Laminate(A,coordinates,fiber_width,fiber_angle,fiber_clearance):
#Coordinates is a 2d(6,2) matrix with points coordinates , starting from left point and going in clock directin
#Anchor point at the middle point  
    anchor=np.zeros([2])
    anchor[0]=int((coordinates[0,0]+coordinates[3,0])/2) # x coordinate
    anchor[1]=int((coordinates[0,1]+coordinates[3,1])/2) # y coordinate
    h=coordinates[1,1]-coordinates[5,1]             # cell height
    w=coordinates[3,0]-coordinates[0,0]             #cell width
    
    if (fiber_angle >= 45 and fiber_angle <= 135 ):              # building through y direction
        t=int(np.ceil(np.sin(np.deg2rad(fiber_angle))*fiber_width/2))
        ratio=np.tan(np.deg2rad(fiber_angle))
        x_clearance=int(fiber_clearance*np.sin(np.deg2rad(fiber_angle)))
        for i in range(int(np.ceil (w/(2*(x_clearance))))): #w/(2*(fiber_clearance+fiber_width
            # print('fiber num'+str(i))
            dl=[x_clearance*i, 0]     # x coordinate clearance
            left_anchor=anchor-dl
            right_anchor=anchor+dl
            # print(anchor)
            # print(left_anchor)
            # print(right_anchor)
            for dy in range(int(np.ceil(h/2+t))):
                dx=int((dy)/ratio)
                # checking if the next stair reaches the border, if it is ,searching for border and setting it as the final stair
                right_index= (np.where(A[
                int(right_anchor[0]-t+dx):int(right_anchor[0]+t+dx),
                int(right_anchor[1]+dy)].flatten() == 1))
                #print('right index size'+str(np.size(right_index)))
                
                right_index = np.asarray(right_index)
                #print('right index='+str(right_index))
                
        
                
            
                
                A[int(left_anchor[0]-t-dx):int(left_anchor[0]+t-dx),int(left_anchor[1]-dy)]=1 #bottom left half
                A[int(right_anchor[0]-t+dx):int(right_anchor[0]+t+dx),int(right_anchor[1]+dy)]=1 #top right half

                if (np.size(right_index) >0 and 
                (right_index[0,0] == 0  )) :

                    # print(' right break on dy='+str(dy))
                    break
            if i>0 : #seperates the calculation for 2 couples half fibers, adding
                for dy in range(1,int(np.ceil(h/2))):
                    dx=int((dy)/ratio)
                    # print(h/2-dy)
                    # print(A.shape)
                    # checking if the next stair reaches the border, if it is ,searching for border and setting it as the final stair
                    left_index= (np.where(A[int(left_anchor[0]-t+dx):int(left_anchor[0]+t+dx),int(left_anchor[1]+dy)].flatten() == 1))
                    #print(np.size(index))
          
                    left_index = np.asarray(left_index)
                    #array_length = len(left_index[0,:])
                    # print(left_index)
                    
            
                    A[int(left_anchor[0]-t+dx):int(left_anchor[0]+t+dx),int(left_anchor[1]+dy)]=1 #top left half
                
                   
                    A[int(right_anchor[0]-t-dx):int(right_anchor[0]+t-dx),int(right_anchor[1]-dy)]=1  #bottom right half
                    if (np.size(left_index) >0 and 
                    (left_index[0,0] == 0  )) :
                        # print('left break on dy='+str(dy))
                        break 
    else : 
         if(fiber_angle<90):
             alpha=90+fiber_angle
         else :
            alpha=fiber_angle-90
         
         t=int(np.ceil(np.sin(np.deg2rad(alpha))*fiber_width/2) )                   
         ratio=np.tan(np.deg2rad(alpha))
         y_clearance=int(fiber_clearance*np.sin(np.deg2rad(alpha)))
         #fibers building loop
         for j in range(int(np.ceil (h/(2*(y_clearance))))):
            #  print('fiber num'+str(j))
             dl=[0, y_clearance*j]     # x coordinate clearance
             bottom_anchor=anchor-dl
             top_anchor=anchor+dl
            #  print(anchor)
            #  print(bottom_anchor)
            #  print(top_anchor)    
             for dx in range(int(np.ceil(w/2+t))): 
                 dy=int((dx)/ratio)
                 # checking if the next stair reaches the border, if it is ,searching for border and setting it as the final stair
                 top_index= (np.where(A[
                    int(top_anchor[0]-dx),
                    int(top_anchor[1]-t+dy):int(top_anchor[1]+t+dy)].flatten() == 1))
                 # print(np.size(top_index))
        
                 top_index = np.asarray(top_index)
                 #array_length = len(top_index[0,:])
        
                 

                 A[int(top_anchor[0]-dx),int(top_anchor[1]-t+dy):int(top_anchor[1]+t+dy)]=1 #top left fiber
                 A[int(bottom_anchor[0]+dx),int(bottom_anchor[1]-t-dy):int(bottom_anchor[1]+t-dy)]=1 # bottom right fiber
                 if (np.size(top_index) >0 and (top_index[0,0] == 0   )) :
                     #   print('top broke at'+str(dx))
                     break
             if j>0 : #adding additional 2 half fibers
                 for dx in range(1,int(np.ceil(w/2+t))):
                     dy=int((dx)/ratio)
                     bottom_index= np.where(A[
                     int(bottom_anchor[0]-dx),
                     int(bottom_anchor[1]-t+dy):int(bottom_anchor[1]+t+dy)].flatten() == 1) # bottom left fiber tip index
                     
                    #  print('bottom index'+str(np.size(bottom_index)))
                     bottom_index = np.asarray(bottom_index)
                     array_length = len(bottom_index[0,:])
                     
                     A[int(top_anchor[0]+dx),int(top_anchor[1]-t-dy):int(top_anchor[1]+t-dy)]=1 #top right fiber
                     A[int(bottom_anchor[0]-dx),int(bottom_anchor[1]-t+dy):int(bottom_anchor[1]+t+dy)]=1 # bottom left fiber
                     if (np.size(bottom_index) >0 and (bottom_index[0,0] == 0   )) :
                         # print('bottom broke at'+str(dx))
                         break

            

    return  A