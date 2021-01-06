import os
import numpy as np
#import cv2
from odbAccess import *
homedir = os.getcwd()

t = 't_82_adj'

odb = openOdb(path=homedir+'/'+ t +'.odb')
new_folder = homedir + '\\' + t + '_frames'
if not (os.path.exists(new_folder)):#if the file does not exist
    os.mkdir(new_folder)#Create the writing directory
##step_name = 'Transition'
step_name ='TANG_DISPL'
frames = odb.steps['TANG_DISPL'].frames#get the number of frames for the for loop
#for f in odb.steps['TANG_DISPL'].frames:
for f in range(0,len(frames)):
    frame = odb.steps[step_name].frames[f]
    CN_F=frame.fieldOutputs['CNORMF']
    Copen=frame.fieldOutputs['COPEN']
    Cpress=frame.fieldOutputs['CPRESS']
    U=frame.fieldOutputs['U']
    CSLIP1 = frame.fieldOutputs['CSLIP1']
    CNAREA=frame.fieldOutputs['CNAREA']
    
    #sets
    nodeset1=odb.rootAssembly.nodeSets['PLATE_CONTACT']
    
    #node_left=odb.rootAssembly.nodeSets['PLATE_LEFT']
    #node_right=odb.rootAssembly.nodeSets['PLATE_RIGHT']
    
    #Get data
    CnFV= CN_F.getSubset(region=nodeset1).values
    CopV= Copen.getSubset(region=nodeset1).values
    CPV= Cpress.getSubset(region=nodeset1).values 
    UV = U.getSubset(region=nodeset1).values
    CSLIP1 = CSLIP1.getSubset(region=nodeset1).values
    CnAV= CNAREA.getSubset(region=nodeset1).values
    #Write coordinates
    file_con = t +'_Fields_' +step_name +'_Frame_' + str(f) + '.dat'
    file_contact = new_folder +'/' + file_con
    file_coordinates = new_folder + '/' + t +'_Coordinates_' +step_name +'_Frame_'+ str(f) + '.dat'
    dispFile = open(file_contact,'w')
    dispFile.write('Ind,Node,N0_real,N0_gap,PRESS,U1, U2,U3, CSLIP1, CNAREA \n')     
    
    
    nodes_array = nodeset1.nodes[0]
    map = np.zeros((len(nodes_array),4))
    
    for iii in range(len(nodes_array)):
        map[iii][0] = nodes_array[iii].label
        map[iii][1] = nodes_array[iii].coordinates[0]
        map[iii][2] = nodes_array[iii].coordinates[1]
        map[iii][3] = nodes_array[iii].coordinates[2]
    
    
    for v in range(len(CnFV)):
        slip1 = CSLIP1[v].data
        cnarea = CnAV[v].data           
        normalforce = np.sqrt(CnFV[v].data[0]**2+CnFV[v].data[1]**2+CnFV[v].data[2]**2)
        press_read = CPV[v].data; u1 = UV[v].data[0];u2 = UV[v].data[1];u3 = UV[v].data[2];n0_print = normalforce        
        if normalforce< 1e-7: n0_print = -np.abs(CopV[v].data)
        dispFile.write('%10d, %10d, %15.10f, %15.10f, %15.10f, %15.10f, %15.10f, %15.10f, %15.10f, %15.10f \n' % \
        (v+1, CnFV[v].nodeLabel, normalforce, n0_print, press_read, u1, u2, u3, slip1, cnarea))
#     
#    print (file_coordinates)
#    print (file_con)
#    
    np.savetxt(file_coordinates,map)
    dispFile.close()
    
    print(int(f)) 
    print (file_con)
    
    dispFile.close()
