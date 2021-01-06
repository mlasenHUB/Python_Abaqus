import os
import numpy as np
from odbAccess import *
homedir = os.getcwd()
t = 't_82'
odb = openOdb(path=homedir+'/'+ t +'.odb')
frames = odb.steps['TANG_DISPL'].frames
file_con = t +'_History.dat'
#file_con = t +'_History.dat'
file_contact = homedir+'/' + file_con
dispFile = open(file_contact,'w')
dispFile.write('Ind_frame,RF1_LEFT,RF1_right,U1_punch,CSLIP1_ave_contact, CNAREA_sum_contact, dt, U2_tip_center, U2_tip_right, U2_tip_left, U2_active_base \n')

#for f in odb.steps['TANG_DISPL'].frames:
for f in range(0,len(frames)):
    
    frame = odb.steps['TANG_DISPL'].frames[f]
#Fields output
    U=frame.fieldOutputs['U']
    #U2=frame.fieldOutputs['U.U2']
    RF=frame.fieldOutputs['RF']
    CSLIP1 = frame.fieldOutputs['CSLIP1']
    COPEN = frame.fieldOutputs['COPEN']
    CNAREA=frame.fieldOutputs['CNAREA']
#regions-sets
    nodeset1=odb.rootAssembly.nodeSets['PLATE_CONTACT']
    node_punch=odb.rootAssembly.nodeSets['PUNCH_RP']
    node_left=odb.rootAssembly.nodeSets['PLATE_LEFT']
    node_right=odb.rootAssembly.nodeSets['PLATE_RIGHT']
    tip_center = odb.rootAssembly.instances['TIP-1'].nodes[1]#any node in tip Does the job cause it is a rigid body
    tip_right = odb.rootAssembly.instances['TIP-2'].nodes[1]#any node in tip Does the job cause it is a rigid body
    tip_left = odb.rootAssembly.instances['TIP-3'].nodes[1]#any node in tip Does the job cause it is a rigid body
    active_bases=odb.rootAssembly.nodeSets['PIEZO_BASES_ACTIVE']
    #tip_center = odb.rootAssembly.nodes[5] for some reason this does not work
    
    
#time
    dt = frame.frameValue
    
#Get Fields from regions - sets
    #UV_contact = U.getSubset(region=nodeset1).values
    UV_punch = U.getSubset(region=node_punch).values
    RF_left = RF.getSubset(region=node_left).values
    RF_right = RF.getSubset(region=node_right).values
    CSLIP_contact = CSLIP1.getSubset(region=nodeset1).values
    COPEN_contact = COPEN.getSubset(region=nodeset1).values
    CNAREA_contact = CNAREA.getSubset(region=nodeset1).values
    
    U_tip_center = U.getSubset(region=tip_center).values
    U_tip_right = U.getSubset(region=tip_right).values
    U_tip_left = U.getSubset(region=tip_left).values
    U_active_bases = U.getSubset(region=active_bases).values
    
    
    #Write coordinates
    #u1 = np.zeros(((len(UV_contact)),1))
    cslip1 = np.zeros(((len(CSLIP_contact)),1))
    cnarea = np.zeros(((len(CNAREA_contact)),1))
    copen = np.zeros(((len(COPEN_contact)),1))
    cslip1_adj = np.empty(((len(CSLIP_contact)),1))
    cslip1_adj[:] = np.nan 
    cnarea_adj = np.empty(((len(CNAREA_contact)),1))
    cnarea_adj[:] = np.nan 
    for n in range(len(CSLIP_contact)):
        #u1[n,0] = UV_contact[n].data[0]
        cslip1[n,0] = CSLIP_contact[n].data
        cnarea[n,0] = CNAREA_contact[n].data
        copen[n,0] = COPEN_contact[n].data
    
    U1_punch = UV_punch[0].data[0]# 0==>U1; 1==>U2;2==>U3
    U2_tip_center = U_tip_center[0].data[1]
    U2_tip_right = U_tip_right[0].data[1]
    U2_tip_left = U_tip_left[0].data[1]
    U2_active_bases = U_active_bases[0].data[1]
    
    cslip1_adj = cslip1[copen<=0]#just consider the values that are in contact copen negative
    cnarea_adj = cnarea[copen<=0]
    aux_cslip1 = np.abs(cslip1_adj)
    aux_cnarea = np.abs(cnarea_adj)
    cslip1_ave =np.mean(aux_cslip1[~np.isnan(aux_cslip1)])
    cnarea_sum =np.sum(aux_cnarea[~np.isnan(aux_cnarea)])
    RF1_left = RF_left[0].data[0]
    RF1_right = RF_right[0].data[0]
    
    
         
    dispFile.write('%10d, %15.10f, %15.10f, %15.10f, %15.10f, %15.10f , %15.10f , %15.10f , %15.10f , %15.10f , %15.10f\n' % \
    (f+1, RF1_left, RF1_right, U1_punch, cslip1_ave, cnarea_sum, dt, U2_tip_center, U2_tip_right, U2_tip_left, U2_active_bases))
    print(int(f))
    
print (file_contact)
dispFile.close()
