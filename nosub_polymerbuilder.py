import numpy as np
import os
import sys
sys.dont_write_bytecode = True


def generate_rectangular_hcp(x_length, y_length, z_layers, particle_diameter, 
                            x_offset=0, y_offset=0, z_offset=0):

    a = particle_diameter
    radius = a / 2
    sqrt3_2 = np.sqrt(3)/2
    
    a1 = np.array([a, 0, 0])
    a2 = np.array([a/2, a*sqrt3_2, 0])
    a3 = np.array([0, 0, a*np.sqrt(8/3)])
    basis = np.array([
        [0, 0, 0],
        [a/2, a*sqrt3_2/3, a*np.sqrt(8/3)/2]
    ])

    n_x = int(np.ceil(x_length / a)) + int(2/sqrt3_2) + 3
    n_y = int(np.ceil(y_length / (a*sqrt3_2))) + int(2/sqrt3_2) + 3
    
    positions = []
    for i in range(-2*n_x, 2*n_x):
        for j in range(-2*n_y, 2*n_y):
            for k in range(z_layers):  
                for atom in basis:
                    pos = i*a1 + j*a2 + k*a3 + atom + np.array([x_offset, y_offset, z_offset])
                    
                    in_x = (x_offset + radius <= pos[0] <= x_offset + x_length - radius)
                    in_y = (y_offset + radius <= pos[1] <= y_offset + y_length - radius)
                    
                    if in_x and in_y:
                        positions.append(pos)
    
    positions = np.round(np.array(positions), decimals=8)
    positions = np.unique(positions, axis=0)
    
    return positions

def randomwalk(bond):

    theta = np.random.rand() * 2 * np.pi
    dz = 2 * np.random.rand() - 1
    dx = np.sqrt(1 - dz**2) * np.cos(theta)
    dy = np.sqrt(1 - dz**2) * np.sin(theta)
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    scale = bond / r
    dx *= scale
    dy *= scale
    dz *= scale
    return np.array([dx, dy, dz])

# def consbcp(x1=1,x2=1,y1=1,y2=1,a=30,b=30,c=30,d=30, run_dir="bulk.lammps"):
def consbcp(x1=1,x2=1,y1=1,y2=1,a=30,b=30,c=30,d=30, run_dir="."):
    mass=[0, 1.0, 1.0, 1.0,1.0]
    #N_polymer

    const = 5000
    L1A = a*x1+1
    L1B = b*x2+1
    L2A = c*y1+1
    L2B = d*y2+1

    if L1A==L1B:
        n1 = int(0.5*const/(L1A+L1B))
        n2 = int(0.5*const/(L2A+L2B))
    else:
        n1 = int(0.5*const*(L2A-L2B)/(L1B*L2A-L1A*L2B))
        n2 = int(0.5*const*(L1B-L1A)/(L1B*L2A-L1A*L2B))

    N1 = np.array([0,0,0,0,0,0,0,0,0,0]) 
    N2 = np.array([0,0,0,0,0,0,0,0,0,0])

    N1[3*(x1-1)+x2]=1
    N2[3*(y1-1)+y2]=1
    N1 = N1*n1
    N2 = N2*n2

    N1_polymer = np.sum(N1) #PMLiTFSI
    N2_polymer = np.sum(N2)
    #L_polymer = LA + LB

    L1_polymer = L1A + L1B
    L2_polymer = L2A + L2B

    #N_arm = N

    N1_arm = N1
    N2_arm = N2

    bcpbackboneAty = 1
    bcpbackboneBty = 2

    bcpbackboneAbty = 1
    bcpbackboneBbty = 1
    c = 0
    bcpbackboneBaty = 1
    bcpbackboneAaty = 1

    bonda = 1.0
    bondb = 1.0
    bondc = 1.0
    density = 0.88

    #N_total =  N_polymer * L_polymer 
    N_total =  N1_polymer * L1_polymer + N2_polymer * L2_polymer

    overall = N_total+1
    volume  = N_total/density # half of each layer
    x = 20
    y=  20
    z=(volume/(x*y))


    resid = 0
    #tmpcount_bcp=L_polymer+1
    #cbcp = np.zeros((N_polymer+1,tmpcount_bcp, 3))
    cbcp1 = np.zeros((N1_polymer+1,L1_polymer+1, 3))
    cbcp2 = np.zeros((N2_polymer+1,L2_polymer+1, 3))

    #alength = np.zeros(N_polymer)
    index = 0
    i=0
    print ("Total number of beads:", N_total)
    print ("Number of BCP:        ", N1_polymer+N2_polymer)
    print ("density:              ", density)
    print ("metric:             %10.4f %10.4f %10.4f\n\n" % (x, y, z))
    
    #position =       generate_rectangular_hcp(21,   21, 2, 1, -0.5, -0.5, 2)
    #x_coords = position[:, 0]
    #filter_mask1 = (x_coords > 0) &  (x_coords <= 10)
    #filter_mask2 = (x_coords > 10) & (x_coords <= 20)
#    filter_mask3 = (x_coords > 20) & (x_coords <= 30)
#    filter_mask4 = (x_coords > 30) & (x_coords <= 40)

    #position_ps1 = position[filter_mask1]
    #position_pmma1 = position[filter_mask2]
#    position_ps2 = position[filter_mask3]
#    position_pmma2 = position[filter_mask4]

    #position_ps1 =   generate_rectangular_hcp(10,   41, 2, 1, -0.5, -0.5, 2)
    #position_pmma1 = generate_rectangular_hcp(10,   41, 2, 1, 9.5,  -0.5, 2)
    #position_ps2 =   generate_rectangular_hcp(10,   41, 2, 1, 19.5, -0.5, 2)
    #position_pmma2 = generate_rectangular_hcp(10.5, 41, 2, 1, 29.5, -0.5, 2)
    #mxz = np.max(position_ps1[:,2])
    mxz = 0
    #print(len(position_pmma1),len(position_ps1))

    for pmmaarm in range(1,4):             #constract backbone of BCP POLYMER SPECIES 1
        for psarm in range(1,4):
            c+=1
            N_diff_type = N1_arm[c]
            for d in range(1,N_diff_type+1):
                i+=1
                cbcp1[i][L1A+1, :] = [-0.5, 0, 0]
                cbcp1[i][1,:] = [0.5,0,0]
                countps=0
                countpmma=0
                for j in range(1,psarm+1):  #loop of arm
                    armlength = int((L1B-1)/psarm)
                    for k in range(1,armlength+1): #loop of armlength
                        countps +=1                 
                        dx, dy, dz = randomwalk(bondb)
                        tmpind = L1A + countps + 1
                        if k ==1:
                            cbcp1[i][tmpind,:] = cbcp1[i][L1A+1,:] + [dx, dy, dz]  
                            if cbcp1[i][tmpind,0] >=0:
                                cbcp1[i][tmpind,0] += -dx
                            if cbcp1[i][tmpind,0] <=-10:
                                cbcp1[i][tmpind,0] += -dx
                        if k > 1:
                            cbcp1[i][tmpind,:] = cbcp1[i][tmpind - 1,:] + [dx, dy, dz]
                            if cbcp1[i][tmpind,0] >=0:
                                cbcp1[i][tmpind,0] += -dx
                            if cbcp1[i][tmpind,0] <=-10:
                                cbcp1[i][tmpind,0] += -dx

                for j in range(1,pmmaarm+1):                #construct A block
                    armlength = int((L1A-1)/pmmaarm)
                    for k in range(1,armlength+1):           #loop of armlength
                        countpmma +=1
                        dx, dy, dz = randomwalk(bondb)
                        tmpind = 1 + countpmma
                        if k == 1:
                            cbcp1[i][tmpind, :] = cbcp1[i][1,:] + [dx, dy, dz]
                            if cbcp1[i][tmpind,0] <=0:
                                cbcp1[i][tmpind,0] += -dx
                            if cbcp1[i][tmpind,0] >=10:
                                cbcp1[i][tmpind,0] += -dx                                
                        if k > 1:
                            cbcp1[i][tmpind, :] = cbcp1[i][tmpind - 1, :] + [dx, dy, dz]
                            if cbcp1[i][tmpind,0] <=0:
                                cbcp1[i][tmpind,0] += -dx
                            if cbcp1[i][tmpind,0] >=10:
                                cbcp1[i][tmpind,0] += -dx    

                z_singlechain_min = np.min(cbcp1[i][:, 2])
                
                cbcp1[i][:,2] -= z_singlechain_min 
                #print(np.min(cbcp[i][:, 2]))
    #print((cbcp[:, :, 2]))
    print(np.min(cbcp1[:, :, 2]))
    i=0
    c=0
    for pmmaarm in range(1,4):             #constract backbone of BCP POLYMER SPECIES 2
            for psarm in range(1,4):
                c+=1
                N_diff_type = N2_arm[c]
                for d in range(1,N_diff_type+1):
                    i+=1
                    cbcp2[i][L2A+1, :] = [-0.5, 0, 0]
                    cbcp2[i][1,:] = [0.5,0,0]
                    countps=0
                    countpmma=0
                    for j in range(1,psarm+1):  #loop of arm
                        armlength = int((L2B-1)/psarm)
                        for k in range(1,armlength+1): #loop of armlength
                            countps +=1                 
                            dx, dy, dz = randomwalk(bondb)
                            tmpind = L2A + countps + 1
                            if k ==1:
                                cbcp2[i][tmpind,:] = cbcp2[i][L2A+1,:] + [dx, dy, dz] 
                                if  cbcp2[i][tmpind,0] >=0:
                                    cbcp2[i][tmpind,0] += -dx
                                if  cbcp2[i][tmpind,0] <=-10:
                                    cbcp2[i][tmpind,0] += -dx
                            if k > 1:
                                cbcp2[i][tmpind,:] = cbcp2[i][tmpind - 1,:] + [dx, dy, dz]
                                if  cbcp2[i][tmpind,0] >=0:
                                    cbcp2[i][tmpind,0] += -dx
                                if  cbcp2[i][tmpind,0] <=-10:
                                    cbcp2[i][tmpind,0] += -dx

                    for j in range(1,pmmaarm+1):                #construct A block
                        armlength = int((L2A-1)/pmmaarm)
                        for k in range(1,armlength+1):           #loop of armlength
                            countpmma +=1
                            dx, dy, dz = randomwalk(bondb)
                            tmpind = 1 + countpmma
                            if k == 1:
                                cbcp2[i][tmpind, :] = cbcp2[i][1,:] + [dx, dy, dz]
                                if  cbcp2[i][tmpind,0] <=0:
                                    cbcp2[i][tmpind,0] += -dx
                                if  cbcp2[i][tmpind,0] >=10:
                                    cbcp2[i][tmpind,0] += -dx

                            if k > 1:
                                cbcp2[i][tmpind, :] = cbcp2[i][tmpind - 1, :] + [dx, dy, dz]
                                if  cbcp2[i][tmpind,0] <=0:
                                    cbcp2[i][tmpind,0] += -dx
                                if  cbcp2[i][tmpind,0] >=10:
                                    cbcp2[i][tmpind,0] += -dx                    
                    z_singlechain_min = np.min(cbcp2[i][:, 2])
                    
                    cbcp2[i][:,2] -= z_singlechain_min 
                    #print(np.min(cbcp[i][:, 2]))
        #print((cbcp[:, :, 2]))
    print(np.min(cbcp2[:, :, 2]))

    # asign the position of polymer
    index = 0
    posind = 0

    for n_p in N1:       
        for j in range(1,int(n_p)+1):                    
            index += 1
            posind+= 1

            zplane = np.random.rand() * z           
            yplane = np.random.rand() * y              
            cbcp1[index][:,1]+=yplane
            cbcp1[index][:,2]+=zplane +mxz


            if posind ==1:
                xplane = 10
                cbcp1[index][:,0]+=xplane
            if posind ==2:
                xplane = 20
                cbcp1[index][:,0] = -cbcp1[index][:,0] + xplane
                posind =0



    index =0
    for n_p in N2:       
        for j in range(1,int(n_p)+1):         
            index += 1
            posind+= 1
            zplane = np.random.rand() * z
            yplane = np.random.rand() * y  
            cbcp2[index][:,1]+=yplane
            cbcp2[index][:,2]+=zplane +mxz

            if posind ==1:
                xplane = 10
                cbcp2[index][:,0]+=xplane
            if posind ==2:
                xplane = 20
                cbcp2[index][:,0] = -cbcp2[index][:,0] + xplane
                posind =0
                
    index = 0
    overall=int(overall)
    #atomid = np.zeros(overall+len(position_pmma1)+len(position_ps1))
    atomid = np.zeros(overall)
    # chainid = np.zeros(overall+len(position_pmma1)+len(position_ps1))
    chainid = np.zeros(overall)
    # charge = np.zeros(overall+len(position_pmma1)+len(position_ps1))
    charge = np.zeros(overall)
    # typeid = np.zeros(overall+len(position_pmma1)+len(position_ps1))
    typeid = np.zeros(overall)
    # position = np.zeros((overall+len(position_pmma1)+len(position_ps1), 3))
    position = np.zeros((overall, 3))
    # imagenumber = np.zeros((overall+len(position_pmma1)+len(position_ps1), 3),dtype=int)
    imagenumber = np.zeros((overall, 3),dtype=int)

    bi = 0

    est =  N1_polymer * (L1_polymer-1) + N2_polymer * (L2_polymer-1) +1
    A1 = np.zeros(est)
    A2 = np.zeros(est)
    B1 = np.zeros(est)
    Btype = np.zeros(est)
    check = np.zeros(est)
    check2 = np.zeros(est)

    i=0
    c=0
    index_chain =0
    for pmmaarm in range(1,4):             #generate the atom info and bond info of BCP
        for psarm in range(1,4):
            c+=1
            N_diff_type = N1_arm[c]
            for d in range(1,N_diff_type+1): ##################################
                i+=1
                cbcptmp=cbcp1[i,:]
                index_chain =0

                index +=1
                index_chain +=1

                atomid[index]=index
                chainid[index]=i
                charge[index]=0
                typeid[index]=1
                position[index,:]=cbcptmp[index_chain,:]
                branchindexA = index

                for j in range(1,pmmaarm+1):                #construct A block
                    armlength = int((L1A-1)/pmmaarm)
                    for k in range(1,armlength+1):           #loop of armlength
                        
                        index+=1                        
                        index_chain +=1
                        
                        atomid[index]=index
                        chainid[index]=i
                        charge[index]=0
                        typeid[index]=1
                        position[index,:]=cbcptmp[index_chain,:]
                        
                        bi+=1
                        A1[bi]=index
                        if k ==1:
                            A2[bi]=branchindexA
                        if k > 1:
                            A2[bi]=index-1 
                        Btype[bi]=1
                
                index +=1
                index_chain +=1

                atomid[index]=index
                chainid[index]=i
                charge[index]=0
                typeid[index]=2
                position[index,:]=cbcptmp[index_chain,:]
                branchindexB = index

                for j in range(1,psarm+1):  #loop of arm
                    armlength = int((L1B-1)/psarm)
                    for k in range(1,armlength+1): #loop of armlength
                        
                        index+=1 
                        index_chain +=1

                        atomid[index]=index
                        chainid[index]=i
                        charge[index]=0
                        typeid[index]=2
                        position[index,:]=cbcptmp[index_chain,:]

                        bi+=1
                        A1[bi]=index
                        if k ==1:
                            A2[bi]=branchindexB
                        if k > 1:
                            A2[bi]=index-1
                        Btype[bi]=1 
                
                bi+=1
                A1[bi]=branchindexA
                A2[bi]=branchindexB
                Btype[bi]=2

    N_chain_record = i
    c=0
    index_chain =0
    for pmmaarm in range(1,4):             #generate the atom info and bond info of BCP
        for psarm in range(1,4):
            c+=1
            N_diff_type = N2_arm[c]
            for d in range(1,N_diff_type+1): ##################################
                i+=1
                cbcptmp=cbcp2[i-N_chain_record,:]
                index_chain =0

                index +=1
                index_chain +=1

                atomid[index]=index
                chainid[index]=i
                charge[index]=0
                typeid[index]=1
                position[index,:]=cbcptmp[index_chain,:]
                branchindexA = index

                for j in range(1,pmmaarm+1):                #construct A block
                    armlength = int((L2A-1)/pmmaarm)
                    for k in range(1,armlength+1):           #loop of armlength
                        
                        index+=1                        
                        index_chain +=1
                        
                        atomid[index]=index
                        chainid[index]=i
                        charge[index]=0
                        typeid[index]=1
                        position[index,:]=cbcptmp[index_chain,:]
                        
                        bi+=1
                        A1[bi]=index
                        if k ==1:
                            A2[bi]=branchindexA
                        if k > 1:
                            A2[bi]=index-1 
                        Btype[bi]=1
                
                index +=1
                index_chain +=1

                atomid[index]=index
                chainid[index]=i
                charge[index]=0
                typeid[index]=2
                position[index,:]=cbcptmp[index_chain,:]
                branchindexB = index

                for j in range(1,psarm+1):  #loop of arm
                    armlength = int((L2B-1)/psarm)
                    for k in range(1,armlength+1): #loop of armlength
                        
                        index+=1 
                        index_chain +=1

                        atomid[index]=index
                        chainid[index]=i
                        charge[index]=0
                        typeid[index]=2
                        position[index,:]=cbcptmp[index_chain,:]

                        bi+=1
                        A1[bi]=index
                        if k ==1:
                            A2[bi]=branchindexB
                        if k > 1:
                            A2[bi]=index-1
                        Btype[bi]=1 
                
                bi+=1
                A1[bi]=branchindexA
                A2[bi]=branchindexB
                Btype[bi]=2

    N_chain_count = chainid[index]
    index_chain=0
    # for i in range(0,len(position_ps1)):
    #     index +=1
    #     index_chain +=1
    #     atomid[index]=index
    #     chainid[index]=N_chain_count+index_chain
    #     charge[index]=0
    #     typeid[index]=4
    #     position[index,:]=position_ps1[i,:]

    # for i in range(0,len(position_ps2)):
    #     index +=1
    #     index_chain +=1
    #     atomid[index]=index
    #     chainid[index]=N_chain_count+index_chain
    #     charge[index]=0
    #     typeid[index]=4
    #     position[index,:]=position_ps2[i,:]

    # for i in range(0,len(position_pmma1)):
    #     index +=1
    #     index_chain +=1
    #     atomid[index]=index
    #     chainid[index]=N_chain_count+index_chain
    #     charge[index]=0
    #     typeid[index]=3
    #     position[index,:]=position_pmma1[i,:]

    # for i in range(0,len(position_pmma2)):
    #     index +=1
    #     index_chain +=1
    #     atomid[index]=index
    #     chainid[index]=N_chain_count+index_chain
    #     charge[index]=0
    #     typeid[index]=3
    #     position[index,:]=position_pmma2[i,:]

    for i in range(1,bi+1):
        check[i] =     (position[int(A1[i]), 0] - position[int(A2[i]), 0])**2 + \
                        (position[int(A1[i]), 1] - position[int(A2[i]), 1])**2 + \
                        (position[int(A1[i]), 2] - position[int(A2[i]), 2])**2
        if check[bi]>1.1:
            print('warning: bond between atom',A1[i],' ',A2[i],'is wrong, bond length is: ',check[i])
    zm=np.min(position[:,2])
    

    #bulk_file = os.path.join(run_dir, 'bulk.lammps')
    #with open(bulk_file, 'w') as file:
    save_path = os.path.join(run_dir, 'bulk.lammps')
    with open(save_path, 'w') as file:
        file.write('\n')
        file.write('#PMMA1 length     {:>8}\n'.format(L1A))
        file.write('#PS1 length       {:>8}\n'.format(L1B))
        file.write('#PMMA2 length     {:>8}\n'.format(L2A))
        file.write('#PS2 length       {:>8}\n'.format(L2B))
        file.write('#Number of BCP   {:>8}\n'.format(N1_polymer+N2_polymer))
        file.write('\n')
        #file.write('{:>8}     atoms\n'.format(N_total+len(position_pmma1)+len(position_ps1)))
        file.write('{:>8}     atoms\n'.format(N_total))
        file.write('{:>8}     bonds\n'.format(bi))
        file.write('{:>8}     angles\n'.format(0))
        file.write('{:>8}     dihedrals\n'.format(0))
        file.write('\n')

        file.write('{:>8}     atom types\n'.format(int(max(typeid))))
        file.write('{:>8}     bond types\n'.format(int(max(Btype))))
        file.write('{:>8}     angle types\n'.format(0))
        file.write('{:>8}     dihedral types\n'.format(0))
        file.write('\n')

        file.write('{:>16.3f} {:>8.3f} xlo xhi\n'.format(0, float(x)))    
        file.write('{:>16.3f} {:>8.3f} ylo yhi\n'.format(0, float(y)))
        file.write('{:>16.3f} {:>8.3f} zlo zhi\n'.format(np.min(position[:,2]-2), 50    )  )

        file.write('\n')
        file.write('Atoms\n')
        file.write('\n')

        # for i in range(1,overall+len(position_pmma1)+len(position_ps1)):
        for i in range(1,overall):
            file.write('{:>8} {:>8} {:>8} {:>8.2f} {:>8.3f} {:>8.3f} {:>8.3f}  \n'.format(int(atomid[i]),
                    int(chainid[i]), int(typeid[i]), float(charge[i]), position[i, 0], position[i, 1], position[i, 2]))

        file.write('\n')
        file.write('Bonds\n')
        file.write('\n')

        for i in range(1,bi+1):
            file.write('{:>8} {:>8} {:>8} {:>8}\n'.format(int(i), int(Btype[i]), int(A1[i]), int(A2[i])))
        
        file.write('\n')
        file.write('Masses\n')
        file.write('\n')
        for i in range(1,int(max(typeid))+1):
            file.write('     %i  %.2f\n'% (i,float(mass[i])))

consbcp()
