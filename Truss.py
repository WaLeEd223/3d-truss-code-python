import numpy
import matplotlib.pyplot as plt 
import math 
from scipy import linalg

numpy.set_printoptions(10, suppress=True)

tn = int(input('Enter the total number of nodes : ')) #total nodes
te = int(input('Enter the total number of Elements : ')) #total elements
node_coordinate_matrix=[] #input matrix for getting coordinates of nodes
xco = [] #x co ordinate of nodes
yco = [] #y co ordinate of nodes
zco = [] #z co ordiaate of nodes
print("Enter nodal co-ordinate matrix")
for i in range(tn):
    row=list(map(int,input().split()))
    node_coordinate_matrix.append(row)  
node_coordinate_matrix.sort() 
for i in range(tn):
    x=node_coordinate_matrix[i][1]
    y=node_coordinate_matrix[i][2]
    z=node_coordinate_matrix[i][3]
    xco.append(x)
    yco.append(y)
    zco.append(z)




# material matrix of truss
material_matrix=[]

while 1:
    m_m_row=int(input("Enter the number of rows of material matrix and then enter the matrix"))
    if m_m_row>te:
        print("No. of rows of material matrix cannot be greater than total number of elements")
    else:
        break

for i in range(m_m_row):
    row=list(map(float,input().split()))
    material_matrix.append(row)      
material_matrix.sort()

# print(xco)
# print(yco)
# print(zco)
# print(material_matrix)
#element node matrix
print("Enter element node matrix with material properties")
element_nodal_point_matrix=[]
for i in range(te):
    row=list(map(float,input().split()))
    element_nodal_point_matrix.append(row)  
element_nodal_point_matrix.sort()

snofel = [] #start node of elements
enofel = [] #end node of elements
lenofel = [] #length of the element
elcon = [] #constant of the element
cosxofel = [] #cos x of element
cosyofel = [] #cos y of element
coszofel = [] #cos z of element

for i in range(te):  
    a = int(element_nodal_point_matrix[i][1])
    b = int(element_nodal_point_matrix[i][2])
    x1 = float(xco[a-1])
    y1 = float(yco[a-1])
    z1 = float(zco[a-1])
   
    x2 = float(xco[b-1])
    y2 = float(yco[b-1])
    z2 = float(zco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    indice = int(element_nodal_point_matrix[i][3]-1)
    A = material_matrix[indice][2]
    E = material_matrix[indice][1]
    con = A*E/l
    cosx = (x2-x1)/l
    cosy = (y2-y1)/l
    cosz = (z2-z1)/l
    snofel.append(a)
    enofel.append(b)
    lenofel.append(l)
    elcon.append(con)
    cosxofel.append(cosx)
    cosyofel.append(cosy)
    coszofel.append(cosz)
    
# ##print(snofel)
# ##print(enofel)
# ##print(lenofel)
# ##print(elcon)
# ##print(cosofel)
# ##print(sinofel)

elstmat = [] #element stiffness matrix

for i in range(te):
    
    cxs=float(cosxofel[i])**2
    cys=float(cosyofel[i])**2
    czs=float(coszofel[i])**2
    cxcy=float(cosxofel[i])*float(cosyofel[i])
    cxcz=float(cosxofel[i])*float(coszofel[i])
    cycz=float(cosyofel[i])*float(coszofel[i])
    
    mat = elcon[i]*numpy.array([[cxs, cxcy, cxcz, -cxs, -cxcy, -cxcz],
                                [cxcy, cys, cycz, -cxcy, -cys, -cycz],
                                [cxcz, cycz, czs, -cxcz, -cycz, -czs],
                                [-cxs, -cxcy, -cxcz, cxs, cxcy, cxcz],
                                [-cxcy, -cys, -cycz, cxcy, cys, cycz],
                                [-cxcz, -cycz, -czs, cxcz, cycz, czs]])

    elstmat.append(mat)
    # print(elstmat)
    # print("\n")


gstmatmap = []                          ## Global stiffness matrix mapping, gstmatmap will be the sqare matrix of tn*
for i in range(te):                     ## do this for each elements
    m = snofel[i]*3                     ## taking the start node of element(i) and multiply by 2
    n = enofel[i]*3                     ## taking the end node of element(i) and multiply by 2
    add = [m-2,m-1,m, n-2,n-1,n]              ## Address of columns and rows of gstmatmap for elemet(i)
                                            # if startnode is 1 and end node is 2 then add=[1,2,3,4,5,6]
                                            # if startnode is 1 and end node is 3 then add=[1,2,3,6,7,8]
    gmat = numpy.zeros((tn*3, tn*3))    ## global stiffness matrix loaded with zeros for element(i)
    elmat = elstmat[i]                  ## taking the element stiffness matrix of element(i)
    for j in range(6):                  
        for k in range(6):              
            a = add[j]-1                ## addressing row of GST matrix for element(i)
            b = add[k]-1                ## addressing column of GST matrix for element(i)
            gmat[a,b] = elmat[j,k]      ## updating the values in GST matrix with EST matrix of element(i)
    gstmatmap.append(gmat)              ## storing the resultant matrix in gstmatmap list
##    print(numpy.around(gmat, 3))

gsm = numpy.zeros((tn*3, tn*3))         ## creating an empyty GSM matrix
for mat in gstmatmap:
    gsm = gsm+mat                       ## adding all the matrix in the gstmatmap list
                                            # this will result in assembled stiffness matrix of the truss structure
GSM=numpy.matrix(gsm, dtype=numpy.float64)
# print('\nGlobal Stiffness Matrix of the Truss\n')
print(numpy.around(GSM, 10))

# #-----------------------Boundry condition and Loading---------------------#

displist = []
forcelist = []
for i in range(tn):
    a = str('u')+str(i+1)
    displist.append(a)
    b = str('v')+str(i+1)
    displist.append(b)
    c = str('w')+str(i+1)
    displist.append(c)
    
    d = str('fx')+str(i+1)
    forcelist.append(d)
    e = str('fy')+str(i+1)
    forcelist.append(e)
    f = str('fz')+str(i+1)
    forcelist.append(f)
    

# print(displist)
# print(forcelist)
    
print('\n\n________________Support Specifications______________\n')

dispmat = numpy.ones((tn*3,1))
tsupn = int(input('Enter the total number of nodes having supports : ')) #total number of supported nodes
supcondition = ['P = pinned',
                'X = x restrained ',
                'Y = y restrained ',
                'Z = z restrained ']
                
for i in range(tsupn):
    supn = int(input('\nEnter the node number of suuport : ')) #supported node
    for a in supcondition:
        print(a)
    condition = str(input('\nEnter the condition of the support : '))
    if condition in['P', 'p']:
        dispmat[supn*3-3, 0] = 0
        dispmat[supn*3-2, 0] = 0
        dispmat[supn*3-1, 0] = 0
    elif condition in['X', 'x']:
        dispmat[supn*3-3, 0] = 0
    elif condition in['Y', 'y']:
        dispmat[supn*3-2, 0] = 0
    elif condition in['Z','z']:
        dispmat[supn*3-1,0] = 0
    else:
        print('Please enter valid entries')
# print(dispmat)


print('\n_________________Loading____________________\n')
forcemat = numpy.zeros((tn*3,1))
tlon = int(input('Enter the total number of loaded nodes : ')) #total number of loaded nodes

for i in range(tlon):
    lon = int(input('\nEnter the node number of Loading : ')) #Loaded node
    fx = float(input('Enter the x load at this node in N : '))
    fy = float(input('Enter the y load at this node in N : '))
    fz = float(input('Enter the z load at this node in N : '))
    forcemat[lon*3-3, 0] = fx
    forcemat[lon*3-2, 0] = fy
    forcemat[lon*3-1, 0] = fz
# print(forcemat)    


# _________________Matrix Reduction_________________###


rcdlist = []
for i in range(tn*3):
    if dispmat[i,0] == 0:
        rcdlist.append(i)
rrgsm = numpy.delete(GSM, rcdlist, 0) #row reduction
crgsm = numpy.delete(rrgsm, rcdlist, 1) #column reduction
rgsm = crgsm #reduced global stiffness matrix
rforcemat = numpy.delete(forcemat, rcdlist, 0) #reduced force mat
rdispmat = numpy.delete(dispmat, rcdlist, 0) #reduced disp mat
# print(rforcemat,rdispmat)
# ###_______________Solving____________________###
# print(linalg.inv(rgsm))
dispresult = numpy.matmul(linalg.inv(rgsm), rforcemat)
# print(dispresult)
rin = 0
for i in range(tn*3):
    if dispmat[i,0] == 1:
        dispmat[i,0] = dispresult[rin,0]
        rin = rin+1

forceresult = numpy.matmul(GSM, dispmat)

# print(forceresult)


##____________________new co ordinates of nodes____________####
print(dispmat)
newxco = []
newyco = []
newzco = []
count = 0
for i in range(tn):
    k = xco[i]+dispmat[count,0]
    newxco.append(k)
    count = count+1
    l = yco[i]+dispmat[count,0]
    newyco.append(l)
    count = count+1
    m = zco[i]+dispmat[count,0]
    newzco.append(m)
    count = count+1


# ###____________________new length of memebers______________####
    
newlenofel = []
for i in range(te):
    a, b = snofel[i], enofel[i]
    x1 = float(newxco[a-1])
    y1 = float(newyco[a-1])
    z1 = float(newzco[a-1])
    x2 = float(newxco[b-1])
    y2 = float(newyco[b-1])
    z2 = float(newzco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    newlenofel.append(l)

print(newlenofel)
print(lenofel)

# ###______________strain in elements_______________________###
    
numpy.set_printoptions(3, suppress=False)

elstrain = numpy.zeros((te,1))
for i in range(te):
    elstrain[i,0] = (newlenofel[i]-lenofel[i])/(lenofel[i])
print('\n***Positive is Tensile\nNegetive is Compressive***\n')

print('\n\nStrain in the elements')
print(elstrain)
numpy.set_printoptions(3, suppress=True)

# ###__________________stress in elements______________________###

elstress = numpy.zeros((te,1))
for i in range(te):
    elstress[i,0] = E * elstrain[i,0]
    
print('\n\nStress in the elements')
print(elstress)

# ###_________________Member forces____________________#########

eforce = numpy.zeros((te,1))
for i in range(te):
    eforce[i,0] = A * elstress[i,0]

print('\n\nForce in the element')
print(eforce)
# xco.append(xco[0])
# yco.append(yco[0])
# newxco.append(newxco[0])
# newyco.append(newyco[0])
# plt.plot(xco,yco,label="Original Truss")
# plt.plot(newxco,newyco,label="Deformed Truss")
# # naming the x axis
# plt.xlabel('x - axis') 
# # naming the y axis 
# plt.ylabel('y - axis') 
# # giving a title to my graph 
# plt.title('Truss plot')

# # show a legend on the plot 
# plt.legend() 
  
# # function to show the plot 
# plt.show() 
