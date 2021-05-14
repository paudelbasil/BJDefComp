
defExport
! Created by Basil J. Paudel            5/13/2021
! This macro will export the nodal coordinates and the deformation results
!         to a CSV file
! Usage:
! First place the defMacro.F in the Ansys working directory
!         *ulib,defMacro,F,
!         *use, defExport

! DEFINE OUTPUT PARAMETERS
fname='Results_T'

! DO NOT CHANGE BELOW THIS LINE

/post1

/uis,msgpop,1
*MSG,NOTE,1 
Exporting deformation results to Results_T**.csv

! GET NUMBER OF STORED DATA SETS TO LOOP
*GET, numSets,active,0,set,nset

! GET NODAL INFORMATION
*GET, maxNode, NODE,0,NUM,MAX
*GET, minNode, NODE,0,NUM,MIN

! FIRST SETUP A TABLE TO STORE DATA
*DIM, XPT,TABLE, maxNode,3
*DIM, XPT_N,ARRAY, maxNode

nodeNum=1
ALLSEL

! NEED TO LOOP OVER ALL DATA SETS TO EXPORT
SET,FIRST

*DO,cntDS,numSets-1,numSets

*do,nodeNum,minNode,maxNode

! GET NODAL COORDINATES
xx=NX(nodeNum)
yy=NY(nodeNum)
zz=NZ(nodeNum)

! GET NODAL DISPLACEMENTS
sux=UX(nodeNum)
suy=UY(nodeNum)
suz=UZ(nodeNum)

*SET,XPT_N(nodeNum),nodeNum
*SET,XPT(nodeNum,1),xx+sux
*SET,XPT(nodeNum,2),yy+suy
*SET,XPT(nodeNum,3),zz+suz

*enddo

! GET INFO ON SOLUTION TIME
*GET,solTime, ACTIVE, 0, SET, TIME

fnum = chrval(cntDS)
fTime = chrval(solTime)
fullfname = strcat(fname,fTime)

! NOW WRITE TO FILE
*CFOPEN,fullfname,'csv'           ! Open file for writing
*VWRITE,'!', 'Nodal Information Dump'
%1C %67C
*VWRITE,'!', 'Total Nodes =', maxNode-minNode
%C %56C %10I
*VWRITE, '!', 'Node', 'X','Y','Z'
%C%8C%20C%20C%20C
*CFWRITE,nblock,3,,maxNode-minNode

! THIS IS A WEIRD WAY OF DOING THIS
*VWRITE,'1i9,3e20', '9e3'
%8C.%3C

*VWRITE, XPT_N(1),XPT(1,1),XPT(1,2),XPT(1,3)
%9I%20.9G%20.9G%20.9G

*CFWRITE,'-1'
*CFWRITE,'! END OF NODAL INFORMATION '

*CFCLOSE  ! Close the File

SET,NEXT    ! NEXT DATA SET
*ENDDO      ! MOVE TO NEXT DATA SET (STORED TIME STEP)

! DONE EXPORTING

/eof
