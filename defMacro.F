
defExport
! Created by Basil J. Paudel            5/13/2021
! This macro will export the nodal coordinates and the deformation results
!         to a CSV file
! Usage:
! First place the defMacro.F in the Ansys working directory
!         *ulib,defMacro,F,
!         *use, defExport

! DEFINE OUTPUT PARAMETERS
fname='Results'

! DO NOT CHANGE BELOW THIS LINE

/post1
set,last

*MSG,UI 
Exporting deformation results to Results.csv

*GET, maxNode, NODE,0,NUM,MAX
*GET, minNode, NODE,0,NUM,MIN

nodeNum=1
allsel

! FIRST SETUP A TABLE TO STORE DATA
*DIM, XPT,TABLE, maxNode,3

*do,nodeNum,minNode,maxNode

! GET NODAL COORDINATES
xx=NX(nodeNum)
yy=NY(nodeNum)
zz=NZ(nodeNum)

! GET NODAL DISPLACEMENTS
sux=UX(nodeNum)
suy=UY(nodeNum)
suz=UZ(nodeNum)

*SET,XPT(nodeNum,0),nodeNum
*SET,XPT(nodeNum,1),xx+sux
*SET,XPT(nodeNum,2),yy+suy
*SET,XPT(nodeNum,3),zz+suz


*enddo

! NOW WRITE TO FILE

*CFOPEN,fname,'csv'           ! Open file for writing
*VWRITE,'Node','X','Y','Z'
(4A10)
*VWRITE, XPT(1,0),XPT(1,1),XPT(1,2),XPT(1,3)
(1f9.1, 3f12.6)

*CFCLOSE  ! Close the File

! DONE EXPORTING

/eof
