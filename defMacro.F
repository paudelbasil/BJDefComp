
defExport
! Created by Basil J. Paudel            5/13/2021
! This macro will export the nodal coordinates and the deformation results
!         to a CSV file
! Usage:
!         *ulib,defMacro,mac,
!         *use, defMacro


/post1
set,last

fname= arg1
*if,fname,eq,0,then
   fname='Results'
*endif

*MSG,UI 
Exporting deformation results to Results.csv

*GET, maxNode, NODE,0,NUM,MAX
*GET, minNode, NODE,0,NUM,MIN

nodeNum=1
allsel

! FIRST SETUP A TABLE
*DIM, XPT,TABLE, maxNode,6


*do,nodeNum,minNode,maxNode
xx=NX(nodeNum)
yy=NY(nodeNum)
zz=NZ(nodeNum)

sux=UX(nodeNum)
suy=UY(nodeNum)
suz=UZ(nodeNum)

*SET,XPT(nodeNum,0),nodeNum
*SET,XPT(nodeNum,1),xx
*SET,XPT(nodeNum,2),yy
*SET,XPT(nodeNum,3),zz
*SET,XPT(nodeNum,4),sux
*SET,XPT(nodeNum,5),suy
*SET,XPT(nodeNum,6),suz


*enddo

! NOT A PREFERRED OPTION, USE CFWRITE
! /output, 'Results','csv','.'
*CFOPEN,fname,'csv'           ! Open file for writing
*VWRITE,'Node','X','Y','Z','UX','UY', 'UZ'
(7A10)
*VWRITE, XPT(1,0),XPT(1,1),XPT(1,2),XPT(1,3),XPT(1,4),XPT(1,5),XPT(1,6)
(1f9.1, 6f10.6)
! %G %G %G %G %G %G %G

!/output
*CFCLOSE  ! Close the File


/eof
