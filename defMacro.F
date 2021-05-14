
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
*DIM, ExportedResults,TABLE, maxNode,6


*do,nodeNum,minNode,maxNode
xx=NX(nodeNum)
yy=NY(nodeNum)
zz=NZ(nodeNum)

sux=UX(nodeNum)
suy=UY(nodeNum)
suz=UZ(nodeNum)

*SET,ExportedResults(nodeNum,0),nodeNum
*SET,ExportedResults(nodeNum,1),xx
*SET,ExportedResults(nodeNum,2),yy
*SET,ExportedResults(nodeNum,3),zz
*SET,ExportedResults(nodeNum,4),sux
*SET,ExportedResults(nodeNum,5),suy
*SET,ExportedResults(nodeNum,6),suz


*enddo

! NOT A PREFERRED OPTION, USE CFWRITE
! /output, 'Results','csv','.'
*CFOPEN,fname,'csv'           ! Open file for writing
*VWRITE,'Node','X','Y','Z','UX','UY', 'UZ'
(7A10)
*VWRITE, ExportedResults(1,0),ExportedResults(1,1),ExportedResults(1,2),ExportedResults(1,3),ExportedResults(1,4),ExportedResults(1,5),ExportedResults(1,6)
(1f9.1, 6f10.6)
! %G %G %G %G %G %G %G

!/output
*CFCLOSE  ! Close the File


/eof
