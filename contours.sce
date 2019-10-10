clear
clc
close


x=[1:1:7]
y=x

scf(1)
filename="p.csv"        
p = csvRead(filename)
pp=p'
contourf(x,y,p,20)
//contourf(x,y,T,8,3:8,"111"," ",rect=[0,0,6,4])
title("p")


scf(2)
filename="u.csv"        
u = csvRead(filename)
uu=u'
contour2d(x,y,u,20)
//contourf(x,y,T,8,3:8,"111"," ",rect=[0,0,6,4])
title("u")


scf(3)
filename="v.csv"        
v = csvRead(filename)
vv=v'
contour2d(x,y,v,20)
//contourf(x,y,T,8,3:8,"111"," ",rect=[0,0,6,4])
title("v")
