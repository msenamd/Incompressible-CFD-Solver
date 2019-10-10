clear
clc
close all

n=30;
      x=[1:1:n];
      y=x;
      
      
      pp=load("p.csv");
      uu_f=load("u_f.csv");
      vv_f=load("v_f.csv");
      uu_b=load("u_b.csv");
      vv_b=load("v_b.csv");
      
      
      p=pp'
      u_f=uu_f'
      v_f=vv_f'
      u_b=uu_b'
      v_b=vv_b'
      
      figure(1)
      contourf(x,y,p,20)
      h=colorbar ();
      colormap ('jet');
      title('p');
      grid on;
      
      figure(2)
      contourf(x,y,u_f,20)
      h=colorbar ();
      colormap ('jet');
      title('u_f');
      grid on;

            
      figure(3)
      contourf(x,y,v_f,20)
      h=colorbar ();
      colormap ('jet');
      title('v_f');
      grid on;

      figure(4)
      contourf(x,y,u_b,6)
      h=colorbar ();
      colormap ('jet');
      title('u_b');
      grid on;

            
      figure(5)
      contourf(x,y,v_b,6)
      h=colorbar ();
      colormap ('jet');
      title('v_b');
      grid on;


  

