clear
clc
close all

n=7;
      x=[1:1:n];
      y=x;
      
      
      p=load("R.csv");
      u=load("F_f.csv");
      v=load("F_b.csv");
      
      R=p'
      F_f=u'
      F_b=v'
      
      figure(1)
      contourf(x,y,R,6)
      h=colorbar ();
      colormap ('jet');
      title('R');
      grid on;
      
      figure(2)
      contourf(x,y,F_f,6)
      h=colorbar ();
      colormap ('jet');
      title('F_f');
      grid on;

            
      figure(3)
      contourf(x,y,F_b,6)
      h=colorbar ();
      colormap ('jet');
      title('F_b');
      grid on;



  

