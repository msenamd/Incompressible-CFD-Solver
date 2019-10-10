clear
clc
close all

      n=7;

      x=[2:1:n-1];
      y=x;
      
      
      pp=load("p.csv");
      uu=load("u.csv");
      vv=load("v.csv");
      psii=load("stream.csv");
      vorr=load("vorticity.csv");
      
      
      p=pp';
      u=uu';
      v=vv';
      psi=psii';
      vor=vorr';

      U_mag=sqrt(u.^2+v.^2);
            
      
      figure(1,"position",[0,0,800,600])
      contourf(x,y,p,12)
      h=colorbar ();
      colormap ("jet");
      xlabel("X","fontsize",20)
      ylabel("Y","fontsize",20)
      title("P","fontsize",20)
      set(gca, "fontsize", 20)
      set(h, "fontsize", 20)
      axis equal
      
      figure(2,"position",[0,0,800,600])
      contourf(x,y,u,12)
      h=colorbar ();
      colormap ("jet");
      xlabel("X","fontsize",20)
      ylabel("Y","fontsize",20)
      title("u","fontsize",20)      
      set(gca, "fontsize", 20)
      set(h, "fontsize", 20)
      axis equal

            
      figure(3,"position",[0,0,800,600])
      contourf(x,y,v,12)
      h=colorbar ();
      colormap ("jet");
      xlabel("X","fontsize",20)
      ylabel("Y","fontsize",20)
      title("v","fontsize",20)     
      set(gca, "fontsize", 20)
      set(h, "fontsize", 20)
      axis equal

      figure(4,"position",[0,0,800,600])
      contourf(x,y,U_mag,12)
      h=colorbar ();
      colormap ("jet");
      xlabel("X","fontsize",20)
      ylabel("Y","fontsize",20)
      title("Umag","fontsize",20)
      set(gca, "fontsize", 20)
      set(h, "fontsize", 20)
      axis equal
      
      figure(5,"position",[0,0,800,600])
      contourf(x,y,vor,20)
      h=colorbar ();
      colormap ("jet");
      xlabel("X","fontsize",20)
      ylabel("Y","fontsize",20)
      title("vorticity","fontsize",20)
      set(gca, "fontsize", 20)
      set(h, "fontsize", 20)
      axis equal
      
      figure(6,"position",[0,0,800,600])
      contourf(x,y,psi,20)
      h=colorbar ();
      colormap ("jet");
      xlabel("X","fontsize",20)
      ylabel("Y","fontsize",20)
      title("stream function","fontsize",20)
      set(gca, "fontsize", 20)
      set(h, "fontsize", 20)
      axis equal
     