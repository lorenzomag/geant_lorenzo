      REAL E_int, E_rec, E_g(15), E_gp(15), cost_g(15), phi_g(15),
     +     cost_gp(15), cost_r, cosp_r, x_r, y_r, z_r, thet_r,
     +	   xstop, ystop, zstop, xint, yint, zint,x,y,xp,yp,z,xtest(10),
     +     ytest(10),
     +	   etest(10), binitenerg, half_thickness, stop_power,beamtof
      INTEGER Nodec, react, recdet, dsssdpos
      COMMON /HISTORY/ E_int, E_rec, E_g, E_gp, cost_g, phi_g,
     +     cost_gp, cost_r, cosp_r, Nodec, react, recdet, x_r, y_r, z_r,
     +	   thet_r, xstop, ystop, zstop, xint, yint, zint,x,y,xp,yp,z,xtest,
     +	   ytest, etest, dsssdpos,binitenerg, half_thickness, stop_power,beamtof
