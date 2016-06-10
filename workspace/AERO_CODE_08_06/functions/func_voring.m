function ind = func_voring(X_CP,Y_CP,Z_CP,vortex_ring,GAMMA,ONOFF)
       
    if vortex_ring.Z_SIGN == -1
       GAMMA = - GAMMA;
    end
    
    [u1,v1,w1] = func_vortexl(X_CP,Y_CP,Z_CP,vortex_ring.X_1,vortex_ring.Y_1,...
        vortex_ring.Z_1,vortex_ring.X_2,vortex_ring.Y_2,vortex_ring.Z_2,GAMMA);
    
    [u2,v2,w2] = func_vortexl(X_CP,Y_CP,Z_CP,vortex_ring.X_2,vortex_ring.Y_2,...
        vortex_ring.Z_2,vortex_ring.X_4,vortex_ring.Y_4,vortex_ring.Z_4,GAMMA);
    
    [u3,v3,w3] = func_vortexl(X_CP,Y_CP,Z_CP,vortex_ring.X_4,vortex_ring.Y_4,...
        vortex_ring.Z_4,vortex_ring.X_3,vortex_ring.Y_3,vortex_ring.Z_3,GAMMA);
    
    [u4,v4,w4] = func_vortexl(X_CP,Y_CP,Z_CP,vortex_ring.X_3,vortex_ring.Y_3,...
        vortex_ring.Z_3,vortex_ring.X_1,vortex_ring.Y_1,vortex_ring.Z_1,GAMMA);     
    
switch ONOFF
    case  0  
ind = [u1+u2+u3+u4,v1+v2+v3+v4,w1+w2+w3+w4];
    case  1
ind = [u2+u4,v2+v4,w2+w4];
    case  2
ind = [u2+u3+u4,v2+v3+v4,w2+w3+w4];
        
end

end