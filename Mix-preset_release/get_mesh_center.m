function yc=get_mesh_center(y,ny)

for j=1:ny-1
    
    yc(j)=0.5*(y(j+1)+y(j));

end

end