function output_mesh(nx,nx1,ny,nz,Lx,Lx1,Ly,Lz,dy1,aa0,beta0)

%nx=1421;
%Lx=11.58*20;
%Lx1=11.58*16;
%nx1=1281;

% Streamwise grid
[beta,xx]=get_mesh_stream(nx,nx1,Lx,Lx1,beta0);

% Get normal direction grid
[aa,yy]=get_mesh_normal(ny,Ly,dy1,aa0);

% Get spanwise direction grid
dz=Lz/(nz-1);
for k=1:nz
    zz(k)=(k-1)*dz -Lz/2;
end

% Generate 3D mesh
for k=1:nz
    for j=1:ny
        for i=1:nx
            x(i,j,k)=xx(i);
            y(i,j,k)=yy(j);
            z(i,j,k)=zz(k);
        end
    end
end

% Output mesh to fortran readable format
fid=fopen('MESH.IN','wb');
if(fid>0)
    fwrite(fid,1,'integer*4');
    fwrite(fid,[nx,ny,nz],'integer*4');
    for i=1:nx
        for j=1:ny
            for k=1:nz
                fwrite(fid,x(i,j,k),'real*8');
            end
        end
    end
    
    for i=1:nx
        for j=1:ny
            for k=1:nz
                fwrite(fid,y(i,j,k),'real*8');
            end
        end
    end    
    
    for i=1:nx
        for j=1:ny
            for k=1:nz
                fwrite(fid,z(i,j,k),'real*8');
            end
        end
    end
    
end
fclose(fid);

end