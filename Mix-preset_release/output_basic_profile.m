function output_basic_profile(u,yc)

GG=zeros(length(yc)+2*4,1);
GG(1)=u(1);
GG(2)=u(1);
GG(3)=u(1);
GG(4)=u(1);
GG(5:end-4)=u(:);
GG(end-3)=u(end);
GG(end-2)=u(end);
GG(end-1)=u(end);
GG(end-0)=u(end);

fid=fopen('PROFILE-BL-CENTER.DAT','wb');
if(fid>0)
    fwrite(fid,size(GG,1),'integer*4');
    for j=1:size(GG,1)
        fwrite(fid,GG,'real*8');
    end
end
fclose(fid);

Ny=size(GG,1)

end