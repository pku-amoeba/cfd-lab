function output_eigenfun(N,eigfun,ph,flag)


FF=zeros(5,N); 
for var=1:4
    FF(var,2:N-1)=eigfun(var:5:end,flag);
    FF(var,1)=FF(var,2);
    FF(var,end)=FF(var,end-1);
end

FF(5,2:N-1)=ph;
FF(5,1)=FF(5,2);
FF(5,end)=FF(5,end-1);

GG=zeros(5,N+2*4);
for var=1:5
    GG(var,1)=FF(var,1);
    GG(var,2)=FF(var,1);
    GG(var,3)=FF(var,1);
    GG(var,4)=FF(var,1);
    GG(var,5:end-4)=FF(var,:);
    GG(var,end-3)=FF(var,end);
    GG(var,end-2)=FF(var,end);
    GG(var,end-1)=FF(var,end);
    GG(var,end-0)=FF(var,end);
end

fid=fopen('EIGENFUN.DAT','wb');
if(fid>0)
    fwrite(fid,size(GG),'integer*4');
    for var=1:size(GG,1)
        for j=1:size(GG,2)
            fwrite(fid,[real(GG(var,j)),imag(GG(var,j))],'real*8');
        end
    end
end
fclose(fid);

% fid1=fopen('EIGENFUN_UR.DAT');
% if(fid1>0)
%     xx=fread(fid1,1,'integer*4');
%     xxx=fread(fid1,1,'integer*4');
%     for j=1:N+2*4
%         UR(j)=fread(fid1,[1 1],'real*8');
%     end
% end
% fclose(fid1);


% plot(imag(eigenfun_beta(4,:)))
% hold on 
% plot(-imag(eigenfun_mbeta(4,:)))

end