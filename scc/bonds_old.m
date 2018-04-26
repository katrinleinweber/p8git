#julia program to read cube.xyz and write bonds in lammps format
#bonds are added between the nearest and the second nearest neighbors

input=ARGS[1];

f=readdlm(input);

D=f[3:end,2:end];
N=convert(Int64,float(ARGS[2]));
natom=size(D,1);

myarray=[1 0 0;
		0 1 0;
		0 0 1;
		1 1 0;
		0 1 1;
		1 0 1;
		-1 1 0;
		-1 0 1;
		0 1 -1];#changed

ff=open("bonds.lmp","w");

bonds=zeros(natom,size(myarray,1));
bcount=0;
jx_up=[1 4 6 ];jy_up=[2 4 5 7 9];jz_up=[3 5 6 8];
jx_dn=[7 8];jy_dn=[];jz_dn=[9];
vector=[1;N;N*N];

for p=1:N
	for n=1:N
		for m=1:N
			
			for j=1:size(myarray,1)
			
				ivec=[m;n-1;p-1];
				i=sum(ivec.*vector);
				lxyz=[myarray[j,1];myarray[j,2];myarray[j,3]]+ivec;
				#Apply pbc forward
				
				if((m==N)&&(!isempty(find(jx_up.==j))))
					ivec[1]=0;
				end
				
				if((n==N)&&(!isempty(find(jy_up.==j))))
					ivec[2]=-1;
				end
				
				if((p==N)&&(!isempty(find(jz_up.==j))))
					ivec[3]=-1;
				end
				
				## Apply pbc backward
				if((m==1)&&(!isempty(find(jx_dn.==j))))
					ivec[1]=N;
				end
				
				if((n==1)&&(!isempty(find(jy_dn.==j))))
					ivec[2]=N-1;
				end
				
				if((p==1)&&(!isempty(find(jz_dn.==j))))
					ivec[3]=N-1;
				end
				
				#lxyz=[myarray[j,1];myarray[j,2];myarray[j,3]]+ivec;
				println("i=",i," ",lxyz," ",ivec )
				bonds[i,j]=sum(vector.*lxyz);
				bcount=bcount+1;
				@printf(ff,"%i  1  %i  %i\n",bcount,i,bonds[i,j]);
			end
			
		end
	end
end

close(ff)


