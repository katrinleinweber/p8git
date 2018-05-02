#julia program to read cube.xyz and write bonds in lammps format
#bonds are added between the nearest and the second nearest neighbors

input=ARGS[1];

f=readdlm(input);

D=f[3:end,2:end];
N=convert(Int64,float(ARGS[2]));
natom=size(D,1);

array1=[1 0 0;
		0 1 0;
		0 0 1;
		1 1 0;
		0 1 1;
		1 0 1;
		-1 1 0;
		1 0 -1;
		0 -1 1];#changed

		array2=-array1;
myarray=[array1;array2];

ff=open("bonds.lmp","w");

bonds=zeros(natom,size(myarray,1));
bcount=0;
vector=[1;N;N*N];

for p=1:N
	for n=1:N
		for m=1:N
			
			for j=1:size(myarray,1)
			
				ivec=[m;n-1;p-1];
				i=sum(ivec.*vector);
				
				
				
				lxyz=[myarray[j,1];myarray[j,2];myarray[j,3]]+ivec;
				
				#Apply pbc forward
				
				if(lxyz[1]==N+1)
					lxyz[1]=1;
				end
				
				if(lxyz[2]==N)
					lxyz[2]=0;
				end
				
				if(lxyz[3]==N)
					lxyz[3]=0;
				end
				
				## Apply pbc backward
				if(lxyz[1]==0)
					lxyz[1]=N;
				end
				
				if(lxyz[2]==-1)
					lxyz[2]=N-1;
				end
				
				if(lxyz[3]==-1)
					lxyz[3]=N-1;
				end
				
				#println("i=",i," ",lxyz," ",ivec )
				bonds[i,j]=sum(vector.*lxyz);
				
				if(bonds[i,j]>i)
					bcount=bcount+1;
					if(sum(myarray[j,:].*myarray[j,:])>1) b=2; else b=1;end
					@printf(ff,"%i  %i  %i  %i\n",bcount,b,i,bonds[i,j]);
					
				end
			end
			
		end
	end
end

close(ff)


