#julia program to read cube.xyz and write bonds in lammps format
#bonds are added between the nearest and the second nearest neighbors

include("pbc.m")

input=ARGS[1];

f=readdlm(input);

D=f[3:end,2:end];
N=convert(Int64,float(ARGS[2]));

natom=size(D,1);
#M=reshape(D,N,N,N);

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

ff=open("bonds1.lmp","w");

bonds=zeros(natom,size(myarray,1));
bcount=0;


for z=1:N
	for y=1:N
		for x=1:N
			i=x+(y-1)*N+(z-1)*N*N;
			
			for j=1:18	
				temp=[x y z]+myarray[j,:];
				
				##PBC

				temp=pbc(temp,N);
				
				bonds[i,j]=temp[1]+(temp[2]-1)*N+(temp[3]-1)*N*N;	
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


