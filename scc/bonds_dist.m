#julia program to read cube.xyz and write bonds in lammps format
#bonds are added between the nearest and the second nearest neighbors

input=ARGS[1];

f=readdlm(input);

D=float(f[3:end,2:end]);
N=convert(Int64,float(ARGS[2]));
natom=size(D,1);

ff=open("bonds2.lmp","w");

bcount=0;
a=4.08;
bd1=a*a;
bd2=bd1*2;

Lx=N*a;
Ly=N*a;
Lz=N*a;

for p=1:natom
	for n=p:natom
		
		dx=min(abs(D[p,1]-D[n,1]),Lx-(abs(D[p,1]-D[n,1]))); 
		dy=min(abs(D[p,2]-D[n,2]),Ly-(abs(D[p,2]-D[n,2]))); 
		dz=min(abs(D[p,3]-D[n,3]),Lz-(abs(D[p,3]-D[n,3]))); 
		
		b=2;
		
		if(abs(dx^2+dy^2+dz^2-bd1)<=0.1) 
			b=1;
			bcount=bcount+1;
			@printf(ff,"%i  %i  %i  %i\n",bcount,b,p,n)
			if(p==1)
			println(sqrt(dx^2+dy^2+dz^2)," ",n); 
			end
		end
		
		#if(abs(dx^2+dy^2+dz^2-bd2)<=0.1) 
		#	bcount=bcount+1;
		#	@printf(ff,"%i  %i  %i  %i\n",bcount,b,p,n);
		#	if(p==1)
		#	println(sqrt(dx^2+dy^2+dz^2)," ",n);
		#	end
		#end
		
	end
end

close(ff)


