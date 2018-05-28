function pbc(x,N)
	n=maximum(size(x));
	
	Y=zeros(Int8,n,1);
	for i=1:n
		y=mod(x[i],N);
		if(y==0) 
			Y[i]=N;
		else 
			Y[i]=y;
		end
	end 
	return Y;
end