load ff.dat;

[n,m]=size(ff)
number=ones(n,1);
ncount=1;
for j=1:n
  for k=1:n
    if(j~=k&&ff(j,2)==ff(k,2)&&ff(j,3)==ff(k,3))
      [ncount j k]
      number(ncount)=j;
      ncount=ncount+1;
    end
  end
end

