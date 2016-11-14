function fposchar(fp, ch, n)

if ( nargin < 3 )
    n = 1;
end

for i=1:n
    jnk=0;
    while(jnk~=ch)
        jnk=fscanf(fp,'%1c',[1 1]);
    end
end
