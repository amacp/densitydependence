clc

k1 = 5;
k2 = 5;
b1 = 0.2;
b2 = 0.6;
d = 0.05;
nmax = max(k1,k2);
S = sparse((nmax+1)*(nmax+1), (nmax+1)*(nmax+1));

for i1 = 0:nmax
    for i2 = 0:nmax
        c = e(i1,i2,nmax);
        a1 = 0; a2 = 0; a3 = 0; a4 = 0;
        for j1 = 0:nmax
            for j2 = 0:nmax
                r = e(j1,j2,nmax);
                %Scenario 1
                if(j1 == i1-1 && j2 == i2)
                   if (S(r,c) ~= 0)
                       S(r,c)= S(r,c) + d*i1;
                   else
                       S(r,c)= d*i1;
                   end
                   a1 = a1+S(r,c);
                end
                %Scenario 2                
                if(j1 == i1 && j2 == i2-1)
                   if (S(r,c) ~= 0)
                       S(r,c)= S(r,c) + d*i2;
                   else
                       S(r,c)= d*i2;
                   end
                   a1 = a1+ S(r,c);
                end 
                %Scenario 3
                if(j1 == i1+1 && j2 == i2)
                   if (S(r,c) ~= 0)
                       if(i1+i2<k1)
                           S(r,c)= S(r,c) + b1*i1*(1-(i1+i2)/k1);
                       else
                           S(r,c)= S(r,c) + 0;
                       end
                   else
                       if(i1+i2<k1)
                           S(r,c)= b1*i1*(1-(i1+i2)/k1);
                       else
                           S(r,c)=0;
                       end
                   end
                   a1 = a1+ S(r,c);
                end
                %Scenario 4
                if(j1 == i1 && j2 == i2+1)
                   if (S(r,c) ~= 0)
                       if(i1+i2<k2)
                           S(r,c)= S(r,c) + b2*i2*(1-(i1+i2)/k2);
                       else
                           S(r,c)= S(r,c) + 0;
                       end
                   else
                       if(i1+i2<k2)
                           S(r,c)= b2*i2*(1-(i1+i2)/k2);
                       else
                           S(r,c)=0;
                       end
                   end
                   a1 = a1+S(r,c);
                end
                      
            end
        end
        S(c,c) = -a1;
    end
end

mc = dtmc(S)


function out = e(i1,i2,nmax)
    out = i1*(nmax+1)+i2+1;
end