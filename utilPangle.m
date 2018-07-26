function pangle = utilPangle(th)
    if(th<0), 
        pangle=2*pi+th; 
    else
        pangle = th;
    end;  
end