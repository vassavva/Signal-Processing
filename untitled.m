m=4;
x=1;
input_x = zeros(1,m+1) ; 
for i = 1:3
    for n = 1:m
            temp = input_x;
            input_x(n+1) = temp(n);
    end
    
    input_x(1) = x;
    x=x+1;
end