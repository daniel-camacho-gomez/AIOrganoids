function activation_output = nnActivation(x,type) 

if strcmp(type,'sigmoid')
activation_output = 1./(1+exp(-x));
end

if strcmp(type,'rel') 
activation_output = max(0,x(:))';  
end

if strcmp(type,'tanh') 
activation_output = tanh(x);  
end


end

