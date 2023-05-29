function decision = nnForwardPropagation(I, h_1, h_2, o_1, o_2) 

w_h    = [h_1 h_2]; 
w_o    = [o_1 o_2];

I      = [I(1) I(2) I(3) 1];

%Hidden layer
H   = I*w_h;   

type = 'tanh'; 
H_a = nnActivation(H,type);  

H_a = [H_a(1) H_a(2) 1];

%Output layer
O    = H_a*w_o; 

type = 'sigmoid'; 
O_a = nnActivation(O,type);

decision = [0 0]; 

[dec,id] = max(O_a);

decision(id) = 1; 


end