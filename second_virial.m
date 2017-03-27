function [ B ] = second_virial( T,eps,sig )
%This function evaluates the second virial coefficient at a specific
%temperature for a given epsilon (in K) and sigma (in Ang)

sum = 0;

for n = 0:30
   
    product = (2^((2*n+1)/2)) / (4*factorial(n)) .* (eps./T).^((2*n+1)/4) .* gamma((2*n-1)/4);
    sum = sum + product;
    
end

B = sum .* (-2*pi*sig.^3) / 3 / 1660.5778811026237; % The last term converts it to kmol/m^3

end

