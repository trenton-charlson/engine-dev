function dydt = biteme(t,y)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2)=-48cos(y(1))-.3(y(2)
end