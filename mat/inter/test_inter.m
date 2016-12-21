load data1.mat
method = 7;
inter_par= interpolateparametarization(xi,yi,method,[])
x = [0.5;0.5]
y = interpolate_val([0.5;0.5], inter_par)