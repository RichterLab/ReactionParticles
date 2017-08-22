%generate_fields.m

num_fields = 1; %number of velocity fields to generate

%specify domain size
Lx = 10;
Ly = 1;

%specify permeability field parameters
var = 1;     %variance of the lognormal K field
k_range = 5; %range of wavenumbers used in generation of permeability fields

%specify velocity field parameters
mean_velocity = 0.001; %mean velocity in x direction (mean u)


avg_vel = zeros(num_fields,2); %column 1 = mean(u), column 2 = mean(v)

for ii = 1:num_fields
    ii
    filenumber = ii;
    
    [K]=permeability(Lx,Ly,var,k_range,filenumber);  %generate permeability K
    [u,v,P]=velocity_function(Lx,Ly,mean_velocity,filenumber); %generate velocity fields
    
    avg_vel(ii,1) = mean(mean(u));
    avg_vel(ii,2) = mean(mean(v));
end

avg_u = mean(avg_vel(:,1))
avg_v = mean(avg_vel(:,2))
