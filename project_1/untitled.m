num_simulations = 100000;  
target_position = [0, 0, 0];  
gps_error = [0, 0.025];  
lidar_error = [0, 0.1];  
euler_angle_error = [0, 0.01745];  
runs = zeros(1, num_simulations); 
for i = 1:num_simulations
    runs(i) = simulate_hit_or_miss(target_position, gps_error, lidar_error, euler_angle_error);
end

histogram(runs, 50)
title("Histogram of Landing Distances from Center of Zone")
xlabel("Feet Away From Center")
set(gca,'ytick',[])
xlim([-2 18])

function distance = simulate_hit_or_miss(target_position, gps_error, lidar_error, euler_angle_error)
   
    offset_x_y = normrnd(0, 0.82, [1, 2]);
    error_x_y = normrnd([gps_error(1), gps_error(1)], [gps_error(2), gps_error(2)]);
    error_x_y = error_x_y + offset_x_y;
    error_z = normrnd(lidar_error(1), lidar_error(2));
    error_pitch_roll_yaw = normrnd([euler_angle_error(1), euler_angle_error(1), euler_angle_error(1)], ...
                                   [euler_angle_error(2), euler_angle_error(2), euler_angle_error(2)]);
    
    % this sucks
    error_x_y(1) = error_x_y(1) + rad2deg(error_pitch_roll_yaw(3)); 
    error_x_y(2) = error_x_y(2) + rad2deg(error_pitch_roll_yaw(1)); 

    landing_position = target_position + [error_x_y, error_z];
    distance_from_target = sqrt(sum((landing_position - target_position).^2)) * 3.28;
    distance = distance_from_target;
end
