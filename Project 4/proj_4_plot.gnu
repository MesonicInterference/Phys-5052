## actual position and velocity
#set autoscale y
#plot "vel.dat",  "pos.dat"

## calculated velocities
set autoscale x
#set autoscale y
plot "rk4_vel.dat"

## all velocities
#set autoscale y
#plot "mod_euler_vel.dat", "rk4_vel.dat", "vel.dat"

## velocity errors
#set logscale y
#plot "mod_euler_vel_err.dat", "rk4_vel_err.dat"

## calculated positions
#set autoscale y
#plot "mod_euler_pos.dat", "rk4_pos.dat"

## all positions
#set autoscale y
#plot "mod_euler_pos.dat", "rk4_pos.dat", "pos.dat"

## position errors
#set logscale y
#plot "rk4_pos_err.dat"

## modified Euler errors
#set logscale y
#plot "mod_euler_vel_err.dat",  "mod_euler_pos_err.dat"

## RK4 errors
#set logscale y
#plot "rk4_vel_err.dat",  "rk4_pos_err.dat"

pause -1