import numpy as np
import matplotlib.pyplot as plt

m = 0.145 #kg, mass of baseball
g = np.array([0, -9.8]) #m/s^2, accel of grav at surface
rho = 1.225 #kg/m^3, density of air at STP
Cd = 0.47 #[unitless], drag coefficient of sphere
A = 0.00414 #m^2, cross-sectional area of baseball

def analytic_projectile_nodrag(t,
                               initial_position=None,
                               initial_velocity=None):
    '''t = array of times
    initial_position = inition position coordinates
    initial_velocity = initial velocity vector'''
    times = t[:,np.newaxis]
    return initial_position +\
           initial_velocity*times + \
           g*0.5*times**2

def deriv_noDrag(v):
    '''v = anything (or velocity vector)'''
    return g

def deriv_withDrag(v, m, g, rho, Cd, A):
    '''v = velocity vector,
    m = mass
    g = gravitational acceleration vector
    rho = density of medium
    Cd = drag coefficient
    A = cross-sectional area'''
    vMag = (v[0]**2+v[1]**2)**0.5
    drag = -0.5*rho*(vMag**2)*Cd*A*(v/vMag)
    return (m*g+drag)/m

def rk4(t,
        timestep = None,
        initial_x = None,
        initial_v = None,
        deriv = None,
        deriv_params = None):
    '''t = array of times
    timestep = time between time values
    initial_x = inition position coordinates
    initial_v = initial velocity vector
    deriv = derivative function
    deriv_params = dict of derivative function parameters'''
    times = t[:,np.newaxis]
    values = [initial_x]
    x = initial_x
    v = initial_v
    for time in times:
        k1 = timestep*v
        l1 = timestep*deriv(v, **deriv_params)
        k2 = timestep*(v+l1/2)
        l2 = timestep*deriv(v+l1/2, **deriv_params)
        k3 = timestep*(v+l2/2)
        l3 = timestep*deriv(v+l2/2, **deriv_params)
        k4 = timestep*(v+l3)
        l4 = timestep*deriv(v+l3, **deriv_params)
        x = x+(1/6)*(k1+2*k2+2*k3+k4)
        values.append(x)
        v = v+(1/6)*(l1+2*l2+2*l3+l4)
    return np.array(values)

if __name__ == "__main__":
    times = np.arange(0, 10, 0.01)

    #calculates baseball's path analytically w/o drag at 45 degrees and 10 m/s
    angle = 45*np.pi/180 #degrees
    speed = 10 #m/s
    x0 = np.array([0,2])
    v0 = speed*np.array([np.cos(angle), np.sin(angle)])
    positions = analytic_projectile_nodrag(times,initial_position=x0,initial_velocity=v0)
    
    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's analytic path w/o drag
    plt.plot(x_pos, y_pos)

    #calculates baseball's path using rk4 w/o drag at 45 degrees and 10 m/s
    angle = 45*np.pi/180 #degrees
    speed = 10 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_noDrag,
                    deriv_params={})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions[1::3]:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's rk4 path w/o drag
    plt.scatter(x_pos, y_pos, color='orange', s=15)
    
    plt.xlabel('Horizontal Distance (m)')
    plt.ylabel('Vertical Distance (m)')
    plt.title("Baseball's Path w/o Drag")
    plt.legend(['Calculated Analytically','Calulated w/RK4'])
    plt.savefig('Baseball_Path_NoDrag.png')
    plt.show()
    plt.clf

    #calculates baseball's path using rk4 with drag at 30 degrees and 85 mph
    angle = 30*np.pi/180 #degrees
    speed = 37.9984 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_withDrag,
                    deriv_params={'m':m, 'g':g, 'rho':rho, 'Cd':Cd, 'A':A})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's path
    plt.plot(x_pos, y_pos, c='blue', alpha=0.5)

    #calculates baseball's path using rk4 with drag at 45 degrees and 85 mph
    angle = 45*np.pi/180 #degrees
    speed = 37.9984 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_withDrag,
                    deriv_params={'m':m, 'g':g, 'rho':rho, 'Cd':Cd, 'A':A})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's path
    plt.plot(x_pos, y_pos, c='red', alpha=0.5)

    #calculates baseball's path using rk4 with drag at 60 degrees and 85 mph
    angle = 60*np.pi/180 #degrees
    speed = 37.9984 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_withDrag,
                    deriv_params={'m':m, 'g':g, 'rho':rho, 'Cd':Cd, 'A':A})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's path
    plt.plot(x_pos, y_pos, c='gold', alpha=0.5)

    #calculates baseball's path using rk4 w/o drag at 30 degrees and 85 mph
    angle = 30*np.pi/180 #degrees
    speed = 37.9984 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_noDrag,
                    deriv_params={})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's path
    plt.plot(x_pos, y_pos, color='blue', linestyle='--')

    #calculates baseball's path using rk4 w/o drag at 45 degrees and 85 mph
    angle = 45*np.pi/180 #degrees
    speed = 37.9984 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_noDrag,
                    deriv_params={})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's path
    plt.plot(x_pos, y_pos, color='red', linestyle='--')

    #calculates baseball's path using rk4 w/o drag at 60 degrees and 85 mph
    angle = 60*np.pi/180 #degrees
    speed = 37.9984 #m/s
    v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
    positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_noDrag,
                    deriv_params={})

    #separates baseball's path into x and y components for graphing
    x_pos = []
    y_pos = []
    for position in positions:
        if position[1] >= 0:
            x_pos.append(position[0])
            y_pos.append(position[1])

    #graphs baseball's path
    plt.plot(x_pos, y_pos, color='gold', linestyle='--')
    
    plt.xlabel('Horizontal Distance (m)')
    plt.ylabel('Vertical Distance (m)')
    plt.title("Baseball's Path With and Without Drag")
    plt.legend(['30 Degrees w/Drag','45 Degrees w/Drag','60 Degrees w/Drag',
                '30 Degrees w/o Drag','45 Degrees w/o Drag','60 Degrees w/o Drag'])
    plt.savefig('Baseball_Path_With_And_Without_Drag.png')
    plt.show()
    plt.clf
    
    #calculates max distance a baseball can be thrown at 108 mph
    maxDistance = 0
    angles = np.linspace(30*np.pi/180,60*np.pi/180,30) #angles at which to test distance
    for angle in angles:
            speed = 48.28032 #m/s
            v0 = speed*np.array([np.cos(angle), np.sin(angle)]) #baseball's initial velocity
            positions = rk4(times,timestep=0.01,initial_x=x0,initial_v=v0,deriv=deriv_withDrag,
                            deriv_params={'m':m, 'g':g, 'rho':rho, 'Cd':Cd, 'A':A})
            for position in positions: #tests if max distance is greater than at previous angles
                    if position[1] >= 0:
                            if position[0] > maxDistance:
                                    maxDistance = position[0]
    print('The maximum distance a baseball can be thrown at 108 mph is', maxDistance, 'meters')
    #found answer is 107.375
