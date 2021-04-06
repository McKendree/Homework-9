import math
import numpy as np
import matplotlib.pyplot as plt

def analytic_solution(t, K=1, P0=1, r=0.01):
    """K = maximum carrying capacity of ecosystem
       P0 = initial population at time t=0
       r = population growth rate (eg, 1% = 0.01
       t = time """
    populations = []
    for time in t:
        numerator = K*P0*math.e**(r*time)
        denominator = K+P0*(math.e**(r*time)-1)
        populations.append(numerator/denominator)
    return np.array(populations)

def dP_dt(P, K=10, r=0.01):
    '''P = population
    K = population carrying capacity
    r = population growth rate'''
    return r*P*(1-P/K)


def forward_euler(timestep = None,
                  max_time = None,
                  initial_time = None,
                  initial_val = None,
                  deriv = None, #pass dP_dt
                  deriv_params = None):
    '''timestep = time between time values
    max_time = the last time calculated for
    initial_time = the first time calculated for
    initial_val = initial value
    deriv = derivative function
    deriv_params = derivative function parameters'''
    times = np.linspace(initial_time, max_time, timestep)
    values = [initial_val]
    val = initial_val
    for time in times:
        val = val+timestep*deriv(val, **deriv_params)
        values.append(val)
    return np.array(values)

def rk4(timestep = None,
        max_time = None,
        initial_time = None,
        initial_val = None,
        deriv = None,
        deriv_params = None):
    '''timestep = time between time values
    max_time = the last time calculated for
    initial_time = the first time calculated for
    initial_val = initial value
    deriv = derivative function
    deriv_params = derivative function parameters'''
    times = np.linspace(initial_time, max_time, timestep)
    values = [initial_val]
    val = initial_val
    for time in times:
        k1 = timestep*deriv(val, **deriv_params)
        k2 = timestep*deriv(val+(1/2)*k1, **deriv_params)
        k3 = timestep*deriv(val+(1/2)*k2, **deriv_params)
        k4 = timestep*deriv(val+k3, **deriv_params)
        val = val+(1/6)*(k1+2*k2+2*k3+k4)
        values.append(val)
    return np.array(values)

if __name__ == "__main__":
    K = 10 #10 billion
    P0 = 1 #billion, population at year 1800
    r = 0.014 #1.4% growth rate
    start_year = 1800
    max_year = 2300
    max_time = max_year - start_year
    years_since_start = np.arange(0, max_time)
    timestep=25

    #calculates population over time period analytically
    analytic_sol = analytic_solution(years_since_start, K=K, P0=P0, r=r)

    #calculates population over time period using Euler's method
    eulersol = forward_euler(initial_val=P0, initial_time=0,
                             timestep=timestep,
                             max_time=max_time, deriv=dP_dt,
                             deriv_params={'K':K, 'r':r})

    #calculates population over time period using rk4 method
    rk4sol = rk4(initial_val=P0, initial_time=0,
                             timestep=timestep,
                             max_time=max_time, deriv=dP_dt,
                             deriv_params={'K':K, 'r':r})
    
    #graphs population growth over time period
    plt.plot(years_since_start+1800, analytic_sol)
    plt.scatter(np.linspace(P0,max_time, timestep)+1800, eulersol[:-1], color='orange')
    plt.scatter(np.linspace(P0,max_time, timestep)+1800, rk4sol[:-1], color='green')
    plt.xlabel('Year')
    plt.ylabel('Population (Billions)')
    plt.title('Population Growth')
    plt.legend(['Analytic Method',"Euler's Method",'RK4 Method'])
    plt.savefig('Population_Growth.png')
    plt.show()
