# -*- coding: utf-8 -*-
"""
by Kelly McQuighan 2017

These tools can be used to visualize different numerical integration schemes,
and to compute the associated error. They can also be used to find the order of
various numerical schemes.

"""
from matplotlib import pyplot as plt
from numpy import *
import numpy as np
import scipy.integrate as inte
import scipy.interpolate as interp
import matplotlib as mpl

mpl.rcParams['font.size'] = 17
colors = ['#0058AF','#FF8000','#D682FF','#00693C','#E02102']
styles = ['-',':','-',':','-']

"""
This function is used to make the plot of what the approximation of the integral
looks like for five different numerical methods: Left Riemann sum, Right Riemann sum,
Midpoint rule, Trapezoid rule, and Simpson's rule.
"""
def plots(func, a,b,n,method,ax):
    
    ax.axvline(0.,color='#666666',linewidth=1)
    ax.axhline(0.,color='#666666',linewidth=1)
    if (a>0):
        xlarge = np.linspace(0.,1.1*b,1000)
    elif (b<0):
        xlarge = np.linspace(1.1*a,0.,1000)
    else:
        xlarge = np.linspace(1.1*a,1.1*b,1000)
    flarge = func(xlarge)
    ax.plot(xlarge,flarge,'b', linewidth=5)
    
    ax.set_xlim([xlarge[0], xlarge[999]])
    smallest = np.min(flarge)
    largest = np.max(flarge)
    
    dx = (b-a)/n
    
    xs = np.linspace(a,b,n+1)
    fxs = func(xs)
    
    if method.lower()=='left':
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fxs[i]], [xs[i+1], fxs[i]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)
    elif method.lower()=='right':
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fxs[i+1]], [xs[i+1], fxs[i+1]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)
    elif method.lower()=='midpoint':
        x = np.linspace(a+dx/2.,b-dx/2.,n)
        fx = func(x)
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fx[i]], [xs[i+1], fx[i]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)
    elif method.lower()=='trapezoid':
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fxs[i]], [xs[i+1], fxs[i+1]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)  
    elif method.lower()=='simpson':
        # note: this implementation keeps the number of grid points the same         
        for i in range(0,n,2):
            lag = interp.lagrange([xs[i], xs[i+1], xs[i+2]], [fxs[i], fxs[i+1], fxs[i+2]])
            section = np.linspace(xs[i], xs[i+2], 100)
            fsec = lag(section)
            ax.fill_between(section,fsec, facecolor='g', edgecolor='g', alpha=0.3, linewidth=3)
        
        x_mid = np.linspace(a+dx,b-dx,n/2)
        vert = np.ones(100)   
        for the_x in x_mid:
            ax.plot(the_x*vert,np.linspace(0,func(the_x),100),'g--', linewidth=3, alpha=0.5)
        
    else:
        print ('ERROR: You have not specified a valid method. Please check for typos.')
    
    if smallest>0:
        ax.set_ylim([0,1.1*largest])
    elif largest<0:
        ax.set_ylim([1.1*smallest,0])
    else:
        ax.set_ylim([1.1*smallest, 1.1*largest])
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')

"""
This function is used to make the plot of what the approximation of the integral
looks like for all five different numerical methods: Left Riemann sum, Right Riemann sum,
Midpoint rule, Trapezoid rule, and Simpson's rule.

It also makes a bar chart showing the size of the error for each method so that
the user can quickly see which method is better for a specific fixed value of n.
"""
def plotArea(f,a,b,n):

    a = eval(a)
    b = eval(b)
    
    if n<1:
        n=1
        print ('ERROR: n must be greater than zero. setting n=1.')
    
    func = eval("lambda x: " + f)
    I = inte.quad(func, a, b)[0]
    
    fig = plt.figure(figsize=(15, 6))
    
    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4)
    ax5 = fig.add_subplot(2,3,5)
    ax6 = fig.add_subplot(2,3,6)
    
    plots(func,a,b,n,"left",ax1)
    plots(func,a,b,n,"right",ax2)
    plots(func,a,b,n,"midpoint",ax3)
    plots(func,a,b,n,"trapezoid",ax4)
    plots(func,a,b,n,"simpson",ax5)
    
    area1 = evalArea(func,a,b,n,"left")
    area2 = evalArea(func,a,b,n,"right")
    area3 = evalArea(func,a,b,n,"midpoint")
    area4 = evalArea(func,a,b,n,"trapezoid")
    area5 = evalArea(func,a,b,n,"simpson")
        
    err1 = np.abs(area1-I)
    err2 = np.abs(area2-I)
    err3 = np.abs(area3-I)
    err4 = np.abs(area4-I)
    err5 = np.abs(area5-I)
    
    if (not check_error(err1)): err1=0;
    if (not check_error(err2)): err2=0;
    if (not check_error(err3)): err3=0;
    if (not check_error(err4)): err4=0;
    if (not check_error(err5)): err5=0;
       
    ax6.bar(range(5),[err1,err2,err3,err4,err5])
    ax6.set_xticks(range(5))
    ax6.set_xticklabels(['left','right','mid','trap','Simp'],rotation=70)
    ax6.axhline(0,color='k',linewidth=1)
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=5.5)
    plt.suptitle('f(x) = '+f+', Area = %.3f, n=%d' %(I,n), fontsize=20, y=1.2)
    ax1.set_title('Method "Left"\n Approximate area:%.5f \n Absolute error: %.2e' %(area1, err1))
    ax2.set_title('Method "Right"\n Approximate area:%.5f \n Absolute error: %.2e' %(area2, err2))
    ax3.set_title('Method "Midpoint"\n Approximate area:%.5f \n Absolute error: %.2e' %(area3, err3))
    ax4.set_title('Method "Trapezoid"\n Approximate area:%.5f \n Absolute error: %.2e' %(area4, err4))
    ax5.set_title('Method "Simpson"\n Approximate area:%.5f \n Absolute error: %.2e' %(area5, err5))
    ax6.set_title('Absolute error for each method\n')

    plt.show()
    
"""
This method plots the approximation three times, each time doubling the number of
gridpoints used in the approximation. It also computes and outputs how the error
decreases.
"""
def plot3Areas(f,a,b,n,method):

    a = eval(a)
    b=eval(b)
        
    func = eval("lambda x: " + f)
    I = inte.quad(func, a, b)[0]
    
    plt.figure(figsize=(15, 4))
    
    ax1 = plt.subplot2grid((4,3), (0, 0),rowspan=3)
    ax2 = plt.subplot2grid((4,3), (0, 1),rowspan=3)
    ax3 = plt.subplot2grid((4,3), (0, 2),rowspan=3)
    ax0 = plt.subplot2grid((4,3),(3,0),colspan=3)
    ax0.axis('off')

    plots(func,a,b,n,method,ax1)
    plots(func,a,b,2*n,method,ax2)
    plots(func,a,b,4*n,method,ax3)
    
    area1 = evalArea(func,a,b,n,method)
    area2 = evalArea(func,a,b,2*n,method)
    area3 = evalArea(func,a,b,4*n,method)
    err1 = np.abs(area1-I)
    err2 = np.abs(area2-I)
    err3 = np.abs(area3-I)
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.suptitle('f(x) = '+f+', Method: '+method+', Area = %.3f' %I, fontsize=20, y=1.4)
    ax1.set_title('n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(n,area1, err1))
    ax2.set_title('n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(2*n,area2, err2))
    ax3.set_title('n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(4*n,area3, err3))
    
    if (not check_error(err1)):
        ax0.text('Using method '+method+' to find the area under f(x) = '+f+' returns no errors, so it does not make sense to compare the errors for different numbers of sub-intervals.',
                 ha='left', va='top', fontsize=20, transform=ax0.transAxes)    
    else:
        ax0.text(0.0, 1., 'When using method '+method+' to compute the area under f(x) = '+f+':\n'+
                 '- In doubling the number of subintervals from n=%d to n=%d the error was decreased by a factor of %.2f\n' %(n, 2*n,err1/err2)+
                '- In doubling the number of subintervals from n=%d to n=%d the error was decreased by a factor of %.2f' %(2*n, 4*n,err2/err3),                                                                                                         
                 ha='left', va='top', fontsize=20, transform=ax0.transAxes) 
        
    plt.show()

"""
This method is currently unused by the Notebook because it refreshes in an awkward way.
Instead plot one method at a time using plot3Areas and a Dropdown box for the method type. 
"""       
def plotAllMethods(f,a,b,n):

    a = eval(a)
    b=float(b)
    
    if n<1:
        n=1
        print ('ERROR: n must be greater than zero. setting n=1.')
    
    plot3Areas(f,a,b,n,"left")
    plot3Areas(f,a,b,n,"right")
    plot3Areas(f,a,b,n,"midpoint")
    plot3Areas(f,a,b,n,"trapezoid")
    plot3Areas(f,a,b,n,"simpson")
    plt.show()

"""
This function approximates the integral using one of five possible numerical methods:
Left Riemann sum, Right Riemann sum, Midpoint Rule, Trapezoid Rule, and Simpson's Rule.
"""     
def evalArea(func,a,b,n, method):
    
    dx = (b-a)/n
    
    if method.lower()=='left':
        x = np.linspace(a,b-dx,n)
        fx = func(x)
        area = np.sum(fx)*dx
    elif method.lower()=='right':
        x = np.linspace(a+dx,b,n)
        fx = func(x)
        area = np.sum(fx)*dx
    elif method.lower()=='midpoint':
        x = np.linspace(a+dx/2.,b-dx/2.,n)
        fx = func(x)
        area = np.sum(fx)*dx
    elif method.lower()=='trapezoid':
        x = np.linspace(a+dx,b-dx,n-1)
        fx = func(x)
        area = dx*(0.5*func(a)+0.5*func(b)+np.sum(fx))
    elif method.lower()=='simpson':
        x_mid = np.linspace(a+dx,b-dx,n/2)
        x_trap = np.linspace(a+2*dx,b-2*dx,n/2-1)
        fx_mid = func(x_mid)
        fx_trap = func(x_trap)
        area = dx/3.*(4*np.sum(fx_mid)+func(a)+func(b)+2*np.sum(fx_trap))
    else:
        print ('ERROR: You have not specified a valid method. Please check for typos.')
    
    return area

"""
Checks if the error is near machine precision. If so it does not make sense to 
compare how the error decreases as teh gridsizes increase. For example, evaluating
a constant function using any of the methods will be exact, so the error should be machine
precision.
"""
def check_error(err):
    epsilon = 7./3 - 4./3 -1
    return (err>100*epsilon)
    
"""
This function compares all five methods by making a log-log plot of the error. 

The slope of each curve is computed and can be used to determine the order of each
numerical method.
"""
def compareMethods(f,a,b):
    
    n = 8
    a = eval(a)
    b=float(b)
    
    if n<1:
        n=1
        print ('ERROR: n must be greater than zero. setting n=1.')
        
    func = eval("lambda x: " + f)
    I = inte.quad(func, a, b)[0]
    
    n = int(n)
    if n<1:
        n=1
        
    errors = np.ones((5,6)) # methods in rows, errors in columns   
    ns = np.zeros((6))
    for i in range(6):
        errors[0,i] = np.abs(I-evalArea(func, a,b,2**i*n, 'left'))
        errors[1,i] = np.abs(I-evalArea(func, a,b,2**i*n, 'right'))
        errors[2,i] = np.abs(I-evalArea(func, a,b,2**i*n, 'midpoint'))
        errors[3,i] = np.abs(I-evalArea(func, a,b,2**i*n, 'trapezoid'))
        errors[4,i] = np.abs(I-evalArea(func, a,b,2**i*n, 'simpson'))
        
        ns[i] = i
        
    fig = plt.figure(figsize=(10,5))
    ax = fig.gca()    
        
    ax.set_xlabel(r'$\log_2(n/4)$')
    ax.set_ylabel(r'$\log_2(error_i/error_1)$')
    ax.set_title('log-log plot of the error versus the number of subintervals.\n f(x) = ' + f)
    
    labels = ['Left Riemann', 'Right Riemann', 'Midpoint', 'Trapezoid','Simpson']
    
    for i in range(5):
        if check_error(errors[i,0]):
            log_errors  = np.log(errors[i,:] / errors[i,0])/np.log(2)
            plt.plot(ns, log_errors,linewidth=6,color=colors[i],linestyle=styles[i],label=labels[i])
            poly1 = np.polyfit(ns, log_errors, 1)
            print ('Using method '+labels[i]+' the slope of the log-log plot is %.2f' %poly1[0])
        else:
            print (r'Using method '+labels[i]+' the errors are less than machine precision for n=8 already!')
        
    plt.legend(loc=3)
    
    plt.axhline(0.0,0,5,color='k',linewidth=1)
    ax.set_xlim([0,5])
    
    plt.show()
