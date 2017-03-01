# -*- coding: utf-8 -*-
"""
Copyright Kelly McQuighan 2016

These tools can be used to visualize different numerical integration schemes,
and to compute the associated error. They can also be used to find the order of
various schemes.

"""
from matplotlib import pyplot as plt
from numpy import *
import numpy as np
import scipy.integrate as inte
import scipy.interpolate as interp
import matplotlib as mpl

mpl.rcParams['font.size'] = 17
colors = ['#0058AF','#FF8000','#D682FF','#00693C','#E02102']
#mpl.rcParams['axes.color_cycle'] = colors
styles = ['-',':','-',':','-']
#mpl.rcParams['axes.prop_cycle'] = cycler('color',colors)+cycler('ls',styles)
####### DO NOT CHANGE ANYTHING BELOW THIS LINE ##########
def plots(func, a,b,n,method,ax):
    
    xlarge = np.linspace(a,b,1000)
    flarge = func(xlarge)
    ax.plot(xlarge,flarge,'b', linewidth=5)
    
    ax.set_xlim([a, b])
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
        print 'ERROR: You have not specified a valid method. Please check for typos.'
    
    if smallest>0:
        ax.set_ylim([0,largest])
    elif largest<0:
        ax.set_ylim([smallest,0])
    else:
        ax.set_ylim([smallest, largest])
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')


def plotArea(f,a,b,n):
    
    n = 2*int(float(n)/2)
    if a[:3]=='np.':
        a = eval(a)
    else:
        a = float(a)
    if b[:3]=='np.':
        b = eval(b)
    else:
        b=float(b)
    
    if n<1:
        n=1
        print 'ERROR: n must be greater than zero. setting n=1.'
    
    func = eval("lambda x: " + f)
    I = inte.quad(func, a, b)[0]
    
    fig = plt.figure(figsize=(15, 6))
    
    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4)
    ax5 = fig.add_subplot(2,3,5)
    
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
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=5.5)
    plt.suptitle('f(x) = '+f+', Area = %.3f' %I, fontsize=20, y=1.2)
    ax1.set_title('With method "Left" and n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(n,area1, err1))
    ax2.set_title('With method "Right" and n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(n,area2, err2))
    ax3.set_title('With method "Midpoint" and n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(n,area3, err3))
    ax4.set_title('With method "Trapezoid" and n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(n,area4, err4))
    ax5.set_title('With method "Simpson" and n=%d \n Approximate area:%.5f \n Absolute error: %.2e' %(n,area5, err5))

    plt.show()
    

def plot3Areas(f,a,b,n,method,I):

    func = eval("lambda x: " + f)
    fig = plt.figure(figsize=(15, 3))
       
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)

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
    plt.show()
    
    if err2==0:
        print 'Using method '+method+' to find the area under f(x) = '+f+' returns no errors, so it does not make sense to compare the errors for different numbers of sub-intervals.'

    else:
        print 'When using method '+method+' to compute the area under f(x) = '+f+':'
        print '- In doubling the number of subintervals from n=%d to n=%d the error was decreased by a factor of %.2f' %(n, 2*n,err1/err2)
        print '- In doubling the number of subintervals from n=%d to n=%d the error was decreased by a factor of %.2f' %(2*n, 4*n,err2/err3)
        
def plotAllMethods(f,a,b,n):
    n = 2*int(float(n)/2)
    if a[:3]=='np.':
        a = eval(a)
    else:
        a = float(a)
    if b[:3]=='np.':
        b = eval(b)
    else:
        b=float(b)
    
    if n<1:
        n=1
        print 'ERROR: n must be greater than zero. setting n=1.'
            
    func = eval("lambda x: " + f)
    I = inte.quad(func, a, b)[0]
    
    plot3Areas(f,a,b,n,"left",I)
    plot3Areas(f,a,b,n,"right",I)
    plot3Areas(f,a,b,n,"midpoint",I)
    plot3Areas(f,a,b,n,"trapezoid",I)
    plot3Areas(f,a,b,n,"simpson",I)
    
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
        print 'ERROR: You have not specified a valid method. Please check for typos.'
    
    return area

def check_error(err):
    epsilon = 7./3 - 4./3 -1
    return (err>100*epsilon)
    
def compareMethods(f,a,b):
    
    n = 8
    if a[:3]=='np.':
        a = eval(a)
    else:
        a = float(a)
    if b[:3]=='np.':
        b = eval(b)
    else:
        b=float(b)
    
    if n<1:
        n=1
        print 'ERROR: n must be greater than zero. setting n=1.'
        
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
            print 'Using method '+labels[i]+' the slope of the log-log plot is %.2f' %poly1[0]
        else:
            print r'Using method '+labels[i]+' the errors are less than machine precision for n=8 already!'
        
    plt.legend(loc=3)
    
    plt.axhline(0.0,0,5,color='k',linewidth=1)
    ax.set_xlim([0,5])
    
    plt.show()