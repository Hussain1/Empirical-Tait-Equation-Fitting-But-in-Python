import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
% This code takes an input of known pressure and volume points and
% uses a mix of linear algebra and brute force to fit parameters
% Use it like TaitEqV('C:\blah\blah\points.xlsx',L,H)
% 'C:\blah\blah\points.xlsx' is path to the xlsx sheet sandwiched in ''
% L and H are two bounds where you expect your Tait Equation C parameter to be
% It will output three Tait Equation Paramters Fit C, B, and A for the form
% V = A + B*ln(1+P/C)
% V, A and B are in cm3/mol. P is in bar.
% The Deciding Factor for the optimal C is the RMSE from the known points
% Make sure that the excel sheet has the FIRST COLUMN as PRESSURE and  
% the SECOND as VOLUME with no additional things anywhere
"""

def TaitFit(Path, L, H): #the equation takes a path to an .xlsx sheet, an expected lower bound and higher bound for C
    FitValues = pd.read_excel(Path, header=None, index_col=None) #read the Excel file and assigned to var FaitValues
    Pset = np.array(FitValues[:][0].copy()) #extract the pressure column
    Vset = np.array(FitValues[:][1].copy()) #extract the volume column
    n = len(Pset) #count how many points are there, which will be useful for the linear regression calculations
    Values_Stack = []  #define a couple of empty lists which we will use to stack values and find the minimum later
    Values_Stack2 = []
    
    for i in range(21): #lets do 20 runs of optimization
        if i == 0: #first run (genesis matrix) will run the fit only based on the bounds inputted by user (1000 steps)
            for C in np.linspace(L, H, 1000):    #here we brute force the C
                xset = np.log(1 + (Pset/C))     #linear regression for this and the next two lines
                B = (np.mean(xset*Vset)-np.mean(Vset)*np.mean(xset))/(np.mean(xset**2)-(np.mean(xset))**2)
                A = np.mean(Vset)-B*np.mean(xset)
                ErrorSet =(Vset-A-B*xset)**2    #calculate the RMSD
                if len(Values_Stack) == 0:   #stack every run into a the matrix Values_Stack
                    Values_Stack = np.array([C,B,A,((sum(ErrorSet)/n)**0.5)])
                else:
                    Values_Stack = np.vstack((Values_Stack,np.array([C,B,A,((sum(ErrorSet)/n)**0.5)])))
        else:
            minErrorIndex = np.argmin(Values_Stack,axis=0)[3]  #extract where the minimum RMSE index is and use it next to optimize on C
            for C in np.linspace(Values_Stack[minErrorIndex,0]*(1-0.01/i), Values_Stack[minErrorIndex,0]*(1+0.01/i), 1000):  #define a shrinking lower bound, shrinking upper bound and shrinking step size  and run the same thing as before
                xset = np.log(1 + (Pset/C))
                B = (np.mean(xset*Vset)-np.mean(Vset)*np.mean(xset))/(np.mean(xset**2)-(np.mean(xset))**2)
                A = np.mean(Vset)-B*np.mean(xset)
                ErrorSet =(Vset-A-B*xset)**2
                if len(Values_Stack2) == 0:  #stack all the optimization runs in the matrix Values_Stack2
                    Values_Stack2 = np.array([C,B,A,((sum(ErrorSet)/n)**0.5)])
                else:
                    Values_Stack2 = np.vstack((Values_Stack2,np.array([C,B,A,((sum(ErrorSet)/n)**0.5)])))
                    
    minErrorIndex2 = np.argmin(Values_Stack2,axis=0)[3]  #extract the index of the minimum RMSE and use it next to extract the what A, B, C and the RMSE are at that row
    [C, B, A, Error] = Values_Stack2[minErrorIndex2].round(8) #Edit how many significant figures you want here
    
    #plotting the results
    x = np.arange(Pset.min(), Pset.max(), 50) #based on the simulation points, define pressure points to use with the fit
    fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(15,6))
    ax[0].plot(x, A+B*np.log(1+x/C), color='blue', label='Tait Equation Correlation') #the pressure points are the x, and they are used in the Tait Equation for y
    ax[0].scatter(Pset, Vset, color='red', label='Simulation Points')  #overlay the simulation points to visually represent how if the the fit went well
    ax[0].set_xlabel('Pressure')
    ax[0].set_ylabel('Volume')
    ax[0].set_title('Tait Equation Correlation With Simulation Poitns Overlayed')
    ax[1].plot(Values_Stack[:,0],Values_Stack[:,3],color='orange') #showing how the value of C correlates with the RMSE
    ax[1].set_xlabel('C')
    ax[1].set_ylabel('RMSE')
    ax[1].set_title('C vs RMSE')
    plt.show()
    
    return print('Your C is {}, your B is {}, your A is {}, and your RMSE is {}'.format(C, B, A, Error))