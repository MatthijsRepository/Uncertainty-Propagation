from dataclasses import dataclass
from typing import Union, Optional
import numpy as np
from datetime import datetime


@dataclass 
class CSVData:
    name: str
    data: dict
    timestep: Optional[float] = None
    time_range: Optional[list] = None


@dataclass
class Uncertainty:
    name: str
    is_relative: bool
    is_onesided: bool
    shape: str
    correlation: Union[float, int, np.ndarray]
    #Mutually exclusive fields
    sigma: Optional[Union[float, int, np.ndarray]] = None
    bound: Optional[Union[float, int, np.ndarray]] = None
    
    def __post_init__(self):
        #Check input consistency
        if not 0 <= self.correlation <= 1:
            raise ValueError(f"Uncertianty correlation should be between 0 and 1, is {self.correlation} for uncertainty {self.name}.")
        if not (self.shape=="rectangular" or self.shape=="normal" or self.shape=="triangular" or self.shape=="U-shaped"):
            raise ValueError(f"Uncertainty distribution shape: '{self.shape}' not recognized.")
        if (self.sigma is None and self.bound is None):
            raise ValueError(f"Either a standard deviation or a bound should be provided for uncertainty {self.name}.")
        
        #Populate sigma or bound
        #First, determine difference factor which is dependent on the distribution shape
        if self.shape.lower()=="normal":
            factor = 3
        elif self.shape.lower()=="rectangular":
            factor = np.sqrt(3)
        elif self.shape.lower()=="triangular":
            factor = np.sqrt(6)
        elif self.shape.lower()=="u-shaped":
            factor = np.sqrt(2)
        #Populate the missing sigma or the missing bound
        if self.sigma is None:
            self.sigma = self.bound / factor
        else:
            self.bound = self.sigma * factor
        #Correct for one-sidedness
        if self.is_onesided:
            self.sigma = self.sigma / 2
        
        
class Variable:
    def __init__(self, name, description, values=None, is_basic=True, is_timesum=False,\
                 equation=None, dependency_names=None, start_time=None, end_time=None):
        self.name = name                    #str: variable name
        self.description = description      #str: variable description
        self.is_basic = is_basic            #bool: defines whether variable is basic or derived
        self.is_timesum = is_timesum        #bool: defines whether variable is a timesum
        #If it is a derived variable then an equation and constituent variables should be passed
        if not self.is_basic:
            if not (equation or dependency_names):
                raise ValueError(f"{self.name} is defined as a derived variable, please include equation and variables")
            self.equation = equation        #str: defines variable equation
            self.dependency_names = dependency_names      #list: lists variables in equation
            self.is_root_consistent = None     #bool: whether the variable is traces consistently to basic variables
            
        self.values = values                #[int, float, array]: variable values
        
        self.start_time = start_time        #datetime object: start time of variable
        self.end_time = end_time            #datetime object: end time of variable
        self.timestep = None                #float: number of seconds per timestep
        
        self.uncertainties = []             #list: lists uncertainties of variable
        self.total_uncertainty = None       #[float, array]: total uncertainty of variable (per timestep)
    
    
    def __str__(self):
        if self.is_basic:
            return (f"Basic variable: {self.name} \nDescription: {self.description}")
        else:
            return (f"Derived variable: {self.name} = {self.equation} \nDescription: {self.description}")
        
    def __add__(self, other):
        return self.values + other.values

    def DefineValues(self, values):
        self.values = values
    
    def AddTimestep(self, step, time_range=None):
        #Check if start and end times are provided
        if not time_range is None:
            self.start_time = time_range[0]
            self.end_time = time_range[1]
        #Direct timestep definition
        if (type(step)==int or type(step)==float):
            self.timestep = step
        #Automatic timestep definition
        elif step=="auto":
            #Check if start and end times are defined
            if (self.start_time==None or self.end_time==None):
               raise ValueError(f"Automatic timestep calculation for variable {self.name} failed: requires specified start and end timestamps.") 
            #Check if the variable has enough values
            if (self.values is None or len(self.values)<2):
                raise ValueError(f"Automatic timestep calculation for variable {self.name} failed: no or too little values are specified (at least 2).")
            #Calculate timestep
            #fmt = "%H:%M:%S"
            #start_time = datetime.strptime(self.start_time, fmt)
            #end_time = datetime.strptime(self.end_time, fmt)
            timestep = (self.end_time - self.start_time) / (len(self.values) - 1)
            self.timestep = timestep.total_seconds()
        else:
            raise ValueError(f"Timestep {str(step)} of type {str(type(step))} is not supported.")
    
    def HasUncertainty(self):
        if len(self.uncertainties) == 0:
            return False
        return True
    
    def AddUncertaintySource(self, uncertainty):
        #Can be a single uncertainty source or a whole list of them
        if type(uncertainty) is list:
            self.uncertainties.extend(uncertainty)
        else:
            self.uncertainties.append(uncertainty)
    
    def TimeSum(self, interval=None, start_time=None, end_time=None):
        #check if time-dependent
        if (self.values is None or len(self.values)<2):
            raise ValueError(f"Time series sum for variable {self.name} can not be done: no or too little values.")
        #check if there is a timestep, and try to automatically calculate one
        if self.timestep is None:
            try:
                self.AddTimeStep("auto")
            except:
                raise ValueError(f"Automatic timestep calculation for timesum of variable {self.name} failed, please provide a timestep or time bounds.")
        #Handle partial timesum calculations
        #
        #if not (interval is None):
            #Check if start and end time are provided
        
        #Calculate total timesum
        #self.timesum = np.sum(self.values) * self.timestep
        #return self.timesum
        return Variable(name=f"Daily timesum of {self.name}", description=f"Daily timesum - {self.description}", values=np.sum(self.values) * self.timestep, is_basic=False, is_timesum=True, equation=f"Timesum('{self.name}')", variables=self)
        

"""
if __name__=="__main__":
    T = Variable("T", "Temperature", is_basic=False, equation="S*A", variables="S, A")
    G = Variable("G", "Global Horizontal Irradiance", values=np.array([1,2,3]), is_basic=True, start_time="00:01:20", end_time="00:02:20")
    
    #Timestep test
    G.AddTimestep("auto")
    print(G.timestep)
    #Uncertainty object test
    
    print(T)
    print(G.TimeSum())
    G.AddUncertaintySource(["a", "cookie"])
    G.AddUncertaintySource(["b", "c"])
    print(G.uncertainties)  
    E = G.TimeSum()
    print(E)
    print(E.values)


#"""








            
            