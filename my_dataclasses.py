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
class UncertaintySource:
    name: str
    is_relative: bool
    is_symmetric: bool
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
        
        #If the uncertainty is relative: apply factor 100 correction to change percentage value to fraction
        if self.is_relative:
            factor *= 100
        
        #Populate the missing sigma or the missing bound
        if self.sigma is None:
            self.sigma = self.bound / factor
        else:
            self.bound = self.sigma * factor
        #Correct for one-sidedness
        if not self.is_symmetric:
            self.sigma = self.sigma / 2
            

@dataclass
class VariableUncertainty: ###!!! Handle some stuff in post-init?
    var_name:                               str
    variable:                               "Variable" #Placeholder!
    direct_uncertainty_sources:             Optional[list] = None                   
    direct_uncertainties:                   Optional[Union[np.ndarray, float]] = None ###!!! EXPLANATIONS
    direct_uncertainty_contributions:       Optional[Union[np.ndarray, float]] = None
    total_direct_uncertainty_contribution:  Optional[Union[np.ndarray, float]] = None
    #total_direct_uncertainty:               Optional[Union[np.ndarray, float]] = None
    
    #dependency_uncertainties:               Optional[dict] = None                      ###!!!
    #dependency_uncertainty_contributions:   Optional[Union[np.ndarray, float]] = None
    
    dependency_uncertainty_names:           Optional[list] = None                       ###!!!
    weighted_dependency_uncertainties:      Optional[np.ndarray] = None
    dependency_uncertainty_contributions:   Optional[Union[np.ndarray, float]] = None
    #total_dependency_uncertainty:           Optional[Union[np.ndarray, float]] = None
    
    
    
    total_uncertainty:                      Optional[Union[np.ndarray, float]] = None
    contribution_split:                     Optional[np.ndarray] = None
    correlation:                            Optional[Union[np.ndarray, float]] = None
    is_calculated:                          Optional[bool] = False
    
    def __post_init__(self):
        self.direct_uncertainty_sources = []
        #self.dependency_uncertainties = {}                                         ###!!!
        
        """ Handle initialization of non-inputs here! """
            
        

        
class Variable:
    def __init__(self, name, description=None, values=None, is_basic=True, aggregation_rule=None, \
                 equation=None, dependency_names=None, first_time=None, last_time=None, \
                     is_timesum=False, timesum_settings=None):
        self.name = name                    #str: variable name
        self.description = description      #str: variable description
        self.is_basic = is_basic            #bool: defines whether variable is basic or derived
        self.aggregation_rule = aggregation_rule      #str: defines the quantity aggregation rule
        
        self.values = values                #[int, float, array, None]: variable values
        self.partial_values = {}            #dict: dictionary of values of the evaluated partial derivatives of this variable with respect to each dependency

        #If it is a derived variable then an equation and constituent variables should be passed
        if not self.is_basic:
            if equation is None:
                raise ValueError(f"{self.name} is defined as a derived variable, please include equation")
        self.equation = equation            #str: defines variable equation
        self.dependency_names = dependency_names      #list: lists variables in equation
        self.dependencies = {}              #dict: dict of variables that are the direct dependencies
        self.is_root_consistent = False     #bool: whether the variable traces consistently to basic variables
        
        self.is_timesum = is_timesum        #bool: defines whether variable is a timesum
        self.timesum_settings = timesum_settings #list of options
        
        self.sympy_symbol_map = None        #dict: dictionary of sympy symbols
        self.sympy_equation = None          #sympy interpretable of the variable equation
        self.executable = None              #executable: sympy-built executable of the variable equation - excluding timesums
        self.partial_executables = None     #dict: dictionary of executables for the partial derivates of the variable for each dependency
        self.calculation_engine = None      #calculation engine: necessary for more complex calculations
                
        self.first_time = None        #datetime object: time of the first datapoint
        self.last_time = None            #datetime object: time of the last datapoint
        if first_time is not None and last_time is not None:
            self.addTimeStep([first_time, last_time])        #float: number of seconds per timestep
        else:
            self.timestep = None                
        
        self.uncertainty = VariableUncertainty(var_name=self.name, variable=self)
        
        #self.direct_uncertainty_sources = []        #list: lists uncertainties of variable
        #self.direct_uncertainties = None            #array: array of direct uncertainty magnitudes, 
        #if not self.is_basic:
        #    self.dependency_uncertainties = {}      #dict: per dependency contains total uncertainty information
        #self.total_uncertainty = None       #[float, array]: total uncertainty of variable (per timestep)
    
    
    def __str__(self):
        if self.is_basic:
            return (f"Basic variable: {self.name} \nDescription: {self.description}")
        else:
            return (f"Derived variable: {self.name} = {self.equation} \nDescription: {self.description}")
    
    def __neg__(self):
        return -self.values
    
    def __add__(self, other):
        #If floats or ints are involved the logic is simple
        if isinstance(other, (int, float)):
            return self.values + other
        if isinstance(self.values, (int, float)) or isinstance(other.values, (int, float)):
            return self.values + other
        #Check array lengths and start/end time comparison
        if len(self.values) == len(other.values):
            if self.timestep != other.timestep:
                print(f"WARNING: summing variables {self.name} and {other.name} with different timesteps, {self.timestep} and {other.timestep}!!")
            return self.values + other.values
        else:
            raise ValueError("Summing variables {self.name} and {other.name} failed. {self.name} has length {len(self.values}, {other.name} has length {len(other.values)}")
    
    def __sub__(self, other):
        #If floats or ints are involved the logic is simple
        if isinstance(other, (int, float)):
            return self.values - other
        if isinstance(self.values, (int, float)) or isinstance(other.values, (int, float)):
            return self.values - other
        #Check array lengths and start/end time comparison
        if len(self.values) == len(other.values):
            if self.timestep != other.timestep:
                print(f"WARNING: subtracting variables {self.name} and {other.name} with different timesteps, {self.timestep} and {other.timestep}!!")
            return self.values - other.values
        else:
            raise ValueError("Subtracting variables {self.name} and {other.name} failed. {self.name} has length {len(self.values}, {other.name} has length {len(other.values)}")
         
    def __mul__(self, other):
        #If floats or ints are involved the logic is simple
        if isinstance(other, (int, float)):
            return self.values * other
        if isinstance(self.values, (int, float)) or isinstance(other.values, (int, float)):
            return self.values * other
        #Check array lengths and start/end time comparison
        if len(self.values) == len(other.values):
            if self.timestep != other.timestep:
                print(f"WARNING: multiplying variables {self.name} and {other.name} with different timesteps, {self.timestep} and {other.timestep}!!")
            return self.values * other.values
        else:
            raise ValueError("Multiplying variables {self.name} and {other.name} failed. {self.name} has length {len(self.values}, {other.name} has length {len(other.values)}")
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __truediv__(self, denominator):
        #If floats or ints are involved the logic is simple
        if isinstance(denominator, (int, float)):
            return self.values / denominator
        if isinstance(self.values, (int, float)) or isinstance(denominator.values, (int, float)):
            return self.values / denominator
        #Check array lengths and start/end time comparison
        if len(self.values) == len(denominator.values):
            if self.timestep != denominator.timestep:
                print(f"WARNING: dividing variables {self.name} and {denominator.name} with different timesteps, {self.timestep} and {denominator.timestep}!!")
            return self.values / denominator.values
        else:
            raise ValueError("Dividing variables {self.name} and {denominator.name} failed. {self.name} has length {len(self.values}, {denominator.name} has length {len(denominator.values)}")
    
    def __rtruediv__(self, numerator):
        #If floats or ints are involved the logic is simple
        if isinstance(numerator, (int, float)):
            return numerator / self.values
        if isinstance(self.values, (int, float)) or isinstance(numerator.values, (int, float)):
            return numerator / self.values
        #Check array lengths and start/end time comparison
        if len(self.values) == len(numerator.values):
            if self.timestep != numerator.timestep:
                print(f"WARNING: dividing variables {self.name} and {numerator.name} with different timesteps, {self.timestep} and {numerator.timestep}!!")
            return numerator / self.values
            raise ValueError("Dividing variables {self.name} and {numerator.name} failed. {self.name} has length {len(self.values}, {numerator.name} has length {len(numerator.values)}")
    
    def __pow__(self, power):
        if not power.is_integer():
            raise ValueError(f"Error, taking power {power} of variable {self.name} is not supported, integer powers only.")
        return self.values**power
        
    def defineValues(self, values):
        self.values = values
    
    def hasValues(self):
        if self.values is None:
            return False
        else:
            return True
        
        
    def hasDirectUncertainties(self):               ###!!! Name
        if len(self.direct_uncertainties) == 0:
            return False
        return True
    
    def addUncertaintySource(self, uncertainty):        ###!!! name
        #Can be a single uncertainty source or a whole list of them
        if type(uncertainty) is list:
            #self.direct_uncertainties.extend(uncertainty)
            self.uncertainty.direct_uncertainty_sources.extend(uncertainty)
        else:
            #self.direct_uncertainties.append(uncertainty)
            self.uncertainty.direct_uncertainty_sources.append(uncertainty)
            
    
    def addTimeStep(self, time_range):
        """ Given the time range (which are the times of first and last datapoint registrations) - calculates timestep length """
        self.first_time = time_range[0]
        self.last_time = time_range[1]
        
        if self.values is None or len(self.values)<2:
            raise ValueError(f"Automatic timestep calculation for variable {self.name} failed: no or too little values are specified (at least 2).")
       
        #Calculate timestep
        timestep = (self.last_time - self.first_time) / (len(self.values) - 1)
        
        self.start_time = self.first_time - timestep/2
        self.end_time = self.last_time + timestep/2
        
        self.timestep = timestep.total_seconds()
        if not self.timestep.is_integer():
            raise ValueError(f"Timestep {self.timestep} is not integer, please ensure timesteps are an integer amount of seconds!")
        self.timestep = int(self.timestep)
    
    
    def executeEquation(self, store_results=True, force_recalculation=False, calculation_engine=None):
        """ Executes equation built by sympy and checks whether everything is initialized correctly """
        if self.equation is None:
            raise ValueError(f"Tried to execute the equation of variable {self.name}, for which no equation is defined.")
        elif self.executable is None:
            raise ValueError(f"Tried to execute the equation of variable {self.name} = {self.equation}, but no equation executable has been built for this variable.")
        
        #If values already defined: do nothing unless forced recalculation is desired
        if self.values is not None:
            if force_recalculation is True:
                print(f"WARNING: executing equation of variable {self.name} while values are already defined!")
            else:
                return self.values
        
        #Check if this variable is a timesum, in which case the executable is the equation INSIDE the timesum
        if self.is_timesum:
            if calculation_engine is None:
                raise ValueError(f"Asked to perform timesum of variable {self.name} while no equation engine is given to this variable")
            calculated_values = calculation_engine.timeSum(self)
        else:
            args = [self.dependencies[dep_name] for dep_name in self.dependency_names]
            calculated_values = self.executable(*args)
        
        #Optionally store the result as the new variable values
        if store_results:
            self.values = calculated_values
        return calculated_values


    def executePartialDerivative(self, dep_name, absolute_values=False, store_results=True, force_recalculation=False):
        """ Executes the partial derivative executable of the variable for a given dependency, optionally stores values in partial_values dictionary """
        #If no forced recalculation and if the values are already calculated we simply return the already calculated values
        if force_recalculation is False and dep_name in self.partial_values:
            return self.partial_values[dep_name]
        #If there is no executable for this dependency we raise an error
        if self.partial_executables is None:
            raise ValueError(f"Tried to evaluate partial derivative of variable {self.name} while partial derivative executables have not been built yet.")
        
        #Get partial executable, pass arguments
        partial_executable = self.partial_executables[dep_name]
        args = [self.dependencies[dep_name] for dep_name in self.dependency_names]
        calculated_values = partial_executable(*args)
        
        #In case of a trivial equation, calculated values will be a Variable object. Here we fix that
        if isinstance(calculated_values, Variable):     ###!!! change this to be handled through an equation engine wrapper
            calculated_values = calculated_values.values
        
        if absolute_values:
            calculated_values = np.abs(calculated_values)
        
        #Optionally store results
        if store_results:
            self.partial_values[dep_name] = calculated_values
        return calculated_values
    
    def executeAllPartials(self, absolute_values=False, store_results=True, force_recalculation=False):
        """ Evaluates all partial derivatives of a variable """
        partial_values = {}
        for dep_name in self.dependency_names:
            partial_values[dep_name] = self.executePartialDerivative(dep_name, absolute_values=absolute_values, store_results=store_results, force_recalculation=force_recalculation)
        return partial_values
        
    








            
            