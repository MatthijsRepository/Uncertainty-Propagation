from dataclasses import dataclass
from typing import Union, Optional
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


@dataclass 
class CSVData:
    data: dict
    timestep: Optional[float] = None
    time_range: Optional[list] = None

    
@dataclass
class ParsedVariableData:
    def __post_init__(self):
        self.csv_pointer = None         #str: tells job handler which CSV column this variable should take its values from
        self.is_hardcoded = False       #bool: tells job handler whether variable value is hardcoded. This means this variable is skipped during resets
        self.var_values = None          #[int, float, array]: variable values
        self.timestep = None            #[float, str]: timestep of variable, or setting ###!!!
        self.description = None         #str: variable description
        self.is_basic_variable = None   #bool: flags if variable is basic or derived
        self.is_rate = None             #bool: flags whether variable is a rate (quantity over time) or a quantity
        self.aggregation_rule = None    #str: aggregation rule of the quantity
        self.equation = None            #str: equation of the variable
        self.uncertainties = []         #list: contains all uncertainty sources
        

@dataclass
class TimeHarmonizationData:
    dep_var_name:       str         #Name of the affected variable
    base_timestep:      float       #Original timestep of the variable
    new_timestep:       float       #New timestep of the variable
    new_start_time:     datetime    #New start time
    new_last_time:      datetime    #New end time
    low_index:          int         #Index to start rebinning in original array
    high_index:         int         #Index to end rebinning in original array
    low_fraction:       int         #If a new bin stretches from bin i to bin i+n, fraction of original bin i to include in new bin
    high_fraction:      int         #Identical as above, but for bin i+n. high_fraction = 1-low_fraction: we assume integer upsample factors
    upsample_factor:    int         #Factor increase of new bin size compared to old bin size
    target_var_name:    Optional[str]   = None         #Name of the derived variable that required harmonization of dependencies
    prune_offset_start: Optional[int]   = None         #Number of new timesteps that are discarded at start to ensure a common start time
    prune_offset_end:   Optional[int]   = None         #Number of new timesteps that are discarded at end to ensure a common end time
    #smuggled_time:      Optional[float] = None        #Amount of seconds that  were smuggled during rebinning
    new_values:         Optional[Union[np.ndarray, float]] = None   #Can be used to store new values
    aggregated_correlation: Optional[np.ndarray] = None             #Can be used to cache aggregated correlations
        
        
    def getFirstTime(self):
        return self.new_start_time + timedelta(seconds=self.new_timestep)
    
    def getTotalOffsetSteps(self):
        if self.prune_offset_start is None:
            low_index_total = self.low_index
        else:
            low_index_total = self.low_index + self.upsample_factor * self.prune_offset_start

        if self.prune_offset_end is None:
            high_index_total = self.high_index
        else:
            high_index_total = self.high_index - self.upsample_factor * self.prune_offset_end
        return low_index_total, high_index_total

@dataclass
class UncertaintySource:
    name:               str
    is_relative:        bool
    is_symmetric:       bool
    shape:              str
    correlation:        Union[float, int, np.ndarray]
    multiplier:         Optional[str] = None
    values:             Optional[Union[float, int, np.ndarray]] = None
    aggregated_values:  Optional[dict] = None
    parent_variable:    Optional[str] = None
    #Mutually exclusive fields
    sigma:              Optional[Union[float, int, np.ndarray]] = None
    bound:              Optional[Union[float, int, np.ndarray]] = None

    def __post_init__(self):
        #Check input consistency
        if not 0 <= self.correlation <= 1:
            raise ValueError(f"Uncertianty correlation should be between 0 and 1, is {self.correlation} for uncertainty {self.name}.")
        if not (self.shape=="rectangular" or self.shape=="normal" or self.shape=="triangular" or self.shape=="U-shaped"):
            raise ValueError(f"Uncertainty distribution shape: '{self.shape}' not recognized.")
        if (self.sigma is None and self.bound is None):
            raise ValueError(f"Either a standard deviation or a bound should be provided for uncertainty {self.name}, neither is provided.")
        if (self.sigma is not None and self.bound is not None):
            raise ValueError(f"Either a standard deviation or a bound should be provided for uncertainty {self.name}, not both!")
        
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
        
        #Check if correlation is either 0 or 1 - the only currently supported values
        if self.correlation not in (0,1):
            raise ValueError(f"Correlation for uncertainty source {self.name} is given as {self.correlation}. Currently only correlations 0 or 1 are supported.")
               
        #Populate the missing sigma if a bound is given
        if self.sigma is None:
            self.sigma = self.bound / factor
            self.bound = None
        
        #If the uncertainty is relative: apply factor 100 correction to change percentage value to fraction
        if self.is_relative:
            self.sigma /= 100

        #Correct for one-sidedness
        if not self.is_symmetric:
            self.sigma = self.sigma / 2
        
        #If a multiplier is present: initialize an executable field
        if self.multiplier is not None:
            self.executable = None
        
    
    def getCorrelationMatrix(self, size):
        M = np.ones((size, size)) * self.correlation
        np.fill_diagonal(M, 1)
        return M
            

@dataclass
class VariableUncertainty: ###!!! Handle some stuff in post-init?
    var_name:                               str
    variable:                               "Variable" #Placeholder!
    
    direct_uncertainty_sources:             Optional[list] = None
    
    root_sources:                           Optional[list] = None
    root_weighted_uncertainties:            Optional[dict] = None
    root_upsample_factors:                  Optional[list] = None
    
    aggregated_weighted_uncertainties:      Optional[np.ndarray] = None

    total_uncertainty:                      Optional[Union[np.ndarray, float]] = None
    correlation:                            Optional[Union[np.ndarray, float]] = None
    
    direct_uncertainties_calculated:        Optional[bool] = False
    total_uncertainty_calculated:           Optional[bool] = False
    is_certain:                             Optional[bool] = None
    
    def __post_init__(self):
        self.direct_uncertainty_sources = []

        #self.dependency_uncertainties = {}                                         ###!!!
        
        """ Handle initialization of non-inputs here! """
    
    def reset(self):
        for source in self.direct_uncertainty_sources:
            source.values = None
        self.root_weighted_uncertainties            = None
        self.root_upsample_factors                  = None
        self.total_uncertainty                      = None
        self.direct_uncertainties_calculated        = False
        self.total_uncertainty_calculated           = False
        self.is_certain                             = None
        
    def rescaleUncertaintySources(self, upsample_factor):
        """ Rescales the variance of all uncertainty sources to account for the partial time-aggregation/rebinning of the uncertainty parent variable """
        for source in self.direct_uncertainty_sources:
            if source.correlation == 1:
                source.sigma *= upsample_factor ; continue
            if source.correlation == 0:
                source.sigma *= np.sqrt(upsample_factor) ; continue
            else:
                u_vec = np.ones(upsample_factor) * source.sigma
                m_corr = np.ones((upsample_factor, upsample_factor)) * source.correlation
                np.fill_diagonal(m_corr, 1)
                source.sigma = u_vec.T @ m_corr @ u_vec
        
    
    def getSourceNames(self):
        """ Returns list of names of all direct uncertainty sources """
        return [u.name for u in self.direct_uncertainty_sources]
    
    def getSource(self, source_name):
        """ Tries to retrieve the uncertainty source with the given name from direct or root sources """
        source_dict = {source.name : source for source in self.direct_uncertainty_sources}
        if source_name in source_dict.keys():
            return source_dict[source_name]
        elif self.root_sources is not None:
            source_dict = {source.name : source for source in self.root_sources}
            if source_name in source_dict.keys():
                return source_dict[source_name]
        else:
            print(f"Was not able to retrieve uncertainty source of name: {source_name} from variable {self.var_name}")
            return None
        
    
    def getAbsoluteRootSplit(self, k):
        """ Gets absolute uncertainty split by root source, multiplied by a coverage factor k """
        if not self.is_calculated:
            raise ValueError("Tried to obtain absolute root split for variable {var_name}, uncertainty is not calculated yet.")
        if self.root_uncertainty_contribution_split is None:
            raise ValueError("Cannot get an absolute root split: please calculate the uncertainty root split first!")
        return np.multiply(self.root_uncertainty_contribution_split, self.total_uncertainty * k)
    
    def getRelativeRootSplit(self, k):
        """ Gets relative uncertainty split by root source, multiplied by a coverage factor k """
        absolute_split = self.getAbsoluteRootSplit(k)
        relative_uncertainty_split = np.divide(absolute_split,
                                               self.variable.values,
                                               out=np.zeros_like(self.root_uncertainty_contribution_split),
                                               where=(self.variable.values != 0))
        relative_uncertainty_split *= 100
        relative_uncertainty_split[np.where(relative_uncertainty_split>20)] = 20
        return relative_uncertainty_split
    
    def plotAbsoluteRootSplit(self, k):
        """ Plots absolute root split of the uncertainty, over time, with coverage factor k """
        absolute_split = self.getAbsoluteRootSplit(k)
        time_axis = self.variable.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        ax.stackplot(time_axis, *absolute_split, labels=self.root_uncertainty_sources)
        
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.xlabel("Time")
        plt.ylabel("Absolute error split")
        plt.show()
        
    def plotRelativeRootSplit(self, k):
        """ Plots relative root split of the uncertainty, over time, with coverage factor k """
        relative_split = self.getRelativeRootSplit(k)
        time_axis = self.variable.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        ax.stackplot(time_axis, *relative_split, labels=self.root_uncertainty_sources)
        
        ax.set_ylim(0,10)
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.xlabel("Time")
        plt.ylabel("Relative error split")
        plt.show()
    
            
        

        
class Variable:
    def __init__(self, name, description=None, values=None, is_basic=True, is_hardcoded=False, is_rate=None, aggregation_rule=None, \
                 equation=None, first_time=None, last_time=None, \
                     is_timesum=False, timesum_settings=None):
        self.name               = name              #str: variable name
        self.description        = description       #str: variable description
        self.is_basic           = is_basic          #bool: defines whether variable is basic or derived
        self.is_hardcoded       = is_hardcoded      #bool: defines whether the value is hard-coded in the input scripts (done for universal constants)
        self.is_rate            = is_rate           #bool: defines whether the quantity is a rate (quantity per unit time) or a quantity
        self.aggregation_rule   = aggregation_rule  #str: defines the quantity aggregation rule - depends on whether variable is extensive or intensive
        if self.is_rate and self.aggregation_rule=="average":
            raise ValueError(f"Variable {self.name} not well-defined: variable is defined as the rate of an intensive quantity over time")
        
        self.values             = values            #[int, float, array, None]: variable values
        self.partial_values     = {}                #dict: dictionary of values of the evaluated partial derivatives of this variable with respect to each dependency

        #If it is a derived variable then an equation and constituent variables should be passed
        if not self.is_basic:
            if equation is None:
                raise ValueError(f"{self.name} is defined as a derived variable, please include equation")
        self.equation           = equation          #str: defines variable equation
        self.dependency_names   = None              #list: lists variables in equation
        self.dependencies       = {}                #dict: dict of variables that are the direct dependencies
        self.is_root_consistent = False             #bool: whether the variable traces consistently to basic variables
        
        self.is_timesum         = is_timesum       #bool: defines whether variable is a timesum
        self.aggregation_step   = None             #float: timestep of the dependency over which the variable is time-integrated, required for uncertainty calculation
        self.non_aggregated_values = None          #array: values that are aggregated over - used in uncertainty calculation
        self.timesum_settings   = timesum_settings #list of options ###!!!
        
        self.sympy_symbol_map   = None      #dict: dictionary of sympy symbols
        self.sympy_equation     = None      #sympy interpretable of the variable equation
        self.executable         = None      #executable: sympy-built executable of the variable equation - excluding timesums
        self.partial_executables = None     #dict: dictionary of executables for the partial derivates of the variable for each dependency
        self.calculation_engine = None      #calculation engine: necessary for more complex calculations
                
        if first_time is not None and last_time is not None:
            self.addTimeStep([first_time, last_time])        
        else:
            self.start_time    = None              #datetime object: start time of the timeseries (== the time of the first datapoint - timestep)
            self.first_time    = None              #datetime object: time of the first datapoint
            self.last_time     = None              #datetime object: time of the last datapoint
            self.timestep      = None              #float:           number of seconds per timestep
        
        self.uncertainty       = VariableUncertainty(var_name=self.name, variable=self)
        
        self.harmonization_cache = None            #Dictionary of TimeHarmonizationData objects, stores how each dependency is time-harmonized to calculate this variable
        
        
    def __str__(self):
        if self.is_basic:
            return (f"Basic variable: {self.name} \nDescription: {self.description}")
        else:
            return (f"Derived variable: {self.name} = {self.equation} \nDescription: {self.description}")
    
    def __len__(self):
        if self.values is None:
            raise ValueError(f"Tried to obtain length of variable {self.name} for which no values are defined (yet).")
        elif isinstance(self.values,(float,int)):
            return 1
        else:
            return len(self.values)
    
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
            raise ValueError(f"Summing variables {self.name} and {other.name} failed. {self.name} has length {len(self.values)}, {other.name} has length {len(other.values)}")
    
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
     
    def reset(self):
        self.values                 = None
        self.partial_values         = {}
        self.aggregation_step       = None
        self.non_aggregated_values  = None
        self.start_time             = None
        self.first_time             = None
        self.last_time              = None
        self.timestep               = None
        self.harmonization_cache    = None
        self.uncertainty.reset()
    
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
            
    def getTimeAxis(self):
        if isinstance(self.values, (float, int)):
            raise ValueError("Cannot construct time axis for variable {self.name}, since the variable is not time-dependent.")
        n = len(self.values)
        return [self.first_time + timedelta(seconds=i * self.timestep) for i in range(n)]
    
    def getTimeData(self):
        if self.timestep is None:
            return None
        else:
            return (self.start_time, self.first_time, self.last_time, self.timestep)
    
    def setTimeData(self, time_data):
        if time_data is None:
            return
        elif len(time_data)==3:
            self.start_time, self.last_time, self.timestep = time_data
            self.first_time = self.start_time + timedelta(seconds=self.timestep)
        else:
            self.start_time, self.first_time, self.last_time, self.timestep = time_data 
            
    def printTimeData(self):
        if self.is_timesum:
            print(f"Variable {self.name} is a timesum, thus it has no time data.")
            return
        elif self.timestep is None:
            print(f"Variable {self.name} has no time data and is therefore either not initialized, or a constant.")
            return
        else:
            print(f"Variable {self.name} is a time series. \nTime range: {self.start_time} - {self.last_time} with timestep {self.timestep}. Total number of datapoints: {len(self.values)}.")
            return
            
    def addTimeStep(self, time_range):
        """ Given the time range (which are the times of first and last datapoint registrations) - calculates timestep length """
        self.first_time = time_range[0]
        self.last_time = time_range[1]
        
        if self.values is None or len(self.values)<2:
            raise ValueError(f"Automatic timestep calculation for variable {self.name} failed: no or too little values are specified (at least 2).")
       
        #Calculate timestep
        timestep = (self.last_time - self.first_time) / (len(self.values) - 1)
        
        self.start_time = self.first_time - timestep
        
        self.timestep = timestep.total_seconds()
        if not self.timestep.is_integer():
            raise ValueError(f"Timestep {self.timestep} is not integer, please ensure timesteps are an integer amount of seconds!")
        self.timestep = int(self.timestep)
    
        
    def executeEquation(self, store_results=True, force_recalculation=False, calculation_engine=None):
        """ Tries to call given or internal equation engine to execute variable equation """
        if calculation_engine is None:
            if self.calculation_engine is None:
                raise ValueError(f"Variable {self.name} cannot compute itself: no calculation engine given.")
            else:
                calculation_engine = self.calculation_engine
        return calculation_engine.executeVariableEquation(self, store_results=store_results, force_recalculation=force_recalculation)
    
    
    def giveReport(self, k=2, decimals=3, short_report=False):
        string = f"\nPresenting report of variable: {self.name}\n"
        if self.equation is not None:
            string += f"{self.name} = {self.equation}\n"
        if self.is_basic:
            string += "Basic Variable\n"
        else:
            string += "Derived Variable\n"
        

        if self.values is None:
            string += "No values calculated for this variable\nEnd of report\n"
            print(string)
            return            
        else:
            is_calculated = True
            report_values = np.round(self.values, decimals=decimals)
            
            if np.isscalar(report_values):
                is_scalar = True
                if isinstance(report_values, np.ndarray):
                    report_values = report_values[0]
            else:
                is_scalar = False
        
        if not self.uncertainty.total_uncertainty_calculated:
            string += f"Values: {report_values}\nUncertainty not calculated\nEnd of report\n"
            print(string)
            return
        else:
            if self.uncertainty.is_certain:
                string += f"Values: {report_values} \nVariable is completely certain\nEnd of report\n"
                print(string)
                return
            report_uncertainties = np.round(self.uncertainty.total_uncertainty * k, decimals=decimals)
            if is_scalar:
                if isinstance(report_uncertainties, np.ndarray):
                    report_uncertainties = report_uncertainties[0]
                string += f"Value: {report_values} +/- {report_uncertainties}\n"
            else:
                string += f"Values:         {report_values} \n"
                string += f"Uncertainties:  {report_uncertainties}\n"
            string += f"Coverage factor: k = {k}\n"
            if not short_report:
                string += "Uncertainty sources contributing to variable uncertainty:\n"
                for source in self.uncertainty.root_sources:
                    string += f"-{source.name}\n"
            string += "End of report\n"
        print(string)
        return
    
    
    
    def plotValuesAndUncertainty(self, k):
        if self.values is None or self.uncertainty.total_uncertainty_calculated is False:
            raise ValueError(f"Cannot plot values and uncertainties of variable {self.name} since they are not calculated yet.")
        if np.isscalar(self.values):
            print("Cannot plot values of variable {self.name} because it is a scalar.")
        
        time_axis = self.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        
        ax.plot(time_axis, self.values)
        ax.fill_between(time_axis, self.values - self.uncertainty.total_uncertainty*k, self.values + self.uncertainty.total_uncertainty*k, 
                        color='red', alpha=0.5)
        
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        plt.xlabel("Time")
        plt.ylabel(self.name + " [$W/m^2$]")
        plt.title("$G_{POA}$ " + f"Assumed values and uncertainty, k={k}")
        plt.show()
        
            
            
                
        
    

        
    








            
            