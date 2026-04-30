from dataclasses import dataclass
from typing import Union, Optional
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


    
@dataclass
class ParsedVariableData:
    """ 
    Used by the input parser to store parsed settings before creating a variable.
    
    csv_pointer: str or None
        Tells job handler which CSV/dataframe column this variable should take its values from.
    is_hardcoded: bool
        Tells job handler whether variable value is hardcoded. This means the variable not affected by an equation tree state reset.
    var_values: float or int or np.ndarray or None
        Value(s) to assign to variable.
    timestep: float or str or None
        Variable timestep, typically not used.
    description: float or None
        Variable description.
    is_basic_variable: bool
        Whether the variable is basic or derived.
    is_maskable: bool
        Whether the variable uncertainty is allowed to be ignored if the variable values are 0.
    is_rate: bool
        Whether the variable is a rate (quantity over time) or a quantity.
    aggregation_rule: str
        Whether variable aggregation is a sum or an average, distinguishes extensive and intensive quantities.
    equation: str
        Variable equation.
    uncertainties: list[UncertaintySource]
        List of uncertainty sources acting on the variable.
    """
    def __post_init__(self):
        self.csv_pointer = None         #str: tells job handler which CSV column this variable should take its values from
        self.is_hardcoded = False       #bool: tells job handler whether variable value is hardcoded. This means this variable is skipped during resets
        self.var_values = None          #[int, float, array]: variable values
        self.timestep = None            #[float, str]: timestep of variable, or setting ###!!!
        self.description = None         #str: variable description
        self.is_basic_variable = None   #bool: flags if variable is basic or derived
        self.is_maskable = True         #bool: flags if masking uncertainty of the variable is valued 0 is allowed for this variable, allowed by default
        self.is_rate = None             #bool: flags whether variable is a rate (quantity over time) or a quantity
        self.aggregation_rule = None    #str: aggregation rule of the quantity
        self.equation = None            #str: equation of the variable
        self.uncertainties = []         #list: contains all uncertainty sources
        

@dataclass
class TimeHarmonizationData:
    """ 
    Stores information on how to bring a dependency dep_var to the exact temporal dimensions (time range, time step) of the target_var.
    This object is populated by the TimeEngine and contains information on the temporal resolution upsample facot, number of timesteps that
    must be pruned at the end of timeseries to ensure identical start and end times. 
    
    Attributes
    ----------
    dep_var_name: Variable
        Name of the variable this harmonization concerns.
    base_timestep: float
        Original timestep (in seconds) of the variable.
    new_timestep: float
        New timestep (in seconds) of the variable.
    new_start_time: Datetime
        New start time of the variable.
    new_last_time: Datetime
        New last time of the variable.
    low_index: int
        Index at which rebinning should start in the original timeseries - to ensure temporal matching of step boundaries.
    high_index: int
        Index at which rebinning should end in the original timeseries - to ensure temporal matching of step boundaries.
    upsample_factor: int
        Difference factor between base timestep and new timestep.
    target_var_name: Variable or None
        Variable to which the temporal resolution is matched.
    prune_offset_start: int
        Number of new timesteps to discard at the front to ensure identical start times of the new timeseries as the target variable.
    prune_offset_end: int
        Number of new timesteps to discard at the end to ensure identical end times of the new timeseries as the target variable.
    new_values: float or np.ndarray or None
        Rebinned timeseries of the variable.
    """
    dep_var_name:       str         #Name of the affected variable
    base_timestep:      float       #Original timestep of the variable
    new_timestep:       float       #New timestep of the variable
    new_start_time:     datetime    #New start time
    new_last_time:      datetime    #New end time
    low_index:          int         #Index to start rebinning in original array
    high_index:         int         #Index to end rebinning in original array
    #low_fraction:       int         #If a new bin stretches from bin i to bin i+n, fraction of original bin i to include in new bin
    #high_fraction:      int         #Identical as above, but for bin i+n. high_fraction = 1-low_fraction: we assume integer upsample factors
    upsample_factor:    int         #Factor increase of new bin size compared to old bin size
    target_var_name:    Optional[str]   = None         #Name of the derived variable that required harmonization of dependencies
    prune_offset_start: Optional[int]   = None         #Number of new timesteps that are discarded at start to ensure a common start time
    prune_offset_end:   Optional[int]   = None         #Number of new timesteps that are discarded at end to ensure a common end time
    #smuggled_time:      Optional[float] = None        #Amount of seconds that  were smuggled during rebinning
    new_values:         Optional[Union[np.ndarray, float]] = None   #Can be used to store new values
        
    def getFirstTime(self):
        """ 
        Returns the new first time of the variable.
        
        Returns
        -------
        Timestamp
            New first time of the variable.
        """
        return self.new_start_time + timedelta(seconds=self.new_timestep)
    
    def getTotalOffsetSteps(self):
        """ 
        Returns the number of datapoints to discard at the lower and higher ends of the original timeseries 
        to match the start and end times of the target variable.
        
        Returns
        -------
        Tuple(int, int)
            Number of original datapoints to discard to ensure common start and end times of the original timeseries and the target variable.
        """
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
    """ 
    Dataclass containing all information of a single uncertainty source, including its name, characteristics and values.
    
    Attributes
    ----------
    name: str
        Name of the uncertainty source
    is_relative: bool
        Whether the uncertainty is defined as relative to (a percentage of) the variable value, or an absolute value.
    is_symmetric: bool
        Whether the uncertainty distribution is symmetric or one-sided.
    shape: str
        Distribution shape of the uncertainty source (normal, rectangular, triangular, u-shaped).
    correlation: int or float
        Temporal autocorrelation of the uncertainty.
    multiplier: str or None
        Function that the uncertainty is multiplied with. Used for the irradiance and zenith angle dependence of the directional response error.
    values: float or int or np.ndarray or None
        Absolute uncertainty value(s) assumed by this uncertainty source.
    parent_variable: Variable or None
        Variable which this uncertainty source acts upon.
    sigma: float or int or np.ndarray or None
        Standard deviation of the distribution. Mutually exclusive with bound.
    bound: float or int or np.ndarray or None
        Uncertainty limit of the distribution (or 99% percentile for normal distributions). Mutually exclusive with sigma.
    
    
    """
    name:               str
    is_relative:        bool
    is_symmetric:       bool
    shape:              str
    correlation:        Union[float, int, np.ndarray]
    multiplier:         Optional[str] = None
    values:             Optional[Union[float, int, np.ndarray]] = None
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
            raise RuntimeError(f"Either a standard deviation or a bound should be provided for uncertainty {self.name}, neither is provided.")
        if (self.sigma is not None and self.bound is not None):
            raise RuntimeError(f"Either a standard deviation or a bound should be provided for uncertainty {self.name}, not both!")
        
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
        """ 
        Constructs correlation matrix of requested size based upon internal correlation information.
        
        Parameters
        ----------
        size: int
          Size of the requested correlation matrix.
        
        Returns
        np.ndarray
            2D correlation matrix.
        """
        M = np.ones((size, size)) * self.correlation
        np.fill_diagonal(M, 1)
        return M
            

@dataclass
class VariableUncertainty:
    """ Dataclass storing all uncertainty information related to a specific variable.
    
    Contains a registry of uncertainty sources acting on the variable and informative booleans about uncertainty state and calculation rules.
    Upon calculation of variable uncertainty, fields in this dataclass are populated with all uncertainty sources influencing this variable's uncertainty,
    the variable's total uncertainty, and numerous intermediate results required for propagating uncertainty uptree.
    
    Parameters
    ----------
    var_name: str
        Name of Variable this instance is owned by.
    variable: Variable
        Pointer to Variable object this instance is owned by.
        
    Attributes
    ----------
    direct_uncertainty_sources: list or None
        List of UncertaintySource objects acting on the variable this instance is owned by.
    root_sources: list or None
        List of UncertaintySource objects acting on this variable, or any variable downtree.
    root_total_upsample_factors: list or None
        List of lists, containing integer total upsample factors, which are the integer factors of a variable's ``timestep`` compared to the root variable's ``timestep``.
    root_local_upsample_factors: list or None
        List of lists, containing integer local upsample factors to root, which are the integer factors of a variable's ``timestep`` compared to the previous variable's ``timestep``.
    root_weighted_uncertainties: dict or None
        Dictionary containing the weighted uncertainties at root temporal resolution for each uncertainty source.
    root_upsample_factors: list or None
        List of ints, containing the total upsample factor between this variable's and the uncertainty root variable's timesteps.
    aggregated_weighted_uncertainties: np.ndarray or None
        Weighted uncertainty values of all uncertainty sources, potentially partially aggregated, in the temporal resolution of this variable, corresponding to the variable's ``values``.
    total_uncertainty: np.ndarray, float or None
        Total uncertainty of this variable, either as a single value or as a timeseries.
    is_masked: bool
        Whether masking of uncertainties where the varaible is zero-valued was applied to calculate this variable uncertainty.
    is_certain: bool
        Whether this variable has uncertainty or if it is completely certain.
    """
    var_name:                               str
    variable:                               "Variable" #Placeholder!
    
    direct_uncertainty_sources:             Optional[list] = None
    
    root_sources:                           Optional[list] = None
    root_total_upsample_factors:            Optional[list] = None
    root_local_upsample_factors:            Optional[list] = None
    root_weighted_uncertainties:            Optional[dict] = None
    root_upsample_factors:                  Optional[list] = None
    
    aggregated_weighted_uncertainties:      Optional[np.ndarray] = None

    total_uncertainty:                      Optional[Union[np.ndarray, float]] = None
    
    direct_uncertainties_calculated:        Optional[bool] = False
    total_uncertainty_calculated:           Optional[bool] = False
    is_masked:                              Optional[bool] = False
    is_certain:                             Optional[bool] = None
    
    def __post_init__(self):
        if self.direct_uncertainty_sources is None:
            self.direct_uncertainty_sources = []
    
    def reset(self):
        """ 
        Resets all stateful information of the variable uncertainty.
        """
        for source in self.direct_uncertainty_sources:
            source.values = None
        self.root_weighted_uncertainties            = None
        self.root_total_upsample_factors            = None
        self.root_local_upsample_factors            = None
        self.total_uncertainty                      = None
        self.direct_uncertainties_calculated        = False
        self.total_uncertainty_calculated           = False
        self.is_certain                             = None
        
    def rescaleUncertaintySources(self, upsample_factor):
        """ Rescales the variance of all uncertainty sources to account for the partial time-aggregation/rebinning of the uncertainty parent variable. Destructive operation.
        
        Parameters
        ----------
        upsample_factor: int
            Integer upsample factor by which the timestep increases.
        """
        for source in self.direct_uncertainty_sources:
            if source.correlation == 1:
                source.sigma *= upsample_factor ; continue
            if source.correlation == 0:
                source.sigma *= np.sqrt(upsample_factor) ; continue
            else:
                ###!!!! Replace with uncertainty_source.getCorrelationMatrix function
                u_vec = np.ones(upsample_factor) * source.sigma
                m_corr = np.ones((upsample_factor, upsample_factor)) * source.correlation
                np.fill_diagonal(m_corr, 1)
                source.sigma = u_vec.T @ m_corr @ u_vec
        
    
    def getDirectSourceNames(self):
        """ Returns list of names of all direct uncertainty sources.
        
        Returns
        -------
        list[str]
            List of names of uncertainty sources acting directly on this variable.
        """
        return [u.name for u in self.direct_uncertainty_sources]
    
    def getSourceNames(self):
        """ Returns list of names of all uncertainty sources. Only works if uncertainty of this variable is calculated.
        
        Returns
        -------
        list[str]
            List of names of uncertainty sources acting on this variable, or any variable downtree.
        """
        return [u.name for u in self.root_sources]
    
    def getSource(self, source_name):
        """ Tries to retrieve the uncertainty source with the given name from direct or root sources.
        
        Parameters
        ----------
        source_name: str
            Name of the uncertainty source to be retrieved.
        
        Returns:
        UncertaintySource or None
            If found, returns the uncertainty source asked for, otherwise returns None and logs failure to console.
        """
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
        
        
class Variable:
    """A physical variable or auxiliary node in the equation tree.
    
    A variable is either a basic variable (at the base of the equation tree) or a derived variable.
    Derived variables are defined through an equation. Variables can also be time-aggregates of other variables.
    Variables can be given uncertainty sources, but are also affected by any uncertainties downtree.
    Most regular arithmetic (including exponents and numpy operations) is supported for variables, and acts on the variable values.
    
    Parameters
    ----------
    name : str
        Variable name.
    description : str
        Optional variable description.
    values : int or float or array-like or None
        Value of the variable, or a timeseries of values, or None if the values are to be calculated.
    is_basic : bool
        Flags whether variable is at the root of the equation tree.
    is_hardcoded : bool
        Flags whether the variable is a hardcoded constant. Ensures value is not wiped upon equation tree reset.
    is_maskable : bool
        Flags whether the uncertainty in this variable can be ignored if its value is 0 (e.g. for irradiance measurement at night).
    is_rate : bool
        Flags whether the variable is a quantity (e.g. energy) or a rate of a quantity over time (e.g. power).
    aggregation_rule : str or None
        Defines variable aggregation rule (sum or average) - whether a variable is an intensive or extensive quantity
    equation : str or None
        Equation defining the variable (required if not a basic variable).
    is_timesum : bool
        Defines whether a variable is an aggregation over time. Note: this flag is exclusive to auxiliary timesum variables created by the equation engine.
    
    Raises
    ------
    RuntimeError
        If the variable is defined as both a rate and an intensive quantity.
    RuntimeError
        If the variable is defined as a derived variable but no equation is provided.
    
    Attributes
    ----------
    values : float or int or array_like or None
        Stored variable values
    equation : str or None
        Equation used by SymPy to calculate the values of derived variables. Required for derived variables.
    dependencies : dict
        Dictionary relating variables in the equation to other Variable instances
    uncertainty : VariableUncertainty
        Object storing all uncertainty information related to the variable.
    harmonization_cache : dict or None
        Dictionary of TimeHarmonizationData objects, storing how dependencies are time-harmonized during calculation. Used when dependencies have different temporal resolution.
    """
    def __init__(self, name, description=None, values=None, is_basic=True, is_hardcoded=False, is_maskable=True, 
                 is_rate=None, aggregation_rule=None, equation=None, is_timesum=False):
        self.name               = name              #str: variable name
        self.description        = description       #str: variable description
        self.is_basic           = is_basic          #bool: defines whether variable is basic or derived
        self.is_hardcoded       = is_hardcoded      #bool: defines whether the value is hard-coded in the input scripts (done for universal constants)
        self.is_maskable        = is_maskable       #bool: defines whether the variable uncertainty is allowed to be ignored if it's value at that time is 0 (e.g. True for G or Pout, False for Temperature)             
        self.is_rate            = is_rate           #bool: defines whether the quantity is a rate (quantity per unit time) or a quantity
        self.aggregation_rule   = aggregation_rule  #str: defines the quantity aggregation rule - depends on whether variable is extensive or intensive
        if self.is_rate and self.aggregation_rule=="average":
            raise RuntimeError(f"Variable {self.name} not well-defined: variable is defined as the rate of an intensive quantity over time")
        
        self.values             = values            #[int, float, array, None]: variable values
        self.partial_values     = {}                #dict: dictionary of values of the evaluated partial derivatives of this variable with respect to each dependency

        #If it is a derived variable then an equation and constituent variables should be passed
        if not self.is_basic:
            if equation is None:
                raise RuntimeError(f"{self.name} is defined as a derived variable, please include equation")
        self.equation           = equation          #str: defines variable equation
        self.dependency_names   = None              #list: lists variables in equation
        self.dependencies       = {}                #dict: dict of variables that are the direct dependencies
        self.is_root_consistent = False             #bool: whether the variable traces consistently to basic variables
        
        self.is_timesum         = is_timesum        #bool: defines whether variable is a timesum
        self.aggregation_step   = None              #float: timestep of the dependency over which the variable is time-integrated, required for uncertainty calculation
        self.non_aggregated_values = None           #array: values that are aggregated over - used in uncertainty calculation
        
        self.sympy_symbol_map   = None              #dict: dictionary of sympy symbols
        self.sympy_equation     = None              #sympy interpretable of the variable equation
        self.executable         = None              #executable: sympy-built executable of the variable equation - excluding timesums
        self.partial_executables = None             #dict: dictionary of executables for the partial derivates of the variable for each dependency
        self.calculation_engine = None              #calculation engine: necessary for more complex calculations
        
        self.start_time         = None              #datetime object: start time of the timeseries (== the time of the first datapoint - timestep)
        self.first_time         = None              #datetime object: time of the first datapoint
        self.last_time          = None              #datetime object: time of the last datapoint
        self.timestep           = None              #float:           number of seconds per timestep
        
        self.uncertainty        = VariableUncertainty(var_name=self.name, variable=self)
        
        self.harmonization_cache = None            #Dictionary of TimeHarmonizationData objects, stores how each dependency is time-harmonized to calculate this variable
        
        
    def __str__(self):
        if self.is_basic:
            return (f"Basic variable: {self.name} \nDescription: {self.description}")
        else:
            return (f"Derived variable: {self.name} = {self.equation} \nDescription: {self.description}")
    
    def __len__(self):
        if self.values is None:
            raise RuntimeError(f"Tried to obtain length of variable {self.name} for which no values are defined (yet).")
        elif isinstance(self.values,(float,int)):
            return 1
        else:
            return len(self.values)
    
    def __neg__(self):
        return -self.values
    
    def __add__(self, other):
        #If floats or ints are involved the logic is simple
        if isinstance(other, (int, float, np.ndarray)):
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
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        #If floats or ints are involved the logic is simple
        if isinstance(other, (int, float, np.ndarray)):
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
     
    def __rsub__(self, other):
        return -1*self.__sub__(other)
     
    def __mul__(self, other):
        #If floats or ints are involved the logic is simple
        if isinstance(other, (int, float, np.ndarray)):
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
        if isinstance(denominator, (int, float, np.ndarray)):
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
        if isinstance(numerator, (int, float, np.ndarray)):
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
    
    def __rpower__(self, base):
        if hasattr(base, "values"):
            base = base.values
        if isinstance(base, (int, float)) or isinstance(self.values, (int, float)):
            return base ** self.values
        else:
            if len(self.values) == len(base):
                return base ** self.values
            else:
                raise ValueError(f"Error: taking power with variable {self.name} failed, trying to take power with arrays of shape {len(base)} and {len(self.values)}.")
    
    def exp(self):
        """ Take exponential of own values """
        return np.exp(self.values)
    
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """ Numpy ufunc overrider """
        if ufunc in (np.add, np.subtract, np.multiply, np.divide, np.power):
            if ufunc == np.add:
                return inputs[1].__radd__(inputs[0])
            if ufunc == np.subtract:
                return inputs[1].__rsub__(inputs[0])
            if ufunc == np.multiply:
                return inputs[1].__rmul__(inputs[0])
            if ufunc == np.divide:
                return inputs[1].__rtruediv__(inputs[0])
            if ufunc == np.power:
                return inputs[1].__rpow__(inputs[0])
        return NotImplemented
    
    
    def reset(self):
        """ 
        Resets stateful information of the variable, including timestep information, uncertainty information and any cached time harmonizations.
        """
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
    
    def hasValues(self):
        """ 
        Method to check whether the variable has any values defined.
        
        Returns
        -------
        bool
            Whether the variable has any defined values.
        """
        if self.values is None:
            return False
        else:
            return True
        
    def hasDirectUncertainties(self):
        """ 
        Method to check whether there are any uncertainty sources acting directly on this variable.
        
        Returns
        -------
        bool
            Whether the variable has any direct uncertainty sources acting on it.   
        """
        return self.uncertainty.direct_uncertainty_sources is not None
    
    def addUncertaintySource(self, uncertainty):
        """ 
        Method to add an UncertaintySource, or a list of them, as an uncertainty acting on this variable.
        
        Parameters
        ----------
        uncertainty: VariableUncertainty or list[VariableUncertainty]
            Uncertainty source(s) to add to the registry of direct uncertainty sources of this variable.
        """
        #Can be a single uncertainty source or a whole list of them
        if type(uncertainty) is list:
            self.uncertainty.direct_uncertainty_sources.extend(uncertainty)
        else:
            self.uncertainty.direct_uncertainty_sources.append(uncertainty)
            
    def getTimeAxis(self, times_only=False):
        """ 
        Routine to obtain a full list of timestamps corresponding to the variable values.
        
        Parameters
        ----------
        times_only: bool
            Whether the times returned are only times, or should include a date as well.
            
        Returns
        -------
        list
            List of timestamps.    
        """
        if isinstance(self.values, (float, int)):
            raise RuntimeError("Cannot construct time axis for variable {self.name}, since the variable is not time-dependent.")
        n = len(self.values)
        if times_only:
            return [(self.first_time + timedelta(seconds=i * self.timestep)).time() for i in range(n)]
        else:
            return [self.first_time + timedelta(seconds=i * self.timestep) for i in range(n)]
    
    def getTimeData(self):
        """ 
        Routine to obtain timedata of the variable, if it has timedata.
        
        Returns
        -------
        None or Tuple(Timestamp, Timestamp, Timestamp, float)
            If variable has no timedata, returns None.
            If variable has timedata, returns start time, first time, last time and timestep (in seconds) of the timeseries.
        """
        if self.timestep is None:
            return None
        else:
            return (self.start_time, self.first_time, self.last_time, self.timestep)
    
    def setTimeData(self, time_data):
        """ 
        Routine to populate the timedata of a variable from a tuple completely describing the timedata. 
        
        Parameters
        ----------
        timedata: None or tuple[Timestamp, Timestamp, float] or tuple[Timestamp, Timestamp, Timestamp, float]
            If timedata is None, all timedata will be set to None.
            If timedata is a tuple of length 3, it is assumed the start_time, end_time and timestep (in seconds) are given; first_time is calculated.
            If timedata is a tuple of length 4, it is assumed start_time, first_time, end_time and timestep (in seconds) are given. 
        """
        if time_data is None:
            self.start_time, self.first_time, self.last_time, self.timestep = None, None, None, None
            return
        elif len(time_data)==3:
            self.start_time, self.last_time, self.timestep = time_data
            self.first_time = self.start_time + timedelta(seconds=self.timestep)
        else:
            self.start_time, self.first_time, self.last_time, self.timestep = time_data 
            
    def printTimeData(self):
        """ 
        Prints information about the temporal characteristics of the variable. 
        """
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
        """ 
        Given the time range (which are the times of first and last datapoint registrations) - calculates timestep length and populates timedata.
        
        Parameters
        ----------
        time_range: tuple[Timestamp, Timestamp]
            tuple containing the first_time and last_time of the variable timeseries   
            
        Raises
        ------
        ValueError
            If timestep is not an integer amount of seconds. Non-integer timestep lengths are assumed to be indicative of missing datapoints.
        """
        self.first_time = time_range[0]
        self.last_time = time_range[1]
        
        if self.values is None or len(self.values)<2:
            raise ValueError(f"Automatic timestep calculation for variable {self.name} failed: no or too little values are specified (at least 2).")
       
        #Calculate timestep
        timestep = (self.last_time - self.first_time) / (len(self.values) - 1)
        
        self.start_time = self.first_time - timestep
        
        self.timestep = timestep.total_seconds()
        if not self.timestep.is_integer():
            raise ValueError(f"Setting timestep for {self.name} failed: timestep {self.timestep} is not integer, please ensure timesteps are an integer amount of seconds!")
        self.timestep = int(self.timestep)
    
        
    def executeEquation(self, store_results=True, force_recalculation=False, calculation_engine=None):
        """ 
        Tries to call given or internal calculation engine to execute variable equation. Conditionally also updates own values (default True). 
        
        Parameters
        ----------
        store_results: bool, default = True
            Whether the results should be stored as the variable values.
        force_recalculation: bool, default = False
            Whether the values should be recalculated if the variable already has values defined.
        calculation_engine: None or CalculationEngine, default = None
            Calculation engine to execute the equation, optional if the variable has an internal calculation engine given.
            
        Raises
        ------
        RuntimeError
            If the variable has no calculation engine and no calculation engine is given.
            
        Returns
        -------
        float or array
            Calculated values of the variable.
        """
        if calculation_engine is None:
            if self.calculation_engine is None:
                raise RuntimeError(f"Variable {self.name} cannot compute itself: no calculation engine given.")
            else:
                calculation_engine = self.calculation_engine
        return calculation_engine.executeVariableEquation(self, store_results=store_results, force_recalculation=force_recalculation)
    
    
    def giveReport(self, k=2, decimals=3, short_report=False):
        """ 
        Prints a report on the equation and state of the variable.
        
        Parameters
        ----------
        k: int or float, default = 2
            Coverage factor at which to present uncertainty bounds. k=2 provides an approximate 95% coverage interval.
        decimals: int, default = 3
            Number of decimals at which to report numerical values.
        short_report: bool, default = False
            Whether to print a list of uncertainty sources acting on this variable.
        """
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
    
    
    
    def plotValues(self):
        """ 
        Plot the values assumed by this variable over time.
        
        Raises
        ------
        RuntimeError
            If the variable values are not defined.
        """
        if self.values is None:
            raise RuntimeError(f"Cannot plot values of variable {self.name} since they are not calculated yet.")
        if np.isscalar(self.values):
            print("Cannot plot values of variable {self.name} because it is a scalar.")
            
        time_axis = self.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        
        ax.plot(time_axis, self.values)
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        plt.xlabel("Time", fontsize=14)
        plt.ylabel(self.name, fontsize=14)
        plt.show()
    
    
    def plotValuesAndUncertainty(self, k=2):
        """ 
        Plot the values and indicate associated uncertainty of the variable over time. 
        
        Parameters
        ----------
        k: int or float, default=2
            Coverage factor at which to show the uncertainty interval.
        
        Raises
        ------
        RuntimeError
            If the variable's values or total uncertainties are not defined.
        """
        if self.values is None or self.uncertainty.total_uncertainty_calculated is False:
            raise RuntimeError(f"Cannot plot values and uncertainties of variable {self.name} since they are not calculated yet.")
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
        plt.ylabel(self.name)
        plt.show()
        
            
            
                
        
    

        
    








            
            