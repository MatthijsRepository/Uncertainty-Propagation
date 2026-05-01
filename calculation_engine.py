from my_dataclasses import Variable #, TimeHarmonizationData
#from time_engine import TimeEngine
#import datetime
import numpy as np
import copy


    
class CalculationEngine:
    """ 
    This engine handles the evaluation of the equation tree through recursive executable execution. Also executes partial derivative executables.
    Engine invokes time engine in case temporal range or resolution of variabales does not match.
    
    Attributes
    ----------
    time_engine: TimeEngine
        Time engine that is invoked to ensure temporal harmony between input data and build time harmonization caches for variables.
    equation_engine: EquationEngine or None
        Equation engine that is invoked to build equation executables if these are not in place.
    """
    def __init__(self, time_engine, equation_engine=None):
        self.time_engine = time_engine
        self.equation_engine = equation_engine
        return
    
    def validateBasicVariables(self, variables, equation_engine=None):
        """ 
        Ensures all basic variables in a given variables dictionary have their values in place.
        If a basic variable has no values, attempts to calculate their values from upstream variables that do have their values defined.
        Invokes an equation engine to build the executables to do so, if these are not already in place.
        
        Parameters
        ----------
        variables: dict[str, Variable]
            Dictionary of variables for which the basic variables must be validated.
        equation_engine: EquationEngine or None
            Equation engine invoked to build equation executables if these are required and not already in place.
        
        Raises
        ------
        ValueError
            If a basic variable has no values, but also has no equation to calculate values with.
        """
        if equation_engine is None:
            equation_engine = self.equation_engine
        
        for var in variables.values():
            if var.is_basic and var.values is None:
                if var.equation is None:
                    raise ValueError(f"Basic variable {var.name} has no defined values and is not defined by an equation, no values can be given to this variable.")
                else:
                    #Populate variable dependencies
                    equation_engine.populateVariableDependencies(var, variables=variables)
                    #Generate variable callable
                    equation_engine.buildVariableExecutable(var)
                    #Calculate variable values, var.values is populated
                    values, timedata, non_aggregated_values, aggregation_step, harmonized_data = self._executeVariableEquation(var)
                    var.values   = values
                    var.timedata = timedata                    
    
    def evaluateVariable(self, var, update_var=True, calculate_dependencies=True, force_recalculation=False, silent=True, indent=""):
        """ 
        Recursively calculates the values of a desired variables.
        Variable ``values`` attribute is updated with the calculated values by default.
        Recursion occurs depth-first, node traversal stops once:
            - A variable with defined values is encountered, if ``force_recalculation=False``.
            - A basic variable is encountered, if ``force_recalculation=True``.
        
        Parameters
        ----------
        var: Variable
            Variable to evaluate values for.
        update_var: bool, default=True
            Whether to update the ``values`` of the variables with the calculated values.
        calculate_dependencies: bool, default=True
            Whether to recursively calculate dependency values if these are not yet defined.
        force_recalculation: bool, default=False
            Whether to recalculate variables whose ``values`` are already defined, unless they are basic variables.
        silent: bool, default=True
            Whether to print the recursive path during execution.
        indent: str, default=""
            Print indent for reading clarity of execution path printing.
        
        Raises
        ------
        ValueError
            If the variable has no equation or executable defined.
        ValueError
            If the ``values`` of a dependency are not defined but ``calculate_dependencies=False`` prevents recursion.
        
        Returns
        -------
        float or np.ndarray
            Calculated values of ``var``.
        """
        if not silent:
            print(f"{indent}Calculating values of variable {var.name} with dependencies {var.dependency_names}")
        
        #Check if the variable equation and executable are in place
        if var.equation is None:
            raise ValueError(f"Tried to execute the equation of variable {var.name}, for which no equation is defined.")
        elif var.executable is None:
            raise ValueError(f"Tried to execute the equation of variable {var.name} = {var.equation}, but no equation executable has been built for this variable.")
        
        #Check if the dependencies have been calculated, and optionally recursively calculate these
        for dep in var.dependencies.values():
            if dep.values is not None:
                continue
            elif calculate_dependencies is False:
                raise ValueError(f"Calculation of variable {var.name} failed: dependency {dep.name} has no values defined and automatic dependency calculation is turned off.")
            else:
                self.evaluateVariable(dep, update_var=update_var, force_recalculation=force_recalculation, silent=silent, indent=f"   {indent}")
        
        #Execute variable equation
        values, timedata, non_aggregated_values, aggregation_step, harmonized_data = \
            self._executeVariableEquation(var, force_recalculation=force_recalculation)
        
        #Potentially update variable
        if update_var:
            var.values = values
            var.setTimeData(timedata)
            var.non_aggregated_values = non_aggregated_values ###!!! Can cause duplication of data!
            var.aggregation_step = aggregation_step
            var.harmonization_cache = harmonized_data
        
        if not silent:
            print(f"{indent}Calculation of {var.name} complete, values: {var.values}")
        return values   
    
    def _executeVariableEquation(self, var, force_recalculation=False):
        """ 
        Helper function that handles calculation of variable's values.
        
        Invokes time engine to ensure dependencies are time-harmonious, executes equation and handles timesums, 
        returns calculated values and time harmonization metadata.
        
        Parameters
        ----------
        var: Variable
            Variable to calculate values of.
        force_recalculation: bool, default=False
            Whether to recalculate the time harmonization cache, even if it is already present.
        
        Returns
        -------
        tuple
            Calculation results and metadata.
            Tuple of structure
            (
                np.ndarray or float,
                tuple,
                np.ndarray or None,
                float or None,
                dict[str, TimeHarmonizationData]
            )
            - array or scalar, calculated values of the variable.
            - tuple containing start, end times and timestep of the variable.
            - array of non-aggregated values if variable is a timesum, None otherwise.
            - float of the timestep used during time-aggregation.
            - dictionary containing ``TimeHarmonizationData`` metadata objects for each dependency of the variable.
        """
        #If values already defined: do nothing unless forced recalculation is desired
        if var.values is not None:
            if force_recalculation is True:
                print(f"WARNING: executing equation of variable {var.name} while values are already defined!")
            else:
                return var.values
        
        args, timedata, harmonized_data = self.time_engine.ensureDependencyTimeHarmony(var, force_recalculation=force_recalculation)

        calculated_values = var.executable(*args)   
        
        #Time aggregation if the variable is a timesum
        aggregation_step      = None #Aggregation step is the timestep of the timesummed data. This value is required for uncertainty calculation and therefore passed to the variable
        non_aggregated_values = None #Calculated values before aggregation - required for uncertainty calculation
        if var.is_timesum:
            non_aggregated_values = copy.deepcopy(calculated_values)
            calculated_values     = self._timeSum(var, calculated_values, timedata) ###!!!
            aggregation_step      = timedata[-1]
            timedata              = None
        
        return calculated_values, timedata, non_aggregated_values, aggregation_step, harmonized_data
    
    def _timeSum(self, var, calculated_values, timedata=None):
        """
        Handles the timesum calculation during variable evaluation.
        
        Parameters
        ----------
        var: Variable
            Variable for which the timesum is calculated.
        timedata: tuple[.., float or int]
            Timedata tuple containing the timestep as last entry.
        
        Returns
        -------
        float
            Time-aggregate of the variable values.
        """
        
        if var.aggregation_rule == "sum":
            calculated_values = np.sum(calculated_values)
        elif var.aggregation_rule == "average":
            calculated_values = np.average(calculated_values)
        else:
            raise ValueError(f"Timesum failed: no correct aggregation rule defined for variable {var.name}.")
            
        #Handle integration: if the variable is extensive AND defined as a rate over time, we can multiply by the timestep to obtain a quantity
        if var.is_rate:
            #Preferably use passed timedata, fallback to variable timedata if not available.
            if timedata is not None:
                timestep = timedata[-1]
            else:
                timestep = var.timestep
            calculated_values *= timestep
        
        return calculated_values
    
    def executePartialDerivative(self, var, dep_name, absolute_values=False, store_results=True, force_recalculation=False, equation_engine=None):
        """ 
        Executes the partial derivative executable of a variable for a given dependency, optionally also stores values in variable itself.
        
        Parameters
        ----------
        var: Variable
            Variable for which partial derivative is to be calculated.
        dep_name: str
            Name of the dependency to which the partial derivative should be taken.
        absolute_values: bool, default=False
            Whether the absolute value of results should be taken. Set to ``True`` for uncertainty calculations.
        store_results: bool, default=True
            Whether calculation results should be stored in ``var.partials_dict``.
        force_recalculation: bool, default=False
            Whether results should be recalculated if they are found to already exist.
        equation_engine: EquationEngine or None
            Equation engine to be used to build partial derivative executables. Optional if calculation engine already has one.
        
        Raises
        ------
        ValueError
            If the partial derivatives executables have not yet been built and no equation engine can be found to resolve this.
        
        Returns
        -------
        np.ndarray or float
            Values of the partial derivative of the variable with respect to the dependency.
        """
        #If no forced recalculation and if the values are already calculated we simply return the already calculated values
        if force_recalculation is False and dep_name in var.partial_values:
            return var.partial_values[dep_name]
        #If there is no executable for this dependency we raise an error
        if var.partial_executables is None:
            #Try to find a calculation engine
            if equation_engine is None:
                equation_engine = self.equation_engine
                if equation_engine is None:
                    raise ValueError(f"Tried to evaluate partial derivative of variable {var.name} while partial derivative executables have not been built yet.")
            #Use equation engine to build partial derivative executables
            equation_engine.buildPartialDerivativeExecutables(var)
                
        #Get partial executable, arguments; calculate values
        partial_executable = var.partial_executables[dep_name]
        args, timedata, harmonized_data = self.time_engine.ensureDependencyTimeHarmony(var, force_recalculation=force_recalculation)
        
        calculated_values = partial_executable(*args)
            
        #If the result is an array it may be of dtype 'object', here we cast the array to float
        if isinstance(calculated_values, np.ndarray):
            calculated_values = calculated_values.astype(float)
        
        #Here we catch shape mismatches in case of trivial derivatives that evaluate to a constant. 
        #In this case the partial derivatives may not match the shape of the dependency in question. These shapes must match for the uncertainty calculation
        #We must extract the expected length from the inputs. The args list follows the same order as the dependency_names list, so we can use the same index
        if np.isscalar(calculated_values):
            target_length = len(args[var.dependency_names.index(dep_name)])
            if target_length > 1:
                calculated_values = np.full(target_length, calculated_values, dtype=float)
        
        #In case of a trivial equation, calculated values will be a Variable object. Here we ensure we return numerical values.
        if isinstance(calculated_values, Variable):     ###!!! change this to be handled through an equation engine wrapper
            calculated_values = calculated_values.values
        
        #Take absolute values; i.e. obtain the sensitivities
        if absolute_values:
            calculated_values = np.abs(calculated_values)
        
        #Optionally store results
        if store_results:
            var.partial_values[dep_name] = calculated_values
        return calculated_values
    
    def executeAllPartials(self, var, absolute_values=False, store_results=True, force_recalculation=False, equation_engine=None):
        """ 
        Evaluates the partial derivatives with respect to each dependency of a variable.
        
        Parameters
        ----------
        var: Variable
            Variable for which the partial derivatives are to be calculated.
        absolute_values: bool, default=False
            Whether the absolute value of results should be taken. Set to ``True`` for uncertainty calculations.
        store_results: bool, default=True
            Whether calculation results should be stored in ``var.partials_dict``.
        force_recalculation: bool, default=False
            Whether results should be recalculated if they are found to already exist.
        equation_engine: EquationEngine or None
            Equation engine to be used to build partial derivative executables. Optional if calculation engine already has one.
        
        Returns
        -------
        dict[str, float or np.ndarray]
            Dictionary of dependency names and the values of their corresponding partial derivatives.
        """
        partial_values = {}
        for dep_name in var.dependency_names:
            partial_values[dep_name] = self.executePartialDerivative(var, dep_name, 
                                                                     absolute_values     = absolute_values, 
                                                                     store_results       = store_results, 
                                                                     force_recalculation = force_recalculation,
                                                                     equation_engine     = equation_engine)
        return partial_values
    
    
    



        
        
        


