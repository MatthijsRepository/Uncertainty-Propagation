from my_dataclasses import Variable, TimeHarmonizationData
from time_engine import TimeEngine
import datetime
import numpy as np


    
class CalculationEngine:
    def __init__(self, variables, time_engine, equation_engine=None):
        self.variables = variables
        self.time_engine = time_engine
        self.equation_engine = equation_engine
        return
    
    def validateBasicVariables(self, equation_engine, variables=None):
        """ Validates whether all basic variables are well-defined, and if not, uses an equation engine to try to calculate the values using the provided equation """
        if variables is None:
            variables = self.variables
        
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
                    self.executeVariableEquation(var, store_results=True)
    
    def evaluateVariable(self, var, update_var=True, calculate_dependencies=True, force_recalculation=False, silent=True, indent=""):
        """ Checks if the values of all dependencies are calculated, optionally calculates dependencies, and then calculates values of specified variable """
        if not silent:
            print(f"{indent}Calculating values of variable {var.name} with dependencies {var.dependency_names}")
        for dep in var.dependencies.values():
            if dep.values is not None:
                continue
            elif calculate_dependencies is False:
                raise ValueError(f"Calculation of variable {var.name} failed: dependency {dep.name} has no values defined and automatic dependency calculation is turned off.")
            else:
                self.evaluateVariable(dep, update_var=update_var, silent=silent, indent=f"   {indent}")
                
        #values = var.executeEquation(store_results=update_var, calculation_engine=self)
        values = self.executeVariableEquation(var, store_results=update_var, force_recalculation=force_recalculation)

        if not silent:
            print(f"{indent}Calculation of {var.name} complete, values: {var.values}")
        return values   
    
    def executeVariableEquation(self, var, store_results=True, force_recalculation=False, integrate=True):
        if var.equation is None:
            raise ValueError(f"Tried to execute the equation of variable {var.name}, for which no equation is defined.")
        elif var.executable is None:
            raise ValueError(f"Tried to execute the equation of variable {var.name} = {var.equation}, but no equation executable has been built for this variable.")
        
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
            non_aggregated_values = calculated_values
            calculated_values     = self.timeSum(var, calculated_values, timedata, integrate=integrate) ###!!!
            aggregation_step      = timedata[-1]
            timedata              = None
        
        #Optionally store the result as the new variable values
        if store_results:
            var.values = calculated_values
            var.setTimeData(timedata)
            var.non_aggregated_values = non_aggregated_values ###!!! Can cause duplication of data!
            var.aggregation_step = aggregation_step
            var.harmonization_cache = harmonized_data

        return calculated_values
    
    
    def timeSum(self, var, calculated_values=None, timedata=None, integrate=True): ###!!! Integrate?
        """ Handles the timesum calculation for a variable
            timesums are dependent on the variable aggregation rules, which can be passed directly at definition or are inferred from dependencies
            aggregation rules work with simple seniority: presence of integrate > add > average """
        if calculated_values is None:
            calculated_values = var.values
        if var.aggregation_rule == "sum":
            calculated_values = np.sum(calculated_values)
        elif var.aggregation_rule == "average":
            calculated_values = np.average(calculated_values)
        else:
            raise ValueError(f"Timesum failed: no correct aggregation rule defined for variable {var.name}.")
            
        #Handle integration: if the variable is extensive AND defined as a rate over time, we can multiply by the timestep to obtain a quantity
        if var.is_rate and integrate:
            if timedata is None:
                timestep = var.timestep
            else:
                timestep = timedata[-1]
            calculated_values *= timestep
        
        return calculated_values
    
    
    def executePartialDerivative(self, var, dep_name, absolute_values=False, store_results=True, force_recalculation=False):
        """ Executes the partial derivative executable of a variable for a given dependency, optionally stores values in partial_values dictionary """
        #If no forced recalculation and if the values are already calculated we simply return the already calculated values
        if force_recalculation is False and dep_name in var.partial_values:
            return var.partial_values[dep_name]
        #If there is no executable for this dependency we raise an error
        if var.partial_executables is None:
            if self.equation_engine is not None:
                self.equation_engine.buildPartialDerivativeExecutables(var)
            else:
                raise ValueError(f"Tried to evaluate partial derivative of variable {var.name} while partial derivative executables have not been built yet.")
        
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
        
        #In case of a trivial equation, calculated values will be a Variable object. Here we fix that
        if isinstance(calculated_values, Variable):     ###!!! change this to be handled through an equation engine wrapper
            calculated_values = calculated_values.values
        
        #Take absolute values; i.e. obtain the sensitivities
        if absolute_values:
            calculated_values = np.abs(calculated_values)
        
        #Optionally store results
        if store_results:
            var.partial_values[dep_name] = calculated_values
        return calculated_values
    
    
    
    def executeAllPartials(self, var, absolute_values=False, store_results=True, force_recalculation=False):
        """ Evaluates all partial derivatives of a variable """
        partial_values = {}
        for dep_name in var.dependency_names:
            partial_values[dep_name] = self.executePartialDerivative(var, dep_name, absolute_values=absolute_values, store_results=store_results, force_recalculation=force_recalculation)
        return partial_values
    
    
    



        
        
        


