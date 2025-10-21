from my_dataclasses import Variable
import numpy as np


def TIMESUM_TEMP(var, options=None):    ###!!!
    print("TIMESUM_TEMP test")
    print(var)
    if isinstance(var, (float,int)):
        return var
    elif isinstance(var.values, (float, int)):
        raise ValueError(f"Timesums of constants are not supported yet (variable {var.name})")
    else:
        return var.TimeSum()
    
    
    
class CalculationEngine:
    def __init__(self, variables):
        self.variables = variables
        return
    
    def add(self, var1, var2):
        return var1 + var2
    
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
                    #Calculate variable values, var.values is populated by default
                    var.executeEquation()
        
    
    def calculateValues(self, var, update_var=True, calculate_dependencies=True):
        """ Checks if the values of all dependencies are calculated, optionally calculates dependencies, and then calculates values of specified variable """
        
        print(f"\nCalculating values of variable {var.name} with dependencies {var.dependency_names}")
        for dep in var.dependencies.values():
            if dep.values is not None:
                continue
            elif calculate_dependencies is False:
                raise ValueError(f"Calculation of variable {var.name} failed: dependency {dep.name} has no values defined and automatic dependency calculation is turned off.")
            else:
                print(f"Dependency detected with no values: {dep.name}. Calculating...")
                self.calculateValues(dep, update_var=update_var)
        
        #Check if the variable is a timesum: in this case we have to perform different actions
        #if var.is_timesum:
        #    values = var.executeEquation(store_results=False)
        #    values = np.sum(values) * var.timestep
        #    if update_var:
        #        var.values = values
        #else:
        values = var.executeEquation(store_results=update_var, calculation_engine=self)
        print(f"Calculation of {var.name} complete, values: {np.sum(var.values)}")
        return values
        
    
    def harmonizeTimeSeries(self, dependencies):
        start_times, end_times, timesteps = []
        for dep in dependencies.values():
            if dep.start_time is not None:
                start_times.append(dep.start_time)
                end_times.append(dep.end_time)
                timesteps.append(dep.timestep)
        
                
    
    
    def timeSum(self, var):
        if not var.is_timesum:
            print(f"WARNING: tried to calculate timesum of variable {var.name}, which is not a timesum!")
        
        
        #Harmonizing dependencies should be done here ###!!!
        args = [var.dependencies[dep_name] for dep_name in var.dependency_names]
        calculated_values = var.executable(*args)
        print("WARNING: timesum calculation placeholder!")
        return 10
        
        
        


