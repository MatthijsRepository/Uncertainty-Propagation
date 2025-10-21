from my_dataclasses import Variable
import datetime
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
        values = var.executeEquation(store_results=update_var, calculation_engine=self)
        print(f"Calculation of {var.name} complete, values: {var.values}")
        return values
        
    
    def harmonizeTimeSeries(self, dependencies):
        start_times, end_times, timesteps = [], [], []
        for dep in dependencies.values():
            if dep.start_time is not None:
                start_times.append(dep.start_time)
                end_times.append(dep.end_time)
                timesteps.append(dep.timestep)
        print(start_times)
        print(end_times)
        print(timesteps)
        #Determine LCM timestep
        new_timestep = np.lcm.reduce(timesteps)
        
        print(max(start_times))
        print(min(end_times))

        print()
        print(new_timestep)
        
                
    def increaseTemporalResolution(self, var, new_timestep, benchmark_time):
        """ Returns array of values for the new time resolution for a given variable, benchmark time is a specified border between 2 bins """
        #Convert new timestep to datetime object
        new_timestep_datetime = datetime.timedelta(seconds=new_timestep)
        #identify index of the benchmark time ###!!!
        print("WARNING: benchamrk times can be inside a timestep ; i.e. a timestep can span 11:55-12:05 and benchmark time is 12:00 - what then?")
        
        #calculate increase factor
        factor = benchmark_time / var.timestep
        #calculate from where we should start binning:
        
        return
        
    
    def timeSum(self, var):
        """ Handles the timesum calculation for a variable
            timesums are dependent on the variable aggregation rules, which can be passed directly at definition or are inferred from dependencies
            aggregation rules work with simple seniority: presence of integrate > add > average """
        if not var.is_timesum:
            print(f"WARNING: tried to calculate timesum of variable {var.name}, which is not a timesum!")
        
        #Harmonizing dependencies should be done here ###!!!
        print("WARNING: time series harmonization not performed")
        
        args = [var.dependencies[dep_name] for dep_name in var.dependency_names]
        calculated_values = var.executable(*args)
        
        #In case of trivial equations like TS('B'), sympy will not return numeric values but the variable itself, here we manually extract the values if this is the case
        if isinstance(calculated_values, Variable):
            calculated_values = calculated_values.values
        
        #Time aggregation happens according to the variable aggregation rule
        if var.aggregation_rule == "integrate":
            #try to extract a timestep from the dependencies
            for dep in var.dependencies.values():           ###!!! not robust!!!
                if dep.timestep is not None:
                    timestep = dep.timestep
                    break
            calculated_values = np.sum(calculated_values) * timestep
        elif var.aggregation_rule == "sum":
            calculated_values = np.sum(calculated_values)
        else:
            calculated_values = np.average(calculated_values)
        return calculated_values
        
        
        


