from my_dataclasses import Variable
import datetime
import numpy as np
from datetime import timedelta


    
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
        
    
    def calculateValues(self, var, update_var=True, calculate_dependencies=True, silent=True, indent=""):
        """ Checks if the values of all dependencies are calculated, optionally calculates dependencies, and then calculates values of specified variable """
        if not silent:
            print(f"{indent}Calculating values of variable {var.name} with dependencies {var.dependency_names}")
        for dep in var.dependencies.values():
            if dep.values is not None:
                continue
            elif calculate_dependencies is False:
                raise ValueError(f"Calculation of variable {var.name} failed: dependency {dep.name} has no values defined and automatic dependency calculation is turned off.")
            else:
                self.calculateValues(dep, update_var=update_var, silent=silent, indent=f"   {indent}")
        values = var.executeEquation(store_results=update_var, calculation_engine=self)
        if not silent:
            print(f"{indent}Calculation of {var.name} complete, values: {var.values}")
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
        
                
    def increaseTemporalResolution(self, var, new_timestep, benchmark_time, smuggle_limit=0):
        """ Returns array of values for the new time resolution for a given variable, benchmark time is a specified border between 2 bins """
        #Check if temporal operations make sense for this variable (e.g. not a float)
        if isinstance(var.values, float):
            raise ValueError(f"Rebinning of variable {var.name} terminated, variable is a constant.")
        
        #Convert new timestep to datetime object
        new_timestep_datetime = datetime.timedelta(seconds=new_timestep)
        #identify index of the benchmark time ###!!!
        print("WARNING: benchmark times can be inside a timestep ; i.e. a timestep can span 11:55-12:05 and benchmark time is 12:00 - handle correct time split/aggregation here")
        
        #calculate increase factor
        factor = new_timestep / var.timestep
        print(var.start_time)
        print(var.first_time)
        #calculate from where we should start binning:
        present_bin_edges = np.linspace(var.start_time, var.end_time, len(var.values)+1)
        print(present_bin_edges)
        
        #We calculate the new start and end times by seeing where bin limits - given the benchmark - fit in the previous timerange
        #Note: we allow to smuggle a a specified amount; to allow manual avoiding of instances where an hour of data is discarded based on a small time mismatch.
        delta_start = (var.start_time - benchmark_time).total_seconds()
        offset_steps = np.floor(delta_start/new_timestep)
        first_edge = benchmark_time + timedelta(seconds=offset_steps * new_timestep)
        if first_edge<var.start_time:
            if (var.start_time - first_edge).total_seconds() > smuggle_limit:
                first_edge += timedelta(seconds=new_timestep)
                
        
        delta_end = (var.end_time - benchmark_time).total_seconds()
        offset_steps_end = np.ceil(delta_end / new_timestep)
        last_edge = benchmark_time + timedelta(seconds=offset_steps_end * new_timestep)
        if last_edge>var.end_time:
            if (last_edge - var.end_time).total_seconds() > smuggle_limit:
                last_edge -= timedelta(seconds=new_timestep)
        
        
        print(var.start_time)
        print(var.end_time)
        print(first_edge)
        print(last_edge)
        #factor = benchmark_time / var.timestep
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
        
        
        


