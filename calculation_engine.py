from my_dataclasses import Variable, TimeHarmonizationData
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
                    #Calculate variable values, var.values is populated
                    self.executeVariableEquation(var, store_results=True)
        
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
                
        #values = var.executeEquation(store_results=update_var, calculation_engine=self)
        values = self.executeVariableEquation(var, store_results=update_var)
        
        
        if not silent:
            print(f"{indent}Calculation of {var.name} complete, values: {var.values}")
        return values   
    
    def executeVariableEquation(self, var, store_results=True, force_recalculation=True):
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
        
        ###!!!
        #Check if the dependencies are already time-harmonious or time-independent and the equation can be executed directly
        is_harmonious, timedata = self._checkDependencyTimeHarmony(var.dependencies)
        if is_harmonious:
            args = [var.dependencies[dep_name] for dep_name in var.dependency_names]
            calculated_values = var.executable(*args)
        else:
            #Check if the harmonization is already performed
            harmonized_data, timedata = self.harmonizeTimeSeries(var.dependencies, var_name=var.name)
            args = [harmonized_data[dep_name].new_values if dep_name in harmonized_data.keys() else var.dependencies[dep_name] for dep_name in var.dependency_names]
            calculated_values = var.executable(*args)
        
        #Time aggregation if the variable is a timesum
        if var.is_timesum:
            calculated_values = self.timeSum_NEW(var, calculated_values) ###!!!
            timedata = None
        
        #Optionally store the result as the new variable values
        if store_results:
            var.values = calculated_values
            var.setTimeData(timedata)
        return calculated_values
    
    
    def _checkDependencyTimeHarmony(self, dependencies):
        """ Checks if a set of dependencies is time-harmonious """
        start_times, end_times, timesteps = [], [], []
        for dep in dependencies.values():
            if dep.timestep is not None:
                start_times.append(dep.start_time)
                end_times.append(dep.end_time)
                timesteps.append(dep.timestep)
        #If no timestep was retrieved we are dealing exclusively with constants:
        if len(timesteps)==0:
            return True, None
        #If all start times, end times and timesteps are the same, the time series are harmonious
        elif len(set(start_times))==1 and len(set(end_times))==1 and len(set(timesteps))==1:
            return True, (start_times[0], end_times[0], timesteps[0])
        #Otherwise the dpeendencies are not time-harmonious
        else:
            return False, None
                    
    def _pruneHarmonizedTimeSeriesTails(self, harmonized_dataset, new_timestep):
        """ Takes a dictionary of harmonized time series data and ensures all start times and end times match perfectly,
            prunes data from the harmonized dataset that do not fall within the common timerange. 
            Returns updated dataset and a tuple of the common start time, end time and timestep """
        datetime_timestep = timedelta(seconds=new_timestep)
        
        #Extract common start, end time
        start_times, end_times = [], []
        for harmonized_data in harmonized_dataset.values():
            start_times.append(harmonized_data.new_start_time)
            end_times.append(harmonized_data.new_end_time)
        common_start_time = max(start_times)
        common_end_time   = min(end_times)
        
        for harmonized_data in harmonized_dataset.values():
            #Prune start
            offset_steps = int((common_start_time - harmonized_data.new_start_time) / datetime_timestep)
            if offset_steps>0:
                harmonized_data.new_values              = harmonized_data.new_values[offset_steps:]
                harmonized_data.new_start_time          += offset_steps * datetime_timestep
                harmonized_data.prune_offset_start      = offset_steps
            #Prune tail
            offset_steps_end = int((harmonized_data.new_end_time - common_end_time) / datetime_timestep)
            if offset_steps_end>0:
                harmonized_data.new_values              = harmonized_data.new_values[:-offset_steps_end]
                harmonized_data.new_end_time            -= offset_steps_end * datetime_timestep
                harmonized_data.prune_offset_end        = offset_steps_end
            #Update the harmonized_data dictionary with pruned datasets
        return harmonized_dataset, (common_start_time, common_end_time, new_timestep)
    
    def harmonizeTimeSeries(self, dependencies, var_name=None, smuggle_limit=0):
        dep_names, start_times, timesteps = [], [], []
        for dep in dependencies.values():
            if dep.timestep is not None:
                dep_names.append(dep.name)
                start_times.append(dep.start_time)
                timesteps.append(dep.timestep)
        
        #Determine LCM timestep
        new_timestep = float(np.lcm.reduce(timesteps))

        #We take as benchmark time the start time of the dataset with the biggest timestep
        benchmark_time = start_times[np.argmax(timesteps)]
        print(f"Warning: benchmark time for harmonization is start time of variable with biggest timestep: {benchmark_time}")
        #benchmark_time = common_start_time
        #print(f"Warning: benchmark time is common start time: {benchmark_time}")
        ###!!!
        
        harmonized_data = {}
        for dep_name in dep_names:
            temp_harmonization_data = self.decreaseTemporalResolution(dependencies[dep_name], new_timestep, benchmark_time=benchmark_time, smuggle_limit=smuggle_limit)
            temp_harmonization_data.target_var_name = var_name
            harmonized_data[dep_name] = temp_harmonization_data
            #each entry is a TimeHarmonizationData object
        harmonized_data, new_timedata = self._pruneHarmonizedTimeSeriesTails(harmonized_data, new_timestep)
        return harmonized_data, new_timedata
    
    def _rebinTimeSeries(self, var, low_index, high_index, low_fraction, high_fraction, factor):
        """ Handles rebinning of variable time series data into a new timeseries of greater granularity 
            Allows for fractional splitting of old bins between two new bins """
        #Note cannot handle full time aggregation
        new_values = np.zeros( int((high_index-low_index)/factor) )
        if len(new_values)==1:
            raise ValueError("Aggregating time series data into a single bin currently not supported. - Timesum calling not implemented yet")

        #Note: new bins can contain fractions of old bins - thus our loop should start on the same bin as the previous iteration ended on
        #Since new bin i can for instance include 2/3 of old bin j, then new bin i+1 should contain 1/3 of old bin j!
        for i in range(len(new_values)):
            start = low_index + i*factor
            if start < 0:
                #First bin handling
                new_values[i] = np.sum(var.values[:start+factor])
                new_values[i] += abs(low_index)-1 * var.values[0]
                new_values[i] += var.values[0] * low_fraction
                new_values[i] += var.values[start+factor] * high_fraction
            elif start + factor-1 >= len(var.values):
                #Last bin handling
                new_values[i] = np.sum(var.values[start:])
                new_values[i] += var.values[start-1] * low_fraction
                new_values[i] += (high_index - len(var.values)-1) * var.values[-1]
                new_values[i] += var.values[-1] * high_fraction
            else:
                #General bin handling
                new_values[i] = np.sum(var.values[start:start+factor-1])
                new_values[i] += var.values[start] * low_fraction
                new_values[i] += var.values[start+factor-1] * high_fraction

        
        #If aggregation rule is average, we take the time average
        if var.aggregation_rule == "average":
            new_values /= factor
        return new_values
        
    def decreaseTemporalResolution(self, var, new_timestep, benchmark_time, smuggle_limit=0):
        """ Returns array of values for the new time resolution for a given variable, benchmark time is a specified border between 2 bins """
        #Check if temporal operations make sense for this variable (e.g. not a float)
        if isinstance(var.values, float):
            raise ValueError(f"Rebinning of variable {var.name} terminated, variable is a constant.")
        
        #We calculate the new start and end times by seeing where bin limits - given the benchmark - fit in the previous timerange
        #Note: we allow to smuggle a specified amount; to allow manual avoiding of instances where an hour of data is discarded based on a small time mismatch.
        #We calculate the difference between the benchmark time and the start time (negative if benchmark time > start time)
        delta_start = (var.start_time - benchmark_time).total_seconds()
        #We calculate how many new timesteps can be taken down from the benchmark time until we are at the supremum step under the actual start time
        offset_steps = np.floor(delta_start/new_timestep)
        #Declare new start time as this time
        new_start_time = benchmark_time + timedelta(seconds=offset_steps * new_timestep)
        #Check whether this start time (at or under the actual start time) falls within the smuggle limit, if not, we increase by 1 timestep so we are within the allowed limit
        if new_start_time<var.start_time:
            if (var.start_time - new_start_time).total_seconds() > smuggle_limit:
                new_start_time += timedelta(seconds=new_timestep)
                
        #Same procedure but opposite logic
        delta_end = (var.end_time - benchmark_time).total_seconds()
        offset_steps_end = np.ceil(delta_end / new_timestep)
        new_end_time = benchmark_time + timedelta(seconds=offset_steps_end * new_timestep)
        if new_end_time>var.end_time:
            if (new_end_time - var.end_time).total_seconds() > smuggle_limit:
                new_end_time -= timedelta(seconds=new_timestep)
        
        #Identify the factor increase
        factor = new_timestep / var.timestep
        if factor.is_integer():
            factor = int(factor)
        else:
            raise ValueError(f"Temporal granularity increase for variable {var.name} failed: can only be performed for integer multiples of the old timestep! Attempted increase factor was {factor}.")
        
        #Identify first index to include in first bin
        low_index = int((new_start_time - var.start_time).total_seconds() // var.timestep)
        low_fraction = 1 - ((new_start_time - var.start_time).total_seconds() / var.timestep - low_index)
        #Identify last index to include in last bin
        high_index = int((new_end_time - var.end_time).total_seconds() // var.timestep) + len(var.values)
        high_fraction = 1-low_fraction  ### WARNING: we assume new timestep is always an integer multiple of the old timestep. If this is not the case, this method does not work

        new_values = self._rebinTimeSeries(var, low_index, high_index, low_fraction, high_fraction, factor)
        
        harmonization_data = TimeHarmonizationData(
            dep_var_name    = var.name,
            base_timestep   = var.timestep,
            new_timestep    = var.timestep * factor,
            new_start_time  = new_start_time,
            new_end_time    = new_end_time,
            low_fraction    = low_fraction,
            high_fraction   = high_fraction,
            upsample_factor = factor,
            new_values      = new_values)
        
        """
        print()
        print(f"Variable: {var.name}")
        print(new_values)
        print(f"Low index, high index: {low_index}, {high_index}")
        print(f"Offset steps start: {offset_steps}, end: {offset_steps_end}")
        #print(offset_steps)
        #print(offset_steps_end)
        print("Old start, end")
        print(var.start_time)
        print(var.end_time)
        print("New start, end")
        print(new_start_time)
        print(new_end_time)
        """
        return harmonization_data
        
    def timeSum_NEW(self, var, calculated_values):
        """ Handles the timesum calculation for a variable
            timesums are dependent on the variable aggregation rules, which can be passed directly at definition or are inferred from dependencies
            aggregation rules work with simple seniority: presence of integrate > add > average """
        #Time aggregation happens according to the variable aggregation rule
        if var.aggregation_rule == "integrate":
            #try to extract a timestep from the dependencies
            for dep in var.dependencies.values():           ###!!! not robust!!!
                if dep.timestep is not None:
                    timestep = dep.timestep   ###!!! WRONG
                    break
            calculated_values = np.sum(calculated_values) * timestep
        elif var.aggregation_rule == "sum":
            calculated_values = np.sum(calculated_values)
        else:
            calculated_values = np.average(calculated_values)
        return calculated_values
    
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
        
        
        


