from my_dataclasses import Variable, TimeHarmonizationData
import datetime
import numpy as np
from datetime import datetime, timedelta


    
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
    
    def executeVariableEquation(self, var, store_results=True, force_recalculation=False):
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
            harmonized_data = None
            args = [var.dependencies[dep_name] for dep_name in var.dependency_names]
            calculated_values = var.executable(*args)
        else:
            #Check if the harmonization is already performed
            if var.harmonization_cache is not None and force_recalculation is False:
                harmonized_data = var.harmonization_cache
                timedata = var.getTimeData()
            #Else: perform time harmonization
            else:
                harmonized_data, timedata = self.harmonizeTimeSeries(var.dependencies, var_name=var.name)
            args = [harmonized_data[dep_name].new_values if dep_name in harmonized_data.keys() else var.dependencies[dep_name] for dep_name in var.dependency_names]
            calculated_values = var.executable(*args)
        
        #Time aggregation if the variable is a timesum
        aggregation_step      = None #Aggregation step is the timestep of the timesummed data. This value is required for uncertainty calculation and therefore passed to the variable
        non_aggregated_values = None #Calculated values before aggregation - required for uncertainty calculation
        if var.is_timesum:
            non_aggregated_values = calculated_values
            calculated_values     = self.timeSum(var, calculated_values, timedata) ###!!!
            aggregation_step      = timedata[2]
            timedata              = None
        
        #Optionally store the result as the new variable values
        if store_results:
            var.values = calculated_values
            var.setTimeData(timedata)
            var.non_aggregated_values = non_aggregated_values ###!!! Can cause duplication of data!
            var.aggregation_step = aggregation_step
            var.harmonization_cache = harmonized_data

        return calculated_values
    
    def timeSum(self, var, calculated_values=None, timedata=None):
        """ Handles the timesum calculation for a variable
            timesums are dependent on the variable aggregation rules, which can be passed directly at definition or are inferred from dependencies
            aggregation rules work with simple seniority: presence of integrate > add > average """
        if calculated_values is None:
            calculated_values = var.values
        #Time aggregation happens according to the variable aggregation rule
        if var.aggregation_rule == "integrate":
            if timedata is None:
                raise ValueError(f"Time integration of variable {var.name} = {var.equation} failed - no timedata was provided and hence no timestep could be extracted.")
            elif isinstance(timedata, float):   #In this case an actual timestep is passed
                calculated_values = np.sum(calculated_values) * timedata
            else:                               #In this case a timedata tuple is passed - timestep is always on index 2
                calculated_values = np.sum(calculated_values) * timedata[2]
        elif var.aggregation_rule == "sum":
            calculated_values = np.sum(calculated_values)
        elif var.aggregation_rule == "average":
            calculated_values = np.average(calculated_values)
        else:
            raise ValueError(f"Timesum failed: no aggregation rule defined for variable {var.name}.")
        return calculated_values
    
    def harmonizeTimeSeries(self, dependencies, var_name=None, smuggle_limit=0):
        """ Given a dictionary of variables, harmonizes the variables to lcm timestep and identical start and end times.
            Returns a dictionary of rebinned values and the various settings, offsets and fractions calculated for rebinning.
            Also returns a tuple containing the new start time, end time and timestep of the harmonized data. """
        dep_names, start_times, timesteps = [], [], []
        for dep in dependencies.values():
            if dep.timestep is not None:
                dep_names.append(dep.name)
                start_times.append(dep.start_time)
                timesteps.append(dep.timestep)
        
        #Determine LCM timestep
        new_timestep = float(np.lcm.reduce(timesteps))

        #We take as benchmark time the start time of the dataset with the biggest timestep ###!!!
        benchmark_time = start_times[np.argmax(timesteps)]
        print(f"Warning: benchmark time for harmonization is start time of variable with biggest timestep: {benchmark_time}")
    
        #Populate the harmonized data dictionary with a TimeHarmonizationData object for each variable in the given set.        
        harmonized_data = {}
        for dep_name in dep_names:
            temp_harmonization_data = self.calculateTimeHarmonizationData(dependencies[dep_name], new_timestep, \
                                                                          benchmark_time=benchmark_time, smuggle_limit=smuggle_limit)
            temp_harmonization_data.new_values = self._rebinTimeSeries(dependencies[dep_name], temp_harmonization_data)
            temp_harmonization_data.target_var_name = var_name
            harmonized_data[dep_name] = temp_harmonization_data
            
        #Prune the datasets such that they all have the same start and end times, necessary for computations to make sense.
        harmonized_data, new_timedata = self._pruneHarmonizedTimeSeriesTails(harmonized_data, new_timestep)
        return harmonized_data, new_timedata
    
    def calculateTimeHarmonizationData(self, var, new_timestep, benchmark_time=None, smuggle_limit=0):
        """ Returns a populated TimeHarmonizationData object, which stores rebinning information such as base timestep, new timestep, 
            new start and end times, bin-offset fraction with respect to original bin borders and other metadata
            Allows for smuggling a bit with the bins: suppose the dataset ends at 23:59:30 and the envisioned new dataset would end at 00:00:00, 
            ...then it allows to extend the original dataset a bit to accomodate this last bin """
        #Check if temporal operations make sense for this variable (e.g. not a float)
        if isinstance(var.values, float):
            raise ValueError(f"Rebinning of variable {var.name} terminated, variable is a constant.")
        #Default benchmark time to 12:00:00
        if benchmark_time is None:
            benchmark_time = datetime.strptime("12:00:00", "%H:%M:%S")
        
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

        return TimeHarmonizationData(
            dep_var_name    = var.name,
            base_timestep   = var.timestep,
            new_timestep    = var.timestep * factor,
            new_start_time  = new_start_time,
            new_end_time    = new_end_time,
            low_index       = low_index,
            high_index      = high_index,
            low_fraction    = low_fraction,
            high_fraction   = high_fraction,
            upsample_factor = factor)
    
    def decreaseVariableTemporalResolution(self, var, new_timestep, benchmark_time=None, update_var=True, smuggle_limit=0):
        """ Can be used to calculate full TimeHarmonizationData object for a variable - and optionally immediately updates this variable """
        harmonization_data = self.calculateTimeHarmonizationData(var, new_timestep, benchmark_time=benchmark_time, smuggle_limit=smuggle_limit)
        harmonization_data.new_values = self._rebinTimeSeries(var, harmonization_data)

        #Update the values of the actual variable
        if update_var:
            var.values = harmonization_data.new_values
            var.setTimeData((harmonization_data.new_start_time, harmonization_data.new_end_time, harmonization_data.new_timestep))
            var.uncertainty.reset()
        return harmonization_data
    
    def _rebinTimeSeries(self, var, harmonization_data):
        """ Handles rebinning of variable time series data into a new timeseries of greater granularity 
            Allows for fractional splitting of old bins between two new bins """
        low_index, high_index       = harmonization_data.low_index, harmonization_data.high_index
        low_fraction, high_fraction = harmonization_data.low_fraction, harmonization_data.high_fraction
        factor                      = harmonization_data.upsample_factor
        
        #Note: function cannot handle full time aggregation - use timesum for that
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

        #If aggregation rule is to average or integrate, we take the time average
        if var.aggregation_rule == "average" or var.aggregation_rule == "integrate":
            new_values /= factor
        return new_values
    
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
        



        
        
        


