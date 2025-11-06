import numpy as np
from datetime import datetime, timedelta
from my_dataclasses import Variable, TimeHarmonizationData


class TimeEngine:
    def __init__(self):
        return
    
    def ensureDependencyTimeHarmony(self, var, force_recalculation=False):
        """ Ensures all dependencies of a variable are time-harmonious, returns the arguments to be passed in the executable, a new tuple of timedata 
            and a harmonized data dictionary of TimeHarmonizationData objects containing rebinning information for each dependency """
        #Check if the dependencies are already time-harmonious or time-independent and the equation can be executed directly
        is_harmonious, timedata = self._checkDependencyTimeHarmony(var.dependencies)
        if is_harmonious:
            harmonized_data = None
            args = [var.dependencies[dep_name] for dep_name in var.dependency_names]
        else:
            #Check if the harmonization is already performed, cached and recalculation not required
            if var.harmonization_cache is not None and force_recalculation is False:
                harmonized_data = var.harmonization_cache
                timedata = var.getTimeData()
            #Else: perform time harmonization
            else:
                harmonized_data, timedata = self.harmonizeTimeSeries(var.dependencies, var_name=var.name)
            args = [harmonized_data[dep_name].new_values if dep_name in harmonized_data.keys() else var.dependencies[dep_name] for dep_name in var.dependency_names]
        return args, timedata, harmonized_data
    
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
        print(f"Warning: benchmark time for harmonization is start time of variable with biggest timestep: {benchmark_time.strftime('%H:%M:%S')}")
    
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
            then it allows to extend the original dataset a bit to accomodate this last bin """
        #Check if temporal operations make sense for this variable (e.g. not a float)
        if isinstance(var.values, float):
            raise ValueError(f"Rebinning of variable {var.name} terminated, variable is a constant.")
        #Default benchmark time to 12:00:00
        if benchmark_time is None:
            benchmark_time = datetime.strptime("12:00:00", "%H:%M:%S")
        if isinstance(benchmark_time, str):
            benchmark_time = datetime.strptime(benchmark_time, "%H:%M:%S")
        
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
        delta_end = (var.last_time - benchmark_time).total_seconds()
        offset_steps_end = np.ceil(delta_end / new_timestep)
        new_last_time = benchmark_time + timedelta(seconds=offset_steps_end * new_timestep)
        if new_last_time>var.last_time:
            if (new_last_time - var.last_time).total_seconds() > smuggle_limit:
                new_last_time -= timedelta(seconds=new_timestep)
        
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
        high_index = int((new_last_time - var.last_time).total_seconds() // var.timestep) + len(var.values)
        high_fraction = 1-low_fraction  ### WARNING: we assume new timestep is always an integer multiple of the old timestep. If this is not the case, this method does not work

        return TimeHarmonizationData(
            dep_var_name    = var.name,
            base_timestep   = var.timestep,
            new_timestep    = var.timestep * factor,
            new_start_time  = new_start_time,
            new_last_time   = new_last_time,
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
            var.setTimeData((harmonization_data.new_start_time, harmonization_data.new_last_time, harmonization_data.new_timestep))
            var.uncertainty.reset()
            var.uncertainty.rescaleUncertaintySources(harmonization_data.upsample_factor)
            print(f"WARNING: hard-rescaled variable {var.name} to a new timestep. Direct uncertainties for this variable will be rescaled, but this action is destructive for uncertainty information.")            
        return harmonization_data
    
    def _rebinTimeSeries(self, var, harmonization_data):
        """ Handles rebinning of variable time series data into a new timeseries of greater granularity 
            Allows for fractional splitting of old bins between two new bins 
            NOTE: this function is almost exactly copied by the uncertainty engine for correlation aggregation """
        low_index, high_index       = harmonization_data.low_index, harmonization_data.high_index
        low_fraction, high_fraction = harmonization_data.low_fraction, harmonization_data.high_fraction
        factor                      = harmonization_data.upsample_factor
        
        new_values = np.zeros( int((high_index-low_index)/factor) )
        if len(new_values)==1:
            raise ValueError("Aggregating time series data into a single bin currently not supported. - Timesum calling not implemented yet") ###!!!

        for i in range(len(new_values)):
            start = low_index + i*factor
            if start < 0:
                #First bin handling
                new_values[i] = np.sum(var.values[:start+factor])
                new_values[i] += (abs(start)-1) * var.values[0]
                new_values[i] += var.values[0] * low_fraction
                new_values[i] += var.values[start+factor] * high_fraction
            elif start + factor >= len(var.values):
                #Last bin handling
                new_values[i] = np.sum(var.values[start+1:])
                new_values[i] += var.values[start] * low_fraction
                new_values[i] += (high_index - len(var.values)) * var.values[-1]
                new_values[i] += var.values[-1] * high_fraction
            else:
                #General bin handling
                new_values[i] = np.sum(var.values[start+1:start+factor])
                new_values[i] += var.values[start] * low_fraction
                new_values[i] += var.values[start+factor] * high_fraction

        #If aggregation rule is to average, we take the time average
        if var.aggregation_rule == "average": 
            new_values /= factor
        return new_values
    
    def _checkDependencyTimeHarmony(self, dependencies):
        """ Checks if a set of dependencies is time-harmonious """
        start_times, last_times, timesteps = [], [], []
        for dep in dependencies.values():
            if dep.timestep is not None:
                start_times.append(dep.start_time)
                last_times.append(dep.last_time)
                timesteps.append(dep.timestep)
        #If no timestep was retrieved we are dealing exclusively with constants:
        if len(timesteps)==0:
            return True, None
        #If all start times, end times and timesteps are the same, the time series are harmonious
        elif len(set(start_times))==1 and len(set(last_times))==1 and len(set(timesteps))==1:
            return True, (start_times[0], last_times[0], timesteps[0])
        #Otherwise the dependencies are not time-harmonious
        else:
            return False, None
                    
    def _pruneHarmonizedTimeSeriesTails(self, harmonized_dataset, new_timestep):
        """ Takes a dictionary of harmonized time series data and ensures all start times and end times match perfectly,
            prunes data from the harmonized dataset that do not fall within the common timerange. 
            Returns updated dataset and a tuple of the common start time, end time and timestep """
        datetime_timestep = timedelta(seconds=new_timestep)
        
        #Extract common start, end time
        start_times, last_times = [], []
        for harmonized_data in harmonized_dataset.values():
            start_times.append(harmonized_data.new_start_time)
            last_times.append(harmonized_data.new_last_time)
        common_start_time = max(start_times)
        common_last_time  = min(last_times)
        
        for harmonized_data in harmonized_dataset.values():
            #Prune start
            offset_steps = int((common_start_time - harmonized_data.new_start_time) / datetime_timestep)
            if offset_steps>0:
                harmonized_data.new_values              = harmonized_data.new_values[offset_steps:]
                harmonized_data.new_start_time          += offset_steps * datetime_timestep
                harmonized_data.prune_offset_start      = offset_steps
            #Prune tail
            offset_steps_end = int((harmonized_data.new_last_time - common_last_time) / datetime_timestep)
            if offset_steps_end>0:
                harmonized_data.new_values              = harmonized_data.new_values[:-offset_steps_end]
                harmonized_data.new_last_time           -= offset_steps_end * datetime_timestep
                harmonized_data.prune_offset_end        = offset_steps_end
            #Update the harmonized_data dictionary with pruned datasets
        return harmonized_dataset, (common_start_time, common_last_time, new_timestep)
        



