import numpy as np
from datetime import datetime, timedelta
from my_dataclasses import Variable, TimeHarmonizationData


class TimeEngine:
    """
    The TimeEngine is an object of helper functions, that may be invoked during calculations to ensure proper time-matching of 
    timeseries of different ranges and resolutions. 
    
    The time engine generally (aside from the `decreaseVariableTemporalResolution` function) does not change a variable's own data, 
    but creates a `TimeHarmonizationData` object containing all information to bring variables to matching resolutions, 
    and potentially has already calculated the new timeseries and included these in this object.
    Calculation and uncertainty engines use the outputs of this engine as inputs for their own calculations and methods.
    """
    def __init__(self):
        return
    
    def ensureDependencyTimeHarmony(self, var, force_recalculation=False):
        """
        For a given variable, ensures the dependencies of a variable are time-harmonious. 
        
        Invoked by the calculation engine to obtain the inputs for a variable's executable upon calculation of a variable.
        Orchestrator function that checks for cached harmonizations and invokes calculations.
        Returns the arguments to be passed into the variable's executable, the timedata (start time, end time, timestep) 
        of the resulting timeseries, and the accompanying `TimeHarmonizationData` metadata objects, which contain
        rebinning information for each dependency.
        
        Parameters
        ----------
        var: Variable
            Variable for which to ensure dependency time harmony.
        force_recalculation: bool, default=False
            Boolean indicating whether time harmonization should be recalculated, in case a cached harmonization is detected.
        Returns
        -------
        tuple[ list[np.array or numeric], tuple[timestamp, timestamp, int], dict[str, TimeHarmonizationData] ]
            Tuple containing the following:
            - A list of timeseries and scalars, which are the numeric data to be given to `var.executable` to calculate `var`'s values.
            - A tuple containing the start time, end time and timestep (in seconds) of the resulting timeseries of `var`. None if it is dependent on constants.
            - A dictionary of dependency names and `TimeHarmonizationData` objects, containing rebinning information.
        """
        #Check if the dependencies are already time-harmonious or time-independent and the equation can be executed directly
        is_harmonious, timedata = self._checkDependencyTimeHarmony(var.dependencies)
        
        #Short-circuit if dependencies are time-harmonious
        if is_harmonious:
            harmonized_data = None
            args = [var.dependencies[dep_name] for dep_name in var.dependency_names]
            return args, timedata, harmonized_data
        
        #Check if the harmonization is already performed, cached and recalculation not required
        if var.harmonization_cache is not None and not force_recalculation:
            harmonized_data = var.harmonization_cache
            timedata = var.getTimeData()
        #Else: perform time harmonization
        else:
            harmonized_data, timedata = self.harmonizeTimeSeries(var.dependencies, var_name=var.name)
        
        args = [harmonized_data[dep_name].new_values 
                if dep_name in harmonized_data.keys() 
                else var.dependencies[dep_name] 
                for dep_name in var.dependency_names]
        
        return args, timedata, harmonized_data
    
    def harmonizeTimeSeries(self, dependencies, var_name=None, benchmark_time=None, smuggle_limit=0):
        """ 
        Main timeseries harmonization handler. For a given set of dependencies, calculates common start, end times,
        the lcm of their timesteps, and the accompanying `TimeHarmonizationData` objects containing rebinning information.
        
        Parameters
        ----------
        dependencies: dict[str, Variable]
            Dictionary of dependency names and the corresponding variables which must be time-harmonized.
        var_name: str or None, default=None
            Name of the calling variable to give to harmonization cache, optional.
        benchmark_time: datetime.datetime or None, default=None
            Timestamp along which the boundary between two bins must be, thereby defining all bins (as timestep is inferred).
            If `None`, this time is set as the start time of the variable with the greatest timestep.
        smuggle_limit: int or float, default=0
            Amount of seconds that is allowed to be smuggled with timebins during harmonization, if two timeseries mismatch by a small amount of time.
            If one timeseries is slightly shorter than the other timeseries, data is duplicated to stretch the timespan over the missing timerange.
            Useful to prevent data with large timesteps to be excluded because data with a short timestep ends a little too soon.
            Recommended to duplicate timeseries manually before execution in such cases, instead of using this smuggling option.
        
        Returns
        -------
        tuple[ dict[str, TimeHarmonizationData], tuple[timestamp, timestamp, int] ]
            Returns a tuple containing
            - A dictionary of `TimeHarmonizationData` objects for each dependency, containing rebinning information.
            - A tuple containing common start and end times, and lcm timestep in seconds.
        """
        dep_names, start_times, timesteps = [], [], []
        for dep in dependencies.values():
            if dep.timestep is not None:
                dep_names.append(dep.name)
                start_times.append(dep.start_time)
                timesteps.append(dep.timestep)
        
        #Determine LCM timestep
        new_timestep = float(np.lcm.reduce(timesteps))

        #We take as benchmark time the start time of the dataset with the biggest timestep
        if benchmark_time is None:
            benchmark_time = start_times[np.argmax(timesteps)]
            #print(f"Warning: benchmark time for harmonization is start time of variable with biggest timestep: {benchmark_time.strftime('%H:%M:%S')}")
    
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
        """ 
        For a given variable, new timestep and benchmark time, builds the corresponding `TimeHarmonizationData` object,
        which contains rebinning information such as timestep upscale factor, number of timesteps to prune, new start and end times.
        
        Allows for smuggling a bit with the timebins, see `TimeEngine.harmonizeTimeSeries`. Allows for smuggling a bit with the bins.
        Suppose the timeseries ends at 23:59:30 and the envisioned new timeseries would end at 00:00:00, then it allows to extend 
        the original timeseries a bit to accomodate this last bin, if the difference in seconds is within the smuggling limit.
        
        Parameters
        ----------
        var: Variable
            Variable for which the harmonization data must be calculated.
        new_timestep: int
            New timestep in seconds.
        benchmark_time: datetime.datetime or None, default=None
            Timestamp along which the boundary between two bins must be, thereby defining all bins (as the timestep is given).
            If `None`, takes 12:00:00.
         smuggle_limit: int or float, default=0
             Amount of seconds that is allowed to be smuggled with timebins during harmonization.
        
        Returns
        -------
        TimeHarmonizationData
            Populated `TimeHarmonizationData` object, containing all information to rebin `var` to the desired temporal resolution and range.
        
        Raises
        ------
        ValueError
            If `var.values` is a scalar, since this function only works on timeseries data.
        ValueError
            If the desired new timestep is not an integer multiple of the variable's existing timestep.
            
        Returns a populated TimeHarmonizationData object, which stores rebinning information such as base timestep, new timestep, 
        new start and end times, bin-offset fraction with respect to original bin borders and other metadata
        Allows for smuggling a bit with the bins: suppose the dataset ends at 23:59:30 and the envisioned new dataset would end at 00:00:00, 
        then it allows to extend the original dataset a bit to accomodate this last bin 
        """
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
    
    def decreaseVariableTemporalResolution(self, var, new_timestep, benchmark_time=None, update_var=False, new_var=False, smuggle_limit=0):
        """ 
        Function that can be used to irreversibly decrease a variable's temporal resolution.
        Either updates the variable irreversibly, or returns a new `Variable` instance containing the new temporal resolution.
        
        Parameters
        ----------
        var: Variable 
            Variablefor which temporal resolution decrease is desired.
        new_timestep: int
            New timestep in seconds.
        benchmark_time: datetime.datetime or None, default=None
            Timestamp along which the boundary between two bins must be, thereby defining all bins (as timestep is inferred).
            If `None`, takes 12:00:00.
        update_var: bool, default=False
            Whether the original variable should be irreversibly updated. Mutually exclusive with `new_var`.
        new_var: bool, default=False
            Whether a new variable should be created containing the new data. Mutually exclusive with `update_var`.
        smuggle_limit: int or float, default=0
            Amount of seconds that is allowed to be smuggled with timebins during harmonization.
        
        Returns
        -------
        Variable or TimeHarmonizationData
            Returns a `Variable` instance if the option `new_var` is chosen.
            Returns `TimeHarmonizationData` if `var` is chosen to be updated.
        
        Raises
        ------
        ValueError
            If both `new_var` and `update_var` have been set as `True`.
        """
        harmonization_data = self.calculateTimeHarmonizationData(var, new_timestep, benchmark_time=benchmark_time, smuggle_limit=smuggle_limit)
        harmonization_data.new_values = self._rebinTimeSeries(var, harmonization_data)
        
        if new_var and update_var:
            raise ValueError("Tried to decrease temporal resolution for a variable while simultaneously trying to place it inside a new variable. Please choose one of these options.")
        if new_var:
            new_variable = Variable(name=f"upsampled_{var.name}", values = harmonization_data.new_values, is_basic=False, equation=f"'{var.name}'",
                                    is_rate=var.is_rate, aggregation_rule=var.aggregation_rule, first_time=harmonization_data.getFirstTime(), last_time=harmonization_data.new_last_time)
            new_variable.dependency_names = [var.name]
            new_variable.dependencies = {var.name : var}
            new_variable.harmonization_cache = {var.name : harmonization_data}
            return new_variable
            
        #Update the values of the actual variable
        if update_var:
            var.values = harmonization_data.new_values
            var.setTimeData((harmonization_data.new_start_time, harmonization_data.new_last_time, harmonization_data.new_timestep))
            var.uncertainty.reset()
            var.uncertainty.rescaleUncertaintySources(harmonization_data.upsample_factor)
            print(f"WARNING: hard-rescaled variable {var.name} to a new timestep. Direct uncertainties for this variable will be rescaled, but this action is destructive for uncertainty information.")            
        return harmonization_data
    
    
    def _rebinTimeSeries(self, var, harmonization_data):
        """ 
        Helper function that handles the rebinning of timeseries data into larger timebins.
        Accounts for variable's aggregation rules. However, will not convert rates to quantities.
        Allows for fractional splitting of an old bin between two new bins, but this will break uncertainty calculations.
        
        Parameters
        ----------
        var: Variable
            Variable whose timeseries must be rebinned.
        harmonization_data: TimeHarmonizationData
            Populated `TimeHarmonizationData` object from which rebinning instructions are taken.
        
        Returns
        -------
        np.ndarray
            Array containing the new timeseries.
        
        Raises
        ------
        NotImplementedError
            If the timeseries must be rebinned to a single bin.
        """
        low_index, high_index       = harmonization_data.low_index, harmonization_data.high_index
        low_fraction, high_fraction = harmonization_data.low_fraction, harmonization_data.high_fraction
        factor                      = harmonization_data.upsample_factor
        
        new_values = np.zeros( int((high_index-low_index)/factor) )
        if len(new_values)==1:
            raise NotImplementedError("Aggregating time series data into a single bin currently not supported. - Timesum calling not implemented yet") ###!!!

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
        if var.aggregation_rule == "average" or var.is_rate: 
            new_values /= factor
        return new_values
    
    def _checkDependencyTimeHarmony(self, dependencies):
        """ 
        Helper function of `ensureDependencyTimeHarmony` that checks for a set of dependencies if they are time-harmonious by default,
        which allows to short-circuit the time harmony calculations.
        
        Parameters
        ----------
        dependencies: dict[str, Variable]
            Dictionary containing the dependencies to check time-harmony for.
        
        Returns
        tuple[bool, tuple[timestamp, timstap, int] or None]
            Returns a tuple containing:
            - A boolean indicating whether all dependencies are time-harmonious.
            - If dependencies are time-harmonious, returns their common start time, end time and timestep in a tuple.
              If dependencies are not time-harmonious, this field will be `None`.
        """
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
        """ 
        Helper function that ensures a set of harmonizations has a common start and end time.
        
        Detects common start and end time of all `TimeHarmonizationData` objects in the given dictionary, and prunes each timeseries
        such that they all span exactly this interval.
        
        Parameters
        ----------
        harmonized_dataset: dict[str, TimeHarmonizationData]
            Dictionary containing each dependency's `TimeHarmonizationData` object.
        new_timestep: int
            New timestep in seconds.
        
        Returns
        -------
        tuple[ dict[str, TimeHarmonizationData], tuple[timestamp, timestamp, int] ]
            Tuple containing updated `TimeHarmonizationData` objects, and a tuple containing the common start time, end time and timestep.
        """
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
        return harmonized_dataset, (common_start_time, common_last_time, new_timestep)
        



