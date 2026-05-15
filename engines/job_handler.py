from engines.input_handler_modules import EquationTreeReader
from engines.equation_engine import EquationEngine
from engines.calculation_engine import CalculationEngine
from engines.uncertainty_engine import UncertaintyEngine
from engines.time_engine import TimeEngine
from engines.datahandler import DataHandler

import numpy as np
from dataclasses import dataclass
from typing import Union, Optional


@dataclass
class RunResult:
    """ 
    Dataclass used to store the results of a single job execution.
    
    Attributes
    ----------
    identifier: object
        Any desired (preferably unique) identifier of the run result.
    succeeded: bool
        Whether the result concerns a successful job execution.
    data: dict[str, object]
        Dictionary containing all data specified to be stored in the ``main`` script.
    fail_code: str, optional
        Optional string describing reason of execution failure.
    """
    identifier: object          #String which stores the identifier of the run
    succeeded: bool             #Stores whether the run succeeded or not
    data: dict                  #Stores results of the run
    fail_code: Optional(str)    #Stores which datacheck caused preprocessing to fail
    
    def __str__(self):
        temp = f"Run identifier: {self.identifier} \nRun succeeded: {self.succeeded}\n"
        if not self.succeeded:
            temp += f"Error code: {self.fail_code}\n"
        temp += f"Data fields: {list(self.data.keys())}"
        return temp
    
    
class Results:
    """ 
    Stores all run results, and the currently staged data to be flushed to next run result. 
    Results are grouped, to allow separation of conceptually distinct result types (e.g. daily results or yearly results).
    Class contains methods to retrieve arrays of desired results by group (e.g. all calculated daily performance ratios).
    
    Attributes
    ----------
    groups: dict[str: dict]
        Each group has the structure:
        {
            "run_results": list[RunResult],
            "run_identifiers": list[object],
            "num_runs": int
        }
    unique_results: dict[str, object]
        Dictionary containing stored unique results, that are only stored if their key does not yet exist in this dictionary.
    staged_data: dict[str, object]
        Dictionary containing all staged results. Populated during execution loop by results specified in ``JobHandler.main`` function.
        Staged data is flushed to a RunResult upon completion of a single job execution.
    """
    def __init__(self):
        self.groups = {
            "default": {
                "run_results"     : [],
                "run_identifiers" : [],
                "num_runs"        : 0}  
            }
        self.unique_results     = {}           #Stores results that only need to be stored once
        self.staged_data        = {}           #Staged data dictionary to be populated in the present run
       
    def _createGroup(self, name):
        """
        Initializes a new group in the ``self.groups`` dictionary, if the group name is not already taken.
        
        Parameters
        ----------
        name: str
            Name of the new group.
        
        Raises
        ------
        KeyError
            If the group name already exists as a key in ``self.groups``.
        """
        #Check if group name exists already
        if self.groups.get(name) is not None:
            raise KeyError("Error creating new results group: group name {name} already taken. Please check result-writing workflow.")
        self.groups[name] = {
            "run_results"     : [],
            "run_identifiers" : [],
            "num_runs"        : 0}
    
    def _getGroup(self, name):
        """
        Retrieves a group from results.

        Parameters
        ----------
        name : str
            Name of the group to be retrieved

        Raises
        ------
        KeyError
            If the requested group does not exist.

        Returns
        -------
        dict[str: dict]
            Empty group dictionary of structure:
            {
                "run_results": [],
                "run_identifiers": [],
                "num_runs": 0
            }
        """
        group = self.groups.get(name)
        if group is None:
            raise KeyError("Results group '{name}' does not exist.")
        return group
    
    def add(self, key, value):
        """ 
        Append a new result to the staged data under a given key, automatically creates a field for the key if it does not exist yet.
        
        Parameters
        ----------
        key: str
            Key under which the data must be stored.
        value: object
            Data to be stored under the desired key.
        """
        if key not in self.staged_data:
            self.staged_data[key] = value
            return
        
        column = self.staged_data[key]
        if not isinstance(column, list):
            column = list(column)
        column.append(value)
                
    def flushRunResult(self, succeeded=True, identifier=None, fail_code=None, group="default"):
        """ 
        Compiles and stores a RunResult from the staged data under desired group, upon which the staged data is reset.
        
        Parameters
        ----------
        succeeded: bool, default=True
            Whether the job execution has succeeded.
        identifier: object, default=None
            Identifier of the run. If ``None``, an identifier will be assigned.
        fail_code: str or None
            Optional explanation for why the job execution failed.
        group: str, default='default'
            Group name under which the results must be stored.
        """
        #Retrieve results group        
        if self.groups.get(group) is None:
            self._createGroup(group)
        group = self.groups.get(group)
        
        #Define identifier is none is given
        if identifier is None:
            identifier = group["num_runs"]
        
        #Compile result from staged data
        result = RunResult(identifier   = identifier,
                           succeeded    = succeeded,
                           data         = self.staged_data,
                           fail_code    = fail_code)
        #Clear staged data
        self.staged_data = {}
        
        #Append result to group
        group["run_results"].append(result)
        group["run_identifiers"].append(identifier)
        group["num_runs"] += 1


    def getResult(self, identifier=None, index=None, group="default"):
        """
        Retrieve result from a group based on identifier, or integer index in the group.

        Parameters
        ----------
        identifier: str or None
            Identifier of the result to be retrieved. Mutually exclusive with using ``index``.
        index: int or None
            Integer index in group of the result to be retrieved. Mutually exclusive with using ``identifier``.
        group: str, default='default'
            Group name in which to look for the result.

        Raises
        ------
        ValueError
            If both an identifier and index are given in the function input.
        KeyError
            If the requested run identifier is not present in the results group.
        IndexError
            If the requested index exceeds the number of results present in the group.

        Returns
        -------
        RunResult
            Requested run result.
        """
        group = self._getGroup(group)
        
        if not ( (identifier is None) ^ (index is None)):
            raise ValueError("Cannot get result: provide either an identifier or an index, not both.")

        #Get index of the identifier        
        if identifier is not None:
            index = group["run_identifiers"].index(identifier)
        
        if index > len(group["run_results"]):
            raise IndexError(f"Cannot get result: index {index} > length results {len(group['run_results'])}")
        
        return group["run_results"][index]
            

    def getUniqueResult(self, name):
        """ 
        Retrieve result from the unique result dictionary.
        
        Parameters
        ----------
        name: str
            Name of the unique result to be returned.
        
        Raises
        ------
        KeyError
            If the requested key does not exist.
        
        Returns
        -------
        object
            The data stored under the given name.
        """
        return self.unique_results[name]
    
    def getResultList(self, name, give_identifiers=False, group="default"):
        """ 
        Returns a list of all results of a given key in a group.
        Optionally also returns their corresponding identifiers, in which case the lists are matched by index.
        
        Parameters
        ----------
        name: str
            Key under which the requested data has been stored.
        give_identifiers: bool, default=False
            Whether a list of corresponding identifiers should be returned as well.
        group: str, default='default'
            Group from which results should be retrieved.
        
        Returns
        -------
        list 
            If ``give_identifiers=False``, returns a list of results.
        tuple[list, list]
            If ``give_identifiers=True``, returns lists of results and corresponding identifiers, matched by index.
        """
        group = self._getGroup(group)
        
        series, identifiers = [], []
        for result in group["run_results"]:
            datapoint = result.data.get(name)
            if datapoint is None:
                continue
            series.append(datapoint)
            identifiers.append(result.identifier)
        
        if len(series) == 0:
            print(f"WARNING: No results with name {name} were found")
        
        if give_identifiers:
            return series, identifiers
        return series
    
    def getResultArray(self, name, decimals=None, give_identifiers=False, group="default"):
        """ 
        Same as the ``results.getResultSeries`` function, but returns the results as a numpy array.
        Has the added functionality that numeric results can be immediately rounded to desired number of decimals.
        
        Parameters
        ----------
        name: str
            Key under which the requested data has been stored.
        decimals: int or None, default=None
            Number of decimals of numeric results that should be returned. Passing None will skip rounding step.
        give_identifiers: bool, default=False
            Whether a list of corresponding identifiers should be returned as well.
        group: str, default='default'
            Group from which results should be retrieved.
        
        Returns
        -------
        list 
            If ``give_identifiers=False``, returns an array of results.
        tuple[list, list]
            If ``give_identifiers=True``, returns array of results and list of corresponding identifiers, matched by index.
        """
        series, identifiers = self.getResultList(name, give_identifiers=True, group=group)
        series = np.asarray(series)
        if not series.dtype is np.object_ and decimals is not None:
            series = np.round(series, decimals=decimals)
        
        if give_identifiers:
            return series, identifiers
        return series
    
    def getAverageResult(self, name, decimals=None, give_identifiers=False, group="default"):
        """ 
        Returns the average of all results of a given key in a group.
        Optionally also returns the identifiers of results included in the average.
        
        Parameters
        ----------
        name: str
            Key under which the requested data has been stored.
        decimals: int or None, default=None
            Number of decimals of numeric results that should be returned. Passing None will skip rounding step.
        give_identifiers: bool, default=False
            Whether a list of corresponding identifiers should be returned as well.
        group: str, default='default'
            Group from which results should be retrieved.
            
        Returns
        -------
        float 
            If ``give_identifiers=False``, returns an average result.
        tuple[float, list]
            If ``give_identifiers=True``, returns average result and list of identifiers included in the average.
        """
        series, identifiers = self.getResultArray(name, decimals=decimals, give_identifiers=True, group=group)
        
        if np.ndim(series) > 1:
            series = np.average(series, axis=0)
        else:
            series = np.average(series)
        
        if decimals is not None:
            series = np.round(series, decimals=decimals)

        if give_identifiers:
            return series, identifiers
        return series
    
    def getFails(self, failcode=None, group="default"):
        """ 
        Gets lists of failcodes and identifiers of all failed runs in a group.
        If a failcode is specified, returned results only include this failcode.
        
        Parameters
        ----------
        failcode: str or None, default=None
            Failcode to filter for. If None, all failcodes are included in result.
        group: str, default=None
            Group from which results should be retrieved.
        """
        group = self._getGroup(group)
        
        fails, identifiers = [], []
        for result in group["run_results"]:
            if result.succeeded:
                continue
            if failcode is None:
                fails.append(result.fail_code)
                identifiers.append(result.identifier)
                continue
            elif result.failcode == failcode:
                fails.append(result.fail_code)
                identifiers.append(result.identifier)
        return fails, identifiers
    
    def getSuccessBooleans(self, as_array=False, failcode=None, group="default"):
        """ 
        Gets lists of all identifiers and an array of booleans on whether the run was a success or not. 
        By passing a failcode, all fails without that fail code will also be considered a succes.
        
        Can be used to identify seasonal depenency of filtering hits. In case of large differences between march-october and october-march, 
        check whether PVLIB handles daylight savings time in correspondence to how the dataset handles it. 
        
        Parameters
        ----------
        as_array: bool, default=False
            Whether the lists of success booleans and identifiers should be returned as arrays.
        failcode: str or None, default=None
            Can be used to only consider a specific failcode a fail.
            If ``None``, all failcodes will be considered a fail.
        group: str, default='default'
            Group from which results should be retrieved.
        
        Returns
        -------
        tuple[list, list]
            If ``as_array=False``. First list contains success booleans, second list contains corresponding identifiers.
        tuple[np.ndarray, np.ndarray]
            If ``as_array=True``. First array contains success booleans, second array contains corresponding identifiers.
        """
        group = self._getGroup(group)
        
        all_identifiers, success_bools = [], []
        
        for result in group["run_results"]:
            all_identifiers.append(result.identifier)
            if result.succeeded:
                success_bools.append(True)
            elif failcode is None:
                success_bools.append(False)
            elif result.failcode == failcode:
                success_bools.append(False)
            else:
                success_bools.append(True)
        
        if as_array:
            return np.array(success_bools), np.array(all_identifiers)
        return success_bools, all_identifiers
        
    
    def summariseFails(self, group="default"):
        """ 
        Lists the number of times a failcode occurs, the number of times the job succeeded, and the total number of job calls.
        
        Parameters
        ----------
        group: str, default='default'
            Group from which the fails should be summarized.
        """
        fails, identifiers = self.getFails(group=group)
        num_runs = self.groups[group]["num_runs"]
        successful_executions = num_runs - len(fails)
        
        unique_fails = list(set(fails))
        print("Summarising calculation failures:")
        print(f"Successful executions: {successful_executions} out of {num_runs} total executions")
        for fail in unique_fails:
            print(f"Error code {fail} occurred {fails.count(fail)} times")
        print()
            
            

class JobHandler:
    """ 
    Main execution orchestrator and central node between jobscript, datastorage and equation tree.
    The JobHandler owns the equation tree, various execution engines, data backend and results storage.
    The JobHandler also ensures engines and equation tree are properly initialized before execution.
    The class furthermore contains wrapper functions for calculation and storage functionality commonly used in execution scripts.
    
    Attributes
    ----------
    variables: dict[str, Variable]
        Central variable registry. Each variable is linked to its dependencies, creating the equation tree.
    main: function, default = None
        Function that is executed upon calling job.execute(). To be defined by user.
    preprocessing: function, default = None
        Function that is executed as a preprocessing step. Used to clean data and check consistency, and potentially build a blacklist of days. To be defined by user. Optional.
    data: DataHandler
        Central data backend containing pandas dataframes, data retrieval routines, cleaning functions and a registry of blacklisted days.
    var_backend_pointers: dict[str, str]
        Dictionary connecting the names in the variable registry to the name of their corresponding data column in the backend, if they have one.
    results: Results
        Central results storage, where the results of each execution are gathered and stored.
    equation_engine: EquationEngine
        Builds equation tree executables. Can also be used for directly accessing complete equation engine functionality.
    calculation_engine: CalculationEngine
        Evaluates variables in the equation tree. Can also be used for directly accessing complete calculation engine functionality.
    uncertainty_engine: UncertaintyEngine
        Evaluates uncertainty in the equation tree. Can also be used for directly accessing complete uncertainty engine functionality.
    time_engine: TimeEngine
        Handles time alignment of data. Can also be used for directly accessing complete time engine functionality.
    equation_tree_reader: EquationTreeReader
        Text parser that can convert an equation tree text file to a dictionary of partially-populated variables. An equation engine is required for proper initialization of equation tree.
    """
    def __init__(self):
        self.data = DataHandler()
        
        self.equation_engine      = None
        self.calculation_engine   = None
        self.uncertainty_engine   = None
        self.time_engine          = None
        
        self.equation_tree_reader = None
        
        self.preprocessing = None
        self.main = None
        
        self.variables = None
        self.derived_variables_names = None
        self.var_backend_pointers = None
        
        self.results = Results()
        
        ##computational control flow booleans
        self.initialized_eq_tree         = False
        self.initialized_engines         = False
        self.ready_for_execution         = False
        self.basic_variables_validated   = False
        self.backend_variables_populated = False
    
    def loadEquationTree(self, filepath):
        """ 
        Load an equation tree from a text file into the JobHandler.
        Compiles the variable registry, populates dependencies, builds executables and checks equation tree consistency.
        In essence, this function executes equation engine functionality to properly prepare the equation tree.
        
        Parameters
        ----------
        filepath: str
            Path to the text file containing the equation tree.
        """
        if self.equation_tree_reader is None:
            self.equation_tree_reader = EquationTreeReader()
        self.variables, self.var_backend_pointers = self.equation_tree_reader.parse(filepath)
        
        #If no variables need to be populated from csv's: set csv_variables_populated flag to True
        if len(self.var_backend_pointers) == 0:
            self.backend_variables_populated = True
        
        self.equation_engine = EquationEngine(self.variables)
        self.derived_variables_names = self.equation_engine.derived_variables
        
        self.equation_engine.checkEquationTreeConsistency(silent=True)
        self.equation_engine.populateEquationTreeDependencies()
        self.equation_engine.populateEquationTreeTimeSumSettings()
        self.equation_engine.buildEquationTreeExecutables()
        
        self.initialized_eq_tree = True
        
    def _prepareEngines(self):
        """ 
        Prepares the calculation, uncertainty and time engines using the loaded variable registry.
        Calculation and uncertainty engines rely on an initialized equation engine with an internal variable registry,
        and are therefore only initialized after the JobHandler has built a variable registry and initialized an equation engine.
        Furthermore, due to some elements and attributes being lazily calculated or created during job runtime, engines sometimes require each others' functionality.
        This function ensures engines are properly initialized.
        
        Raises
        ------
        RuntimeError
            If no equation tree has been initialized, since an equation engine is given to the CalculationEngine and UncertaintyEngine.
        """
        if not self.initialized_eq_tree:
            raise RuntimeError("Cannot initialize engines, since no equation tree appears to be loaded to the job handler. Please check your operations.")
        self.time_engine = TimeEngine()
        self.calculation_engine = CalculationEngine(time_engine        = self.time_engine, 
                                                    equation_engine    = self.equation_engine)
        self.uncertainty_engine = UncertaintyEngine(equation_engine    = self.equation_engine, 
                                                    calculation_engine = self.calculation_engine)
        self.initialized_engines = True
    
    
    def _resetVariableRegistry(self):
        """ 
        Resets the data of all variables in the equation tree, whose values are not explicitly hard-coded in the equation tree input. 
        """
        for var in self.variables.values():
            if not var.is_hardcoded:
                var.reset()
        
        self.basic_variables_validated = False
        if len(self.var_backend_pointers) != 0:
            self.backend_variables_populated = False
        
    def _populateVariablesFromBackend(self, day=None, blacklist=[]):
        """ 
        This function populates those variables in the equation tree with the requested window from their respective columns in the backend.
        Upon calling this function, the state of variables that are not explicitly hardcoded in the equation tree input will be reset.
        Variables are populated with data of a specific day, or with all data from their backend. 
        Days given as blacklisted, will be masked: their data will be included as zeroes in the equation tree.
        
        Parameters
        ----------
        day: datetime.date or None
            Day that is to be loaded into equation tree. If ``None``, all data for this variable will be loaded. 
        blacklist: list[datetime.date]
            Data of blacklisted days loaded into equation tree will be included exclusively as zeroes.
        """
        #Clean the variable registry
        self._resetVariableRegistry()
            
        for var_name, column_name in self.var_backend_pointers.items():
            #column_name is either CSV.column_name or CSV.coupled_name.column_name
            #column_name is of the form "CSV.___", we strip the first 4 characters
            column_name = column_name[4:]
            #Split coupled and column names
            parts = column_name.split(".")
            if len(parts)>1:
                column_name  = parts[1]
                coupled_name = parts[0]
            else:
                column_name  = parts[0]
                coupled_name = None
            
            data, start_time, end_time, timestep = self.data.getColumn(name         = column_name, 
                                                                       coupled_name = coupled_name, 
                                                                       day          = day, 
                                                                       blacklist    = blacklist)
            var = self.variables[var_name]
            var.values     = data
            var.setTimeData((start_time, end_time, timestep))
        
        self.backend_variables_populated = True

    def _validateBasicVariables(self):
        """ 
        Wrapper for CalculationEngine function of the same name.
        Validation of basic variables ensures all basic variables have populated values, or their values can be inferred from populated variables uptree and an equation.
        """
        #if self.variables is None:
        #    raise ValueError("Validation of basic variables failed: no existing variable registry found.")
        #if not self.backend_variables_populated:
        #    print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
            
        self.calculation_engine.validateBasicVariables(variables=self.variables, equation_engine=self.equation_engine)
        self.basic_variables_validated = True
        return
    
    def _prepareForExecution(self):
        """ 
        Function ensuring the state of the JobHandler is ready for execution and any preprocessing function given to it is executed.
        """
        if not self.initialized_engines:
            self._prepareEngines()
        
        if self.preprocessing is not None:
            self.preprocessing(self)
        
        if self.main is None:
            raise ValueError("No main jobscript is provided to the job handler. Please provide a main function under JobHandler.main")
        
        self.ready_for_execution = True
        
    def _blacklistResolver(self, day, blacklist_mode):
        """ 
        Helper function that resolves execution blacklist setting. 
        Tells the execution whether it can continue, and if so, with which blacklist of dates.
        
        Parameters
        ----------
        day: datetime.date or None
            Day that is executed, or None if the entire dataset is executed at once.
        blacklist_mode: str
            Either 'fail', 'mask' or 'ignore'. 
            -Fail: the execution fails if a blacklisted day is among the days executed.
            -Mask: execution will proceed, but all variable values during this day will be set to 0.
            -Ignore: blacklisting is ignored entirely.
        
        Raises
        ------
        ValueError
            If a blacklist_mode other than 'fail', 'mask' or 'ignore' is passed.
        
        Returns
        -------
        tuple[Bool, str or None, list[datetime.date]]
            Boolean expressing whether execution is allowed to proceed.
            String with explanation if execution is not allowed to proceed. ``None`` otherwise.
            List containing the blacklisted days to be be masked during execution, if ``mode=='mask'``, otherwise empty list.
        """
        #Handle the 'fail' blacklist mode
        if blacklist_mode == "fail":
            #Day is None
            if day is None and len(self.data.blacklist)>0:
                fail_code = "Any day blacklisted"
                return False, fail_code, []
            #Day is not None
            preprocessing_error = self.data.blacklist.get(day)
            if preprocessing_error is not None:
                fail_code = preprocessing_error[0]
                return False, fail_code, []
            #Else: passed blacklist checks
            return True, None, []
        
        #Handle the 'mask' blacklist mode
        elif blacklist_mode == "mask":
            return True, None, list(self.data.blacklist.keys())
        #Handle the 'ignore' blacklist mode
        elif blacklist_mode == "ignore":
            return True, None, []
        else:
            raise ValueError(f"Error: blacklist execution mode {blacklist_mode} not recognized.")
    
    
    def execute(self, day=None, identifier=None, blacklist_mode="fail", results_group="default"):
        """ 
        Main job execution function.
        Calls internal preparation functions, optionally handles preprocessing, resolves blacklisting handling,
        prepares equation tree state, executes user defined main function script.
        Results are stored in the user-specified results group.
        The user-defined ``main`` functions must return ``(bool, str or None)``.
        
        Parameters
        ----------
        day: datetime.date or None, defualt = None
            Day that is to be executed, or None if the entire dataset is executed at once.
        identifier: object, default = None
            Any deisred identifier of the run, used to access the result of the execution from the results.
            If no identifier is given, the results handler will assign an integer identifier corresponding to the number of runs in the result group.
        blacklist_mode: str, default = 'fail'
            String defining the blacklist handling setting during execution.
            -Fail: the execution fails if a blacklisted day is among the days executed.
            -Mask: execution will proceed, but all variable values during this day will be set to 0.
            -Ignore: blacklisting is ignored entirely.
        results_group: str, default = 'default'
            Group to which the results of the run are flushed. Allows for easy separation of run results when compiling result arrays.
            
        Returns
        -------
        tuple [bool, str or None]
            Boolean indicates whether the execution was successful or aborted,
            optional string can be used to describe the reason for abortion.
        """
        #Ensure engines are staged for execution
        if not self.ready_for_execution:
            self._prepareForExecution()
        
        #Resolve blacklisting of days
        succeeded, fail_code, blacklist = self._blacklistResolver(day, blacklist_mode)
        if not succeeded:
            self.results.flushRunResult(succeeded=False, identifier=identifier, fail_code=fail_code, group=results_group)
            return
        
        #Populate and validate equation tree
        self._populateVariablesFromBackend(day=day, blacklist=blacklist)
        self._validateBasicVariables()
        
        #Execute job
        succeeded, fail_code = self.main(self, identifier=identifier)

        #Create run result
        self.results.flushRunResult(succeeded=succeeded, identifier=identifier, fail_code=fail_code, group=results_group)
        return
    
    def _resolve_arg(self, arg):
        """
        Resolve string-based variable references to their runtime values.

        If "arg" is a string starting with ``'var.'``, it is interpreted as a
        reference to a variable and its attributes. The reference is resolved
        by looking up the variable in ``self.variables`` and traversing the
        specified attributes.
    
        Parameters
        ----------
        arg : object
            Argument to resolve. If a string of the form ``'var.<name>.<attr>...'``,
            it is interpreted as a variable reference.
    
        Returns
        -------
        object
            Resolved value if ``arg`` is a variable reference, otherwise ``arg`` unchanged.
        """
        if isinstance(arg, str) and arg.startswith("var."):
            parts = arg.split(".")
            obj = self.variables[parts[1]]
            for attr in parts[2:]:
                obj = getattr(obj, attr)
            return obj
        return arg
        
    def resolve_args(self, args):
        """
        Resolve argument(s) by replacing variable references with their values.
        Applies _resolve_arg to each passed argument.
    
        Parameters
        ----------
        args : object or tuple of object
            Argument or tuple of arguments. Strings of the form
            ``'var.<name>.<attr>...'`` are interpreted as variable references.
    
        Returns
        -------
        object or tuple of object
            Resolved argument(s). If ``args`` is a tuple, a tuple of resolved
            values is returned; otherwise a single resolved value.
        """
        if not isinstance(args, tuple):
            return self._resolve_arg(args)
        return tuple(self._resolve_arg(arg) for arg in args)
    
    def addDataFrame(self, df):
        """ 
        Adds a pandas dataframe to the data backend.
        
        Parameters
        ----------
        df: pandas.DataFrame
            Dataframe to add to backend.
        """
        self.data.addDataFrame(df)
        
    #################################################################
    
    def store(self, name, arg):
        """ 
        Method to add any desired object to the results of a job execution run, under given name.
        To be used when writing results of the ``main`` script.
        Applies ``_resolve_arg`` to the passed argument.
        
        Parameters
        ----------
        name: str
            Dictionary key under which the desired data will be kept in the run result.
        arg: object
            Object to store. If a string of the form ``'var.<name>.<attr>...'``,
            it is interpreted as a variable reference, to be evaluated at runtime.
        """
        arg = self._resolve_arg(arg)
        self.results.add(name, arg)
    
    def storeUniqueResult(self, name, arg, replace=False):
        """ 
        Store a result globally under the given name, instead of under the result of a single run.
        Result is stored under the given name. If the name is already present in the global result, ``arg`` will not be stored unless replacement is specified.
        Applies ``_resolve_arg`` to the passed argument.
        
        Parameters
        ----------
        name: str
            Dictionary key under which the desired data will be stored.
        arg: object
            Object to store. If a string of the form ``'var.<name>.<attr>...'``,
            it is interpreted as a variable reference, to be evaluated at runtime.
        replace: bool, default = False
            Whether the unique result should be overwritten if it already exists.
        """
        if name in self.results.unique_results.keys() and not replace:
            return
        arg = self._resolve_arg(arg)
        self.results.unique_results[name] = arg
    
    def evaluateVariable(self, var, *args, **kwargs):
        """ 
        Wrapper for the ``CalculationEngine`` function of the same name.
        If ``var`` is a string, it is attempted to be retrieved from the internal variable registry.
        
        Parameters
        ----------
        var : str or Variable
            Variable name or Variable object to evaluate.
        *args
            Positional arguments passed to the calculation engine.
        **kwargs
            Keyword arguments passed to the calculation engine.
           
       Raises
       ------
       ValueError
           If ``var`` is a string and does not exist in the variable registry.
        """        
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate non-existing variable '{var}'.")
        self.calculation_engine.evaluateVariable(var, *args, **kwargs)
        
    def evaluateAllVariables(self):
        """ 
        Evaluates all variables in internal variable registry.
        """
        for var_name in self.derived_variables_names:
            self.calculation_engine.evaluateVariable(self.variables[var_name])
            
    def prepareAllDirectUncertainties(self):
        """ 
        Wrapper for the UncertaintyEngine function of the same name.
        """
        self.uncertainty_engine.prepareAllDirectUncertainties()
    
    def prepareDownTreeDirectUncertainties(self, var):
        """ 
        Wrapper for the UncertaintyEngine function of the same name.
        If ``var`` is a string, it is attempted to be retrieved from the internal variable registry.
        
        Parameters
        ----------
        var : str or Variable
            Variable name or Variable object to evaluate.
            
        Raises
        ------
        ValueError
            If ``var`` is a string and does not exist in the variable registry.
        """
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.prepareDownTreeDirectUncertainties(var)
    
    def calculateTotalUncertainty(self, var, *args, **kwargs):
        """ 
        Wrapper for the UncertaintyEngine function of the same name.
        If ``var`` is a string, it is attempted to be retrieved from the internal variable registry.
        
        Parameters
        ----------
        var : str or Variable
            Variable name or Variable object to evaluate.
        *args
            Positional arguments passed to the UncertaintyEngine.
        **kwargs
            Keyword arguments passed to the UncertaintyEngine.
            
        Raises
        ------
        ValueError
            If ``var`` is a string and does not exist in the variable registry.
        """
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.calculateTotalUncertainty(var, *args, **kwargs)
    
