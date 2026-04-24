from input_handler_modules import EquationTreeReader
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine
from datahandler import DataHandler

import numpy as np
from dataclasses import dataclass
from typing import Union, Optional


@dataclass
class RunResult:
    """ Stores the results, either calculation outputs or the error code, for a single run """
    identifier: Union(str, int)  #String which stores the identifier of the run
    succeeded: bool              #Stores whether the run succeeded or not
    data: dict                   #Stores results of the run
    fail_code: Optional(str)    #Stores which datacheck caused preprocessing to fail
    
    def __str__(self):
        temp = f"Run identifier: {self.identifier} \nRun succeeded: {self.succeeded}\n"
        if not self.succeeded:
            temp += f"Error code: {self.fail_code}\n"
        temp += f"Data fields: {list(self.data.keys())}"
        return temp
    
    
class Results:
    """ Stores all run results, as well as the staged data to be turned into the next run result """
    def __init__(self):
        self.groups = {
            "default": {
                "run_results"     : [],
                "run_identifiers" : [],
                "num_runs"        : 0}
            }
        
        #self.num_runs           = 0             #Stores the number of runs this result object contains
        #self.run_identifiers    = []            #For each run, can be used to an identifier (such as the date)
        #self.run_results        = []            #Stores RunResult object for each run
        self.unique_results     = {}            #Stores results that only need to be stored once
        #self.averages_effective_lengts = {}     #Stores effective lengths N for data that is stored as average: A -> (A*(N-1) + value)/N for the N'th result
        self.staged_data        = {}           #Staged data dictionary to be populated in the present run
       
    def _createGroup(self, name):
        ###!!!
        #Check if group name exists already
        if self.groups.get(name) is not None:
            raise ValueError("Error creating new results group: group name {name} already taken. Please check result-writing workflow.")
        self.groups[name] = {
            "run_results"     : [],
            "run_identifiers" : [],
            "num_runs"        : 0}
    
    def _getGroup(self, name):
        ###!!!
        group = self.groups.get(name)
        if group is None:
            raise ValueError("Results group '{name}' does not exist.")
        return group
    
    def add(self, key: str, value):
        """ Append a new result to the staged data for a given key, automatically creates a field for the key if it does not exist yet """
        if key not in self.staged_data:
            self.staged_data[key] = value
            return
        
        column = self.staged_data[key]
        if not isinstance(column, list):
            column = list(column)
        column.append(value)
                
    def createRunResult(self, succeeded=True, identifier=None, fail_code=None, group="default"):
        """ Compiles a RunResult object for the given run, resets the staged data dictionary """
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
        ###!!!
        group = self._getGroup(group)
        
        if not ( (identifier is None) ^ (index is None)):
            raise ValueError("Cannot get result: provide either an identifier or an index, not both.")

        #Get index of the identifier        
        if identifier is not None:
            index = np.where(np.array(group["run_identifiers"]) == identifier)[0][0]
        
        if index > len(group["run_results"]):
            raise ValueError(f"Cannot get result: index {index} > length results {len(group['run_results'])}")
        
        return group["run_results"][index]
            

    def getUniqueResult(self, name):
        """ Retrieve results from the unique result dictionary """
        return self.unique_results[name]
    
    def getResultList(self, name, give_identifiers=False, group="default"):
        """ Creates a list with all results corresponding to the given name from the list of runresult objects, also returns the identifiers. """
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
        else:
            return series
    
    def getResultArray(self, name, decimals=5, give_identifiers=False, group="default"):
        """ Same as the getResultSeries function, but returns the results as a numpy array """
        series, identifiers = self.getResultList(name, give_identifiers=True, group=group)
        series = np.asarray(series)
        if not series.dtype is np.object_:
            series = np.round(series, decimals=decimals)
        
        if give_identifiers:
            return series, identifiers
        else:
            return series
    
    def getAverageResult(self, name, decimals=5, give_identifiers=False, group="default"):
        """ Gets the average result over all runs """
        series, identifiers = self.getResultArray(name, decimals=decimals, give_identifiers=True, group=group)
        
        if np.ndim(series) > 1:
            series = np.average(series, axis=0)
        else:
            series = np.average(series)
            
        series = np.round(series, decimals=decimals)

        if give_identifiers:
            return series, identifiers
        else:
            return series
    
    def getFails(self, failcode=None, give_identifiers=True, group="default"):
        """ Gets a list of failcodes and identifiers of all failed runs """
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
        
        if give_identifiers:
            return fails, identifiers
        return fails
    
    def getSuccessBooleans(self, as_array=False, failcode=None, group="default"):
        """ Gets lists of all identifiers and an array of booleans on whether the run was a success, optionally filter for failcodes.
            Can be used to identify seasonal depenency of filtering hits. In case of large differences between march-october and october-march, 
            Check whether PVLIB handles daylight savings time in correspondence to how the dataset handles it. """
        group = self._getGroup(group)
        
        all_identifiers, success_bools = [], []
        
        for result in group["run_results"]:
            all_identifiers.append(result.identifier)
            if result.succeeded:
                success_bools.append(True)
                continue
            elif failcode is None:
                success_bools.append(False)
                continue
            elif result.failcode == failcode:
                success_bools.append(False)
            else:
                success_bools.append(False)
        
        if as_array:
            return np.array(success_bools), np.array(all_identifiers)
        return success_bools, all_identifiers
        
    
    def summariseFails(self, group="default"):
        """ Lists the number of times a failcode occurs, the number of times the job succeeded, and the total number of job calls """
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
    def __init__(self):
        self.data = DataHandler()
        
        self.equation_engine    = None
        self.calculation_engine = None
        self.uncertainty_engine = None
        self.time_engine        = None
        
        self.preprocessing = None
        self.main = None
        
        self.variables = None
        self.derived_variables_names = None
        self.var_backend_pointers = None
        
        self.blacklist = {}
    
        self.results = Results()
        
        ##computational control flow booleans
        self.initialized_eq_tree         = False
        self.initialized_engines         = False
        self.ready_for_execution         = False
        self.basic_variables_validated   = False
        self.backend_variables_populated = False
    
    def loadEquationTree(self, filepath):
        """ Load an equation tree from a text file into the job handler. 
        Compiles the variable registry, populates dependencies, builds executables and checks equation tree consistency """
        self.equation_tree_reader = EquationTreeReader()
        self.variables, self.var_backend_pointers = self.equation_tree_reader.parse(filepath)
        del self.equation_tree_reader
        
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
        
    def prepareEngines(self):
        """ Prepares the calculation, uncertainty, and time engine using the loaded variable registry """
        if not self.initialized_eq_tree:
            raise ValueError("Cannot initialize engines, since no equation tree appears to be loaded to the job handler. Please check your operations.")
        self.time_engine = TimeEngine()
        self.calculation_engine = CalculationEngine(variables          = self.variables, 
                                                    time_engine        = self.time_engine, 
                                                    equation_engine    = self.equation_engine)
        self.uncertainty_engine = UncertaintyEngine(variables          = self.variables,
                                                    equation_engine    = self.equation_engine, 
                                                    calculation_engine = self.calculation_engine, 
                                                    time_engine        = self.time_engine)
        self.initialized_engines = True
    
    
    def resetVariableRegistry(self):
        """ Resets all variables that are not explicitly hard-coded in the equation tree input """
        for var in self.variables.values():
            if not var.is_hardcoded:
                var.reset()
        
        self.basic_variables_validated = False
        if len(self.var_backend_pointers) != 0:
            self.backend_variables_populated = False
        
    def populateVariablesFromBackend(self, day=None, blacklist=[], reset_registry=True):
        """ For each variable in the backend_pointers dictionary this function will populate the variables with the requested window from their data backend """
        if reset_registry:
            self.resetVariableRegistry()
            
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
                                                                       as_array     = True, 
                                                                       blacklist    = blacklist)
            var = self.variables[var_name]
            var.values     = data
            var.setTimeData((start_time, end_time, timestep))
        
        self.backend_variables_populated = True

    def validateBasicVariables(self):
        """ Wrapper for calculation engine function of the same name, also updates the relevant flag """
        #if self.variables is None:
        #    raise ValueError("Validation of basic variables failed: no existing variable registry found.")
        #if not self.backend_variables_populated:
        #    print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
            
        self.calculation_engine.validateBasicVariables(equation_engine=self.equation_engine, variables=self.variables)
        self.basic_variables_validated = True
        return
    
    def prepareForExecution(self):
        ###!!!
        if not self.initialized_engines:
            self.prepareEngines()
        
        if self.preprocessing is not None:
            self.preprocessing(self)
        self.blacklist = self.data.compileBlacklist()
        
        if self.main is None:
            raise ValueError("No main jobscript is provided to the job handler. Please provide a main function under JobHandler.main")
        
        self.ready_for_execution = True
        
    def _blacklistResolver(self, day, blacklist_mode):
        ###!!!        
        #Handle the 'fail' blacklist mode
        if blacklist_mode == "fail":
            #Day is None
            if day is None and len(self.blacklist)>0:
                fail_code = "Any day blacklisted"
                return False, fail_code, []
            #Day is not None
            preprocessing_error = self.blacklist.get(day)
            if preprocessing_error is not None:
                fail_code = preprocessing_error[0]
                return False, fail_code, []
            #Else: passed blacklist checks
            return True, None, []
        
        #Handle the 'mask' blacklist mode
        elif blacklist_mode == "mask":
            return True, None, list(self.blacklist.keys())
        #Handle the 'ignore' blacklist mode
        elif blacklist_mode == "ignore":
            return True, None, []
        else:
            raise ValueError(f"Error: blacklist execution mode {blacklist_mode} not recognized.")
    
    
    def execute(self, day=None, identifier=None, blacklist_mode="fail", results_group="default"):
        """ Main job execution function, handles correct order of operations for preprocessing, execution and storage of results """
        #Ensure engines are staged for execution
        if not self.ready_for_execution:
            self.prepareForExecution()
        
        #Resolve blacklisting of days
        succeeded, fail_code, blacklist = self._blacklistResolver(day, blacklist_mode)
        if not succeeded:
            self.results.createRunResult(succeeded=False, identifier=identifier, fail_code=fail_code, group=results_group)
            return
        
        #Populate and validate equation tree
        self.populateVariablesFromBackend(day=day, blacklist=blacklist)
        self.validateBasicVariables()
        
        #Execute job
        succeeded, fail_code = self.main(self, identifier=identifier)

        #Create run result
        self.results.createRunResult(succeeded=succeeded, identifier=identifier, fail_code=fail_code, group=results_group)
        return
    
    def _resolve_arg(self, arg):
        """ Replaces function argument string referring to function attribute by the value of this attribute at time of calling """
        if isinstance(arg, str) and arg.startswith("var."):
            parts = arg.split(".")
            obj = self.variables[parts[1]]
            for attr in parts[2:]:
                obj = getattr(obj, attr)
            return obj
        return arg
        
    def resolve_args(self, args):
        """ Resolves all function arguments such that strings are replaced by attributes they refer to """
        if not isinstance(args, tuple):
            return self._resolve_arg(args)
        return tuple(self._resolve_arg(arg) for arg in args)
    
    def addDataFrame(self, df):
        """ Adds a pandas dataframe to the data backend """
        self.data.addDataFrame(df)
        
    #################################################################
    
    def store(self, name, arg):
        """ job task to store attribute 'arg' under name 'name' each job call """
        arg = self._resolve_arg(arg)
        self.results.add(name, arg)
    
    def storeUniqueResult(self, name, arg, replace=False):
        """ Store a single unique result instead of a full timeseries of the result """
        if name in self.results.unique_results.keys() and not replace:
            return
        arg = self._resolve_arg(arg)
        self.results.unique_results[name] = arg
    
    def evaluateVariable(self, var, *args, **kwargs):
        """ Wrapper for the calculation engine function of the same name """
        #if not self.basic_variables_validated:
        #    self.validateBasicVariables()
        
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate non-existing variable '{var}'.")
        self.calculation_engine.evaluateVariable(var, *args, **kwargs)
        
    def evaluateAllVariables(self):
        """ Wrapper for the calculation engine function of the same name """
        #if not self.basic_variables_validated:
        #    self.validateBasicVariables()
        for var_name in self.derived_variables_names:
            self.calculation_engine.evaluateVariable(self.variables[var_name])
            
    def prepareAllDirectUncertainties(self):
        """ Wrapper for the uncertainty engine function of the same name """
        self.uncertainty_engine.prepareAllDirectUncertainties()
    
    def prepareDownTreeDirectUncertainties(self, var):
        """ Wrapper for the uncertainty engine function of the same name """
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.prepareDownTreeDirectUncertainties(var)
    
    def calculateTotalUncertainty(self, var, *args, **kwargs):
        """ Wrapper for the uncertainty engine function of the same name """
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.calculateTotalUncertainty(var, *args, **kwargs)
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    