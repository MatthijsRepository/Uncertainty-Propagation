from input_handler_modules import EquationTreeReader, CSVHandler, PandasCSVHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine
from datahandler import DataHandler
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Union, Optional


@dataclass
class RunResult:
    """ Stores the results, either calculation outputs or the error code, for a single run """
    identifier: Union(str, int)  #String which stores the identifier of the run
    succeeded: bool              #Stores whether the run succeeded or not
    data: dict                   #Stores results of the run
    fail_code: Optional(str)    #Stores which datacheck caused preprocessing to fail
    
    
class Results:
    """ Stores all run results, as well as the staged data to be turned into the next run result """
    def __init__(self):
        self.num_runs           = 0             #Stores the number of runs this result object contains
        self.run_identifiers    = []            #For each run, can be used to an identifier (such as the date)
        self.run_results        = []            #Stores RunResult object for each run
        self.unique_results     = {}            #Stores results that only need to be stored once
        #self.averages_effective_lengts = {}     #Stores effective lengths N for data that is stored as average: A -> (A*(N-1) + value)/N for the N'th result
        self.staged_data        = {}           #Staged data dictionary to be populated in the present run
        
    def add(self, key: str, value):
        """ Append a new result to the staged data for a given key, automatically creates a field for the key if it does not exist yet """
        if key not in self.staged_data:
            self.staged_data[key] = value
            return
        
        column = self.staged_data[key]
        if not isinstance(column, list):
            column = list(column)
        column.append(value)
                
    def createRunResult(self, succeeded=True, identifier=None, fail_code=None):
        """ Compiles a RunResult object for the given run, resets the staged data dictionary """
        if identifier is None:
            identifier = self.num_runs
        
        result = RunResult(identifier   = identifier,
                           succeeded    = succeeded,
                           data         = self.staged_data,
                           fail_code    = fail_code)
        self.run_results.append(result)
        
        self.run_identifiers.append(identifier)
        self.staged_data = {}
        self.num_runs    += 1
    
    def getUniqueResult(self, name):
        """ Retrieve results from the unique result dictionary """
        return self.unique_results[name]
    
    def getResultSeries(self, name, give_identifiers=False):
        """ Creates a list with all results corresponding to the given name from the list of runresult objects, also returns the identifiers. """
        series, identifiers = [], []
        for result in self.run_results:
            datapoint = result.data.get(name)
            if datapoint is None:
                continue
            series.append(datapoint)
            identifiers.append(result.identifier)
        
        if len(series) == 0:
            print(f"WARNING: No result with name {name} were found")
        
        if give_identifiers:
            return series, identifiers
        return series
    
    def getResultArray(self, name, decimals=5, give_identifiers=False):
        """ Same as the getResultSeries function, but returns the results as a numpy array """
        series, identifiers = self.getResultSeries(name, give_identifiers=True)
        series = np.asarray(series)
        if not series.dtype is np.object_:
            series = np.round(series, decimals=decimals)
        
        if give_identifiers:
            return series, identifiers
        return series
    
    def getAverageResult(self, name, decimals=5, give_identifiers=False):
        """ Gets the average result over all runs """
        series, identifiers = self.getResultArray(name, decimals=decimals, give_identifiers=True)
        
        if np.ndim(series) > 1:
            series = np.average(series, axis=0)
        else:
            series = np.average(series)
            
        series = np.round(series, decimals=decimals)

        if give_identifiers:
            return series, identifiers
        return series
    
    def getFails(self, failcode=None, give_identifiers=True):
        """ Gets a list of failcodes and identifiers of all failed runs """
        fails, identifiers = [], []
        for result in self.run_results:
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
    
    def getSuccessBooleans(self, as_array=False, failcode=None):
        """ Gets lists of all identifiers and an array of booleans on whether the run was a success, optionally filter for failcodes.
            Can be used to identify seasonal depenency of filtering hits. In case of large differences between march-october and october-march, 
            Check whether PVLIB handles daylight savings time in correspondence to how the dataset handles it. """
        all_identifiers, success_bools = [], []
        
        for result in self.run_results:
            all_identifiers.append(result.identifier)
            if result.succeeded:
                success_bools.append(True)
                continue
            elif failcode is None:
                success_bools.append(False)
                continue
            elif result.failcode == failcode:
                success_bools.append(False)
        
        if as_array:
            return np.array(success_bools), np.array(all_identifiers)
        return success_bools, all_identifiers
        
    
    def summariseFails(self):
        """ Lists the number of times a failcode occurs, the number of times the job succeeded, and the total number of job calls """
        fails, identifiers = self.getFails()
        successful_executions = self.num_runs - len(fails)
        
        unique_fails = list(set(fails))
        print("Summarising calculation failures:")
        print(f"Successful executions: {successful_executions} out of {self.num_runs} total executions")
        for fail in unique_fails:
            print(f"Error code {fail} occurred {fails.count(fail)} times")
        print()
            
            

import time

class JobHandler:
    def __init__(self):
        #self.CSV_handler = CSVHandler() ###!!! DEPRECATED
        
        self.t_prepro  = 0
        self.t_treepop = 0
        self.t_main    = 0
        
        self.data = DataHandler()
        
        self.equation_engine = None
        self.calculation_engine = None
        self.uncertainty_engine = None
        self.time_engine = None
        
        self.preprocessing = None
        self.main = None
        
        self.variables = None
        self.derived_variables_names = None
        self.var_csv_pointers = None
        self.csv_data = []
    
        self.results = Results()
        
        ##computational control flow booleans
        self.initialized_eq_tree       = False
        self.initialized_engines       = False
        self.has_csv_data              = False
        self.basic_variables_validated = False
        self.csv_variables_populated   = False
        
        ##
        return
    
    def loadEquationTree(self, filepath):
        """ Load an equation tree from a text file into the job handler. 
        Compiles the variable registry, populates dependencies, builds executables and checks equation tree consistency """
        self.equation_tree_reader = EquationTreeReader()
        self.variables, self.var_csv_pointers = self.equation_tree_reader.parse(filepath)
        del self.equation_tree_reader
        
        #If no variables need to be populated from csv's: set csv_variables_populated flag to True
        if len(self.var_csv_pointers) == 0:
            self.csv_variables_populated = True
        
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
        if len(self.var_csv_pointers) != 0:
            self.csv_variables_populated   = False
    
    
    
    def populateVariablesFromCSV(self, day=None, reset_registry=True):
        """ For each variable in the csv_pointers dictionary this function will populate the variables with the requested window from their data backend """
        if reset_registry:
            self.resetVariableRegistry()
            
        for var_name, column_name in self.var_csv_pointers.items():
            #Column name is of the form "CSV.name", we strip the first 4 characters
            column_name = column_name[4:]
            self.variables[var_name].values = self.data.getColumn(column_name, day=day, as_array=True)
            self.variables[var_name].addTimeStep(self.data.getTimeRange(name=column_name, day=day))
        self.csv_variables_populated = True
            
    
    def validateBasicVariables(self):
        """ Wrapper for calculation engine function of the same name, also updates the relevant flag """
        #if self.variables is None:
        #    raise ValueError("Validation of basic variables failed: no existing variable registry found.")
        #if not self.csv_variables_populated:
        #    print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
            
        self.calculation_engine.validateBasicVariables(equation_engine=self.equation_engine, variables=self.variables)
        self.basic_variables_validated = True
        return
    
    
    
    def runDay(self, day):
        self.populateVariablesFromCSV(day=day)
        
    
    
    def execute(self, day, identifier=None):
        """ Main job execution function, handles correct order of operations for preprocessing, storage of results and reinitializing the variable registry after completion """
        if self.main is None:
            raise ValueError("No main jobscript is provided to the job handler. Please provide a main function under JobHandler.main")
        if not self.initialized_engines:
            self.prepareEngines()
        
        #if not self.has_csv_data:
        #    print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
        
        t0 = time.time()
        #Perform the data preprocessing
        if self.preprocessing is not None:
            success, fail_code = self.preprocessing(self, identifier=identifier)
            #If preprocessing failed, we log this in the results, we also clear the loaded csv data
            if not success:
                self.results.createRunResult(succeeded=False, identifier=identifier, fail_code=fail_code)
                self.csv_data = []
                self.has_csv_data = False
                return
        self.t_prepro += time.time()-t0
        
        t0 = time.time()
        
        #Populate variables using loaded CSV data, and subsequently dump the csv data
        self.populateVariablesFromCSV(day=day)
        #self.csv_data = []
        #self.has_csv_data = False
        
        #Validate basic variables
        self.validateBasicVariables()
        
        self.t_treepop += time.time() - t0
        
        t0 = time.time()
        #Execute job
        success, fail_code = self.main(self, identifier=identifier)
        if not success:
            self.results.createRunResult(succeeded=False, identifier=identifier, fail_code=fail_code)
            self.csv_data = []
            self.has_csv_data = False
            return
        self.t_main += time.time()-t0
        
        
        #Create run result
        self.results.createRunResult(succeeded=True, identifier=identifier)
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
    
    #def _resolve_args(self, args):
    #    if not isinstance(args, tuple):
    #        return self._resolve_arg(args)
    #    return tuple(self._resolve_arg(arg) for arg in args)
    
    def resolve_args(self, args):
        """ Resolves all function arguments such that strings are replaced by attributes they refer to """
        if not isinstance(args, tuple):
            return self._resolve_arg(args)
        return tuple(self._resolve_arg(arg) for arg in args)
    
    def addDataFrame(self, df):
        """ Adds a pandas dataframe to the data backend """
        self.data.addDataFrame(df)
        
    
    def addCSVData(self, csv):
        """ Adds CSVData object or list of CSVData objects to internal registry """
        self.csv_data.append(csv)
        self.has_csv_data = True
    
    def getResult(self, name):
        """ get result of name 'name' from the result storage """
        return self.results.get(name)
    
    def printResult(self, name):
        """ print result of name 'name' from the result storage """
        print(self.results.get(name))
    
    #################################################################
    
    
    
    
    def callPreprocessingStep(self, method_name, column_name, *args, **kwargs):
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                method = getattr(csv, method_name, None)
                if method is None:
                    raise ValueError(f"CSVData does not have a method named {method_name}.")
                return method(column_name, *args, **kwargs), f"{column_name}_{method_name}"
        if not found:
            raise ValueError(f"Tried to perform NaN to zenith comparison for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    
    
    def checkForExtremeValues(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.compareNaNToZenith function, used for data consistency checking """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                return csv.checkForExtremeValues(column_name, *args, **kwargs), f"{column_name}_Extreme_Values"
        if not found:
            raise ValueError(f"Tried to perform NaN to zenith comparison for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def compareNaNToZenith(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.compareNaNToZenith function, used for data consistency checking """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                return csv.compareNaNToZenith(column_name, *args, **kwargs), f"{column_name}_NaN_to_Zenith"
        if not found:
            raise ValueError(f"Tried to perform NaN to zenith comparison for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def compareNonZeroToZenith(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.compareNonZeroToZenith function, used for data consistency checkingb """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                return csv.compareNonZeroToZenith(column_name, *args, **kwargs), f"{column_name}_Nonzero_to_Zenith"
        if not found:
            raise ValueError(f"Tried to perform nonzero to zenith comparison for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def interpolateNaN(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.interpolateNaN function for data preprocessing """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.interpolateNaN(column_name, *args, **kwargs)
                break
        if not found:
            raise ValueError(f"Tried to perform NaN interpolation for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
            
    def interpolateExtremeValues(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.interpolateExtremeValues function for data preprocessing """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.interpolateExtremeValues(column_name, *args, **kwargs)
                break
        if not found:
            raise ValueError(f"Tried to perform extreme value interpolation for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")

    def cleanNegatives(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.cleanNegatives function for data preprocessing """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.cleanNegatives(column_name, *args, **kwargs)
                break
        if not found:
            raise ValueError(f"Tried to perform negative cleaning for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def cleanNaN(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.cleanNaN function for data preprocessing """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.cleanNaN(column_name, *args, **kwargs)
                break
        if not found:
            raise ValueError(f"Tried to perform negative cleaning for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def cleanNaNAtNight(self, column_name, *args, **kwargs):
        """ Wrapper for CSVData.cleanNaNAtNight function for data preprocessing """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.cleanNaNAtNight(column_name, *args, **kwargs)
                break
        if not found:
            raise ValueError(f"Tried to perform NaN at night cleaning for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def cleanAllNaN(self, *args, **kwargs):
        """ Executes nan cleaning for all loaded csv's """
        for csv in self.csv_data:
            csv.cleanAllNaN(*args, **kwargs)
    
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
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    