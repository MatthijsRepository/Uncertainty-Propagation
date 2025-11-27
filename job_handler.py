from input_handler_modules import EquationTreeReader, CSVHandler, PandasCSVHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Union, Optional


@dataclass
class RunResult:
    identifier: Union(str, int)  #String which stores the identifier of the run
    succeeded: bool              #Stores whether the run succeeded or not
    data: dict                   #Stores results of the run
    fail_code: Optional(str)    #Stores which datacheck caused preprocessing to fail
    
    
class Results:
    def __init__(self):
        self.num_runs           = 0             #Stores the number of runs this result object contains
        self.run_identifiers    = []            #For each run, can be used to an identifier (such as the date)
        self.run_results        = []            #Stores RunResult object for each run
        #self.averages_effective_lengts = {}     #Stores effective lengths N for data that is stored as average: A -> (A*(N-1) + value)/N for the N'th result
        self.staged_data        = {}           #Staged data dictionary to be populated in the present run
        
    def add(self, key: str, value):
        """ Append a new result for a given key, automatically creates a field for the key if it does not exist yet """
        if key not in self.staged_data:
            self.staged_data[key] = value
            return
        
        column = self.staged_data[key]
        if not isinstance(column, list):
            column = list(column)
        column.append(value)
                
    
    def createRunResult(self, succeeded=True, identifier=None, fail_code=None):
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
        
    def getResultSeries(self, name):
        series, identifiers = [], []
        for result in self.run_results:
            datapoint = result.data.get(name)
            if datapoint is None:
                continue
            series.append(datapoint)
            identifiers.append(result.identifier)
        return series, identifiers
    
    def getResultArray(self, name):
        series, identifiers = self.getResultSeries(name)
        return np.asarray(series), identifiers
    
    def getFails(self, failcode=None):
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
        return fails, identifiers
    
    def summariseFails(self):
        fails, identifiers = self.getFails()
        successful_executions = self.num_runs - len(fails)
        
        unique_fails = list(set(fails))
        print("Summarising calculation failures:")
        print(f"Successful executions: {successful_executions} out of {self.num_runs} total executions")
        for fail in unique_fails:
            print(f"Error code {fail} occurred {fails.count(fail)} times")
        print()
            
            


class JobHandler:
    def __init__(self):
        self.CSV_handler = CSVHandler()
        
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
        return
    
    def loadEquationTree(self, filepath):
        self.equation_tree_reader = EquationTreeReader()
        self.variables, self.var_csv_pointers = self.equation_tree_reader.parse(filepath)
        del self.equation_tree_reader
        
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
    
    
    def getCSVColumn(self, name):
        if not self.has_csv_data:
            raise ValueError(f"Tried to extract column {name} from job handler CSV data, but no CSV data is loaded in handler.")
        for csv in self.csv_data:
            if name in csv.data.keys():
                return csv.data[name]
            raise ValueError(f"Tried to extract column {name} from job handler CSV data, but loaded csv's do not contain {name}.")
    
    def populateVariablesFromCSV(self, reset_registry=True):
        """ For each variable in the csv_pointers dictionary this function attempts to couple the referenced column name to a column name of our processed CSVs """
        if reset_registry:
            self.resetVariableRegistry()
            
        for var_name, column_name in self.var_csv_pointers.items():
            #Column name is of the form "CSV.name", we strip the first 4 characters
            column_name = column_name[4:]
            
            found = False
            for csv in self.csv_data:
                if column_name in csv.data.keys():
                    found = True
                    self.variables[var_name].values = np.asarray(csv.data[column_name])
                    if csv.time_range is not None:
                        self.variables[var_name].addTimeStep(csv.time_range)
                    break
            if not found:
                if not self.has_csv_data:
                    raise ValueError(f"Failed to populate variable {var_name}: no csv data loaded to jobhandler.")
                else:
                    raise ValueError(f"CSV data does not contain data named {column_name}.")
        self.csv_variables_populated = True
            

    def populateVariablesFromCSV_OLD(self, CSV_data, reset_registry=True):
        """ For each variable in the csv_pointers dictionary this function attempts to couple the referenced column name to a column name of our processed CSVs """
        if reset_registry:
            self.resetVariableRegistry()
        
        if not isinstance(CSV_data, (list, np.ndarray)):
            CSV_data = [CSV_data]
        
        for var_name, column_name in self.var_csv_pointers.items():
            #Column name is of the form "CSV.name", we strip the first 4 characters
            column_name = column_name[4:]
            
            found = False
            for csv in CSV_data:
                if column_name in csv.data.keys():
                    found = True
                    self.variables[var_name].values = np.asarray(csv.data[column_name])
                    if csv.time_range is not None:
                        self.variables[var_name].addTimeStep(csv.time_range)
                    break
            if not found:
                raise ValueError(f"CSV data does not contain data named {column_name}.")
        self.csv_variables_populated = True

    
    def readFromCSV(self, filepaths, metadata, clean_nan=True, reset_registry=True):
        """ Populate variables from given csv's using the functions in the CSVHandler function
            Note: this function does not use pandas """
        if reset_registry:
            self.resetVariableRegistry()
        
        CSV_data = []
        #First we read out all CSV data and store it as CSVData objects
        for i, filepath in enumerate(filepaths):
            args = metadata[i]
            temp_csv_data = self.CSV_handler.compileCSVData(filepath, *args)
            #Optional cleaning of NaN
            if clean_nan:
                temp_csv_data.cleanAllNaN()
            CSV_data.append(temp_csv_data)
        
        #populate variables, set flag to True
        self.populateVariablesFromCSV_OLD(CSV_data)
    
    def validateBasicVariables(self):
        """ Wrapper for calculation engine function of the same name, also updates the relevant flag """
        #if self.variables is None:
        #    raise ValueError("Validation of basic variables failed: no existing variable registry found.")
        #if not self.csv_variables_populated:
        #    print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
            
        self.calculation_engine.validateBasicVariables(equation_engine=self.equation_engine, variables=self.variables)
        self.basic_variables_validated = True
        return
    
    
    def execute(self, identifier=None):
        if self.main is None:
            raise ValueError("No main jobscript is provided to the job handler. Please provide a main function under JobHandler.main")
        if not self.initialized_engines:
            self.prepareEngines()
        
        if not self.has_csv_data:
            print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
        
        #Perform the data preprocessing
        if self.preprocessing is not None:
            success, fail_code = self.preprocessing(self)
            #If preprocessing failed, we log this in the results, we also clear the loaded csv data
            if not success:
                self.results.createRunResult(succeeded=False, identifier=identifier, fail_code=fail_code)
                self.csv_data = []
                self.has_csv_data = False
                return
        
        #Populate variables using loaded CSV data, and subsequently dump the csv data
        self.populateVariablesFromCSV()
        self.csv_data = []
        self.has_csv_data = False
        
        #Validate basic variables
        self.validateBasicVariables()
        
        #Execute job
        success, fail_code = self.main(self)
        if not success:
            self.results.createRunResult(succeeded=False, identifier=identifier, fail_code=fail_code)
            self.csv_data = []
            self.has_csv_data = False
            return
        
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
    
    def cleanAllNaN(self, *args, **kwargs):
        """ Executes nan cleaning for all loaded csv's """
        for csv in self.csv_data:
            csv.cleanAllNaN(*args, **kwargs)
    
    #################################################################
    
    def store(self, name, arg):
        """ job task to store attribute 'arg' under name 'name' each job call """
        arg = self._resolve_arg(arg)
        self.results.add(name, arg)
    
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
    

        
        

    


"""


if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\equation_tree.txt"
    
    CSV_metadata_1 = (";", True, ["Time", "zenith", "G", "-", "Pout", "-"], "%H:%M:%S")
    CSV_metadata_2 = (";", True, ["Time", "-", "-", "T", "-", "-"], "%H:%M:%S")
    metadata = [CSV_metadata_1, CSV_metadata_2]
    filepaths = ["C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv", "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"]
    
    
    job = JobHandler()
    job.loadEquationTree(inputfile)
    
    job.variables
    
    
    job.addTask(job.evaluateVariable, "var.PR")
    job.addTask(job.evaluateVariable, "var.PR_temp_corr")
    
    job.addTask(job.calculateTotalUncertainty, "var.PR")
    job.addTask(job.calculateTotalUncertainty, "var.PR_temp_corr")
    
    job.addTask(job.store, "PR values", "var.PR.values")
    job.addTask(job.storeAsAverage, "PR temp corr values", "var.PR_temp_corr.values")
    
    job.addTask(job.store, "PR uncertainty split", "var.PR.uncertainty.aggregated_weighted_uncertainties")
    job.addTask(job.store, "PR temp corr uncertainty split", "var.PR_temp_corr.uncertainty.aggregated_weighted_uncertainties")
    
    job.addTask(job.storeAsAverage, "PR uncertainty split average", "var.PR.uncertainty.aggregated_weighted_uncertainties")
    
    
    
    
    
    for _ in range(3):
        job.readFromCSV(filepaths, metadata)
        job.executeJob()
        
    print("PR values")
    job.printResult("PR values")
    
    #print()
    #print("PR uncertainty split")
    #job.printResult("PR uncertainty split")
    
    print()
    print("PR uncertainty split average")
    job.printResult("PR uncertainty split average")
    """
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    