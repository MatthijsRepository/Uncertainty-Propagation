from input_handler_modules import EquationTreeReader, CSVHandler, PandasCSVHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine
import numpy as np
import pandas as pd
from dataclasses import dataclass



@dataclass
class Task:
    job_handler: object
    func: callable
    args: tuple = ()
    
    def __str__(self):
        return (f"Task: {self.func.__name__}{self.args}")
    
    def execute(self):
        args = self.job_handler.resolve_args(self.args)
        return self.func(*args)

class Results:
    def __init__(self):
        self.data = {}
        self.averages_effective_lengths = {}

    def add(self, key: str, value):
        """ Append a new result for a given key, automatically creates a field for the key if it does not exist yet """
        if key not in self.data:
            self.data[key] = [value]
        else:
            self.data[key].append(value)   
        
    def addAsAverage(self, key: str, value):
        """ Append a new result to an average over all results in this field,   A -> (A*(N-1) + value)/N   for the N'th result"""
        if key not in self.data:
            self.data[key] = value
            self.averages_effective_lengths[key] = 1
        else:
            total = (self.data[key] * self.averages_effective_lengths[key] + value)
            self.averages_effective_lengths[key] += 1
            self.data[key] = total / (self.averages_effective_lengths[key])
        
    def get(self, key: str):
        """ Get result as itself or as an array if it is list-like """
        if np.isscalar(self.data[key]):
            return self.data[key]
        else:
            return np.array(self.data[key], dtype=float)


class JobHandler:
    def __init__(self):
        self.CSV_handler = CSVHandler()
        
        self.equation_engine = None
        self.calculation_engine = None
        self.uncertainty_engine = None
        self.time_engine = None
        
        self.variables = None
        self.derived_variables_names = None
        self.var_csv_pointers = None
        self.csv_data = []
    
        self.data_init_job = [] ###!!!
        self.job = []
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
    
    def executeJob(self):
        """ Executes staged preprocessing and job for loaded csv data """
        if not self.initialized_engines:
            self.prepareEngines()
        
        if not self.has_csv_data:
            print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
        
        #Perform the data preprocessing
        for task in self.data_init_job:
            out = task.execute()
            #out is None for data cleaning tasks, and a boolean for data consistency checks
            #if out is False, then the data consistency check failed and we break off our job execution
            if out is False:
                print(f"Preprocessing {task} failed, breaking off execution.")
                return
            
        
        #Populate variables using loaded CSV data, and subsequently dump the csv data
        self.populateVariablesFromCSV()
        self.csv_data = []
        self.has_csv_data = False
        
        #Validate basic variables
        self.validateBasicVariables()
        
        #Execute job
        for task in self.job:
            task.execute()
        return True
    
    def resetJob(self):
        """ Resets the job list to an empty list """
        self.job = []
        
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
    
    def addPreprocessingTask(self, func, *args):
        """ Adds task to preprocessing job list """
        self.data_init_job.append(Task(self, func, args))
    
    def addTask(self, func, *args):
        """ Adds task to job list """
        self.job.append(Task(self, func, args))
    
    
    def getResult(self, name):
        """ get result of name 'name' from the result storage """
        return self.results.get(name)
    
    def printResult(self, name):
        """ print result of name 'name' from the result storage """
        print(self.results.get(name))
    
    #################################################################
    
    def compareNaNToZenith(self, column_name, *args):
        """ Wrapper for CSVData.compareNaNToZenith function for preprocessing job """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                return csv.compareNaNToZenith(column_name, *args)
        if not found:
            raise ValueError(f"Tried to perform NaN to zenith comparison for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def compareNonZeroToZenith(self, column_name, *args):
        """ Wrapper for CSVData.compareNonZeroToZenith function for preprocessing job """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                return csv.compareNonZeroToZenith(column_name, *args)
        if not found:
            raise ValueError(f"Tried to perform nonzero to zenith comparison for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def interpolateNaN(self, column_name, *args):
        """ Wrapper for CSVData.interpolateNaN function for preprocessing job """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.interpolateNaN(column_name, *args)
                break
        if not found:
            raise ValueError(f"Tried to perform NaN interpolation for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")

    def cleanNegatives(self, column_name, *args):
        """ Wrapper for CSVData.cleanNegatives function for preprocessing job """
        found = False
        for csv in self.csv_data:
            if column_name in csv.data.keys():
                found = True
                csv.cleanNegatives(column_name, *args)
                break
        if not found:
            raise ValueError(f"Tried to perform negative cleaning for {column_name}, but {column_name} was not found as an entry in the loaded csv data.")
    
    def cleanAllNaN(self, *args):
        """ Executes nan cleaning for all loaded csv's """
        for csv in self.csv_data:
            csv.cleanAllNaN(*args)
    
    #################################################################
    
    def store(self, name, arg):
        """ job task to store attribute 'arg' under name 'name' each job call """
        self.results.add(name, arg)
    
    def storeAsAverage(self, name, arg):
        """ job task to store the average of 'arg' under name 'name' over all job calls """
        self.results.addAsAverage(name, arg)
    
    def evaluateVariable(self, var):
        """ Wrapper for the calculation engine function of the same name """
        #if not self.basic_variables_validated:
        #    self.validateBasicVariables()
        
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate non-existing variable '{var}'.")
        self.calculation_engine.evaluateVariable(var)
        
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
    
    def calculateTotalUncertainty(self, var):
        """ Wrapper for the uncertainty engine function of the same name """
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.calculateTotalUncertainty(var)
    

        
        

    


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
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    