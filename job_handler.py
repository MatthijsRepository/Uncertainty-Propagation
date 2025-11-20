from input_handler_modules import EquationTreeReader, CSVHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine
import numpy as np
from dataclasses import dataclass



@dataclass
class Task:
    job_handler: object
    func: callable
    args: tuple = ()
    
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
    
        self.job = []
        self.results = Results()
        
        ##computational control flow booleans
        self.initialized_eq_tree       = False
        self.initialized_engines       = False
        self.basic_variables_validated = False
        self.csv_variables_populated   = False
        return
    
    def LoadEquationTree(self, filepath):
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
        for var in self.variables.values():
            if not var.is_hardcoded:
                var.reset()
        
        self.basic_variables_validated = False
        if len(self.var_csv_pointers) != 0:
            self.csv_variables_populated   = False
    
    def readFromCSV(self, filepaths, metadata, reset_registry=True):
        if reset_registry:
            self.resetVariableRegistry()
        
        CSV_data = []
        #First we read out all CSV data and store it as CSVData objects
        for i, filepath in enumerate(filepaths):
            args = metadata[i]
            temp_csv_data = self.CSV_handler.compileCSVData(filepath, *args)
            CSV_data.append(temp_csv_data)
        
        #Then for each variable in the csv_pointers dictionary we attempt to couple the referenced column name to a column name of our processed csv's
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
    
    def validateBasicVariables(self):
        if self.variables is None:
            raise ValueError("Validation of basic variables failed: no existing variable registry found.")
        if not self.csv_variables_populated:
            print("WARNING: trying to perform calculations while no CSV data appears to be loaded. Crash may occur.")
            
        self.calculation_engine.validateBasicVariables(equation_engine=self.equation_engine, variables=self.variables)
        self.basic_variables_validated = True
        return
    
    def preJobInitialization(self):
        if not self.initialized_engines:
            self.prepareEngines()
        self.validateBasicVariables()
        #IDEA:
        #replace NaN, do datacleaning, etc.
        #ValidateBasicVariables
        #... other things
        return
    
    def executeJob(self):
        self.preJobInitialization()
        for task in self.job:
            task.execute()
    
    def resetJob(self):
        self.job = []
        
    def _resolve_arg(self, arg):
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
        if not isinstance(args, tuple):
            return self._resolve_arg(args)
        return tuple(self._resolve_arg(arg) for arg in args)
    
    
    def addTask(self, func, *args):
        #if callable(func) and not hasattr(func, "__self__"):
        #    try:
        #        #Try to bind the passed function to JobHandler
        #        func = func.__get__(self, self.__class__)
        #    except:
        #        #We simply append the function to our job as-is
        #        pass
        #print(args)
        #args = self._resolve_args(args)
        #print(args)
        self.job.append(Task(self, func, args))
    
    
    #################################################################
    
    def store(self, name, arg):
        self.results.add(name, arg)
    
    def storeAsAverage(self, name, arg):
        self.results.addAsAverage(name, arg)
        
    def getResult(self, name):
        return self.results.get(name)
    
    def printResult(self, name):
        print(self.results.get(name))
        
    
    def evaluateVariable(self, var):
        #if not self.basic_variables_validated:
        #    self.validateBasicVariables()
        
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate non-existing variable '{var}'.")
        self.calculation_engine.evaluateVariable(var)
        
    def evaluateAllVariables(self):
        #if not self.basic_variables_validated:
        #    self.validateBasicVariables()
        for var_name in self.derived_variables_names:
            self.calculation_engine.evaluateVariable(self.variables[var_name])
            
    def prepareAllDirectUncertainties(self):
        self.uncertainty_engine.prepareAllDirectUncertainties()
    
    def prepareDownTreeDirectUncertainties(self, var):
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.prepareDownTreeDirectUncertainties(var)
    
    def calculateTotalUncertainty(self, var):
        if isinstance(var, str):
            var = self.variables.get(var)
            if var is None:
                raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var}'.")
        self.uncertainty_engine.calculateTotalUncertainty(var)
    

        
        

    





if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\equation_tree.txt"
    
    CSV_metadata_1 = (";", True, ["Time", "zenith", "G", "-", "Pout", "-"], "%H:%M:%S")
    CSV_metadata_2 = (";", True, ["Time", "-", "-", "T", "-", "-"], "%H:%M:%S")
    metadata = [CSV_metadata_1, CSV_metadata_2]
    filepaths = ["C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv", "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"]
    
    
    
    job = JobHandler()
    job.LoadEquationTree(inputfile)

    
    
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
    
    
    
    
    
    
    CSV_1_metadata = (delimiter, has_header, structure_list, timeformat)
    CSV_2_metadata = (delimiter, has_header, structure_list, timeformat)
    
    CSV_metadata = [CSV_1_metadata, CSV_2_metadata]
    

    CSV_1_characteristic = ...
    CSV_2_chara...
    
    CSV_characteristics = [...]
    
    identifier_list = [....]
    
    for identifier in identifier_list:        
        clean job_handler variables
        
        for i, characteristic in enumerate(CSV_characteristics):
            
        
        
    
    
    ###########
    
    jobhandler.executeJob(self, clean_data = True):
        #if clean_data: self.clean_data()
        
        load csvdata to variables
        
        execute jobscript
        
        store results (part of jobscript?)
        
        if clean_data: self.clean_data()
        return
    
    ########
    
    
    """
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    