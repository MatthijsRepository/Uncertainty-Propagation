from input_handler_modules import EquationTreeReader, CSVHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine
import numpy as np




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
        
        ##computational control flow booleans
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
        

    
    def prepareEngines(self):
        self.time_engine = TimeEngine()
        self.calculation_engine = CalculationEngine(variables          = self.variables, 
                                                    time_engine        = self.time_engine, 
                                                    equation_engine    = self.equation_engine)
        self.uncertainty_engine = UncertaintyEngine(variables          = self.variables,
                                                    equation_engine    = self.equation_engine, 
                                                    calculation_engine = self.calculation_engine, 
                                                    time_engine        = self.time_engine)
    
        
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
            print(args)
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
    
    def evaluateVariable(self, var_name):
        if not self.basic_variables_validated:
            self.validateBasicVariables()
        
        var = self.variables.get(var_name)
        if var is None:
            raise ValueError(f"Tried to evaluate non-existing variable '{var_name}'.")
            
        self.calculation_engine.evaluateVariable(var)
        
    def evaluateAllVariables(self):
        if not self.basic_variables_validated:
            self.validateBasicVariables()
        
        for var_name in self.derived_variables_names:
            self.calculation_engine.evaluateVariable(self.variables[var_name])
            
    def prepareAllDirectUncertainties(self):
        self.uncertainty_engine.prepareAllDirectUncertainties()
    
    def prepareDownTreeDirectUncertainties(self, var_name):
        var = self.variables.get(var_name)
        if var is None:
            raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var_name}'.")
        self.uncertainty_engine.prepareDownTreeDirectUncertainties(var)
    
    def calculateTotalUncertainty(self, var_name):
        var = self.variables.get(var_name)
        if var is None:
            raise ValueError(f"Tried to evaluate uncertainty for non-existing variable '{var_name}'.")
        self.uncertainty_engine.calculateTotalUncertainty(self, var_name)
        
        

    





if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\equation_tree.txt"
    
    CSV_metadata_1 = (";", True, ["Time", "zenith", "G", "-", "Pout", "-"], "%H:%M:%S")
    CSV_metadata_2 = (";", True, ["-", "-", "-", "T", "-", "-"], "%H:%M:%S")
    metadata = [CSV_metadata_1, CSV_metadata_2]
    filepaths = ["C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv", "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"]
    
    
    job = JobHandler()
    
    
    job.LoadEquationTree(inputfile)

    
    job.resetVariableRegistry()
        
    
    job.readFromCSV(filepaths, metadata)
    
    
    
    
    """ 
    
    CSV_1_metadata = (delimiter, has_header, structure_list, timeformat)
    CSV_2_metadata = (delimiter, has_header, structure_list, timeformat)
    
    CSV_metadata = [CSV_1_metadata, CSV_2_metadata]
    

    CSV_1_characteristic = ...
    CSV_2_chara...
    
    CSV_characteristic = [...]
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    