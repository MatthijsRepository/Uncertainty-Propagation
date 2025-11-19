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
        self.var_csv_pointers = None
        return
    
    def readEquationTree(self, filepath):
        self.equation_tree_reader = EquationTreeReader()
        self.variables, self.var_csv_pointers = self.equation_tree_reader.parse(filepath)
        del self.equation_tree_reader
        
    def initializeEquationTree(self):
        if self.variables is None:
            raise ValueError("Cannot initialize equation tree: variable registry is empty.")
        
        self.equation_engine = EquationEngine(self.variables)
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
      
    
    def populateVariablesFromCSV(self, filepaths, metadata, reset_registry=True):
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
    
    





if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\equation_tree.txt"
    
    CSV_metadata_1 = (";", True, ["Time", "zenith", "G", "-", "Pout", "-"], "%H:%M:%S")
    CSV_metadata_2 = (";", True, ["-", "-", "-", "T", "-", "-"], "%H:%M:%S")
    metadata = [CSV_metadata_1, CSV_metadata_2]
    filepaths = ["C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv", "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"]
    
    
    job_handler = JobHandler()
    
    
    job_handler.readEquationTree(inputfile)
    job_handler.initializeEquationTree()
    
    job_handler.resetVariableRegistry()
        
    
    job_handler.populateVariablesFromCSV(filepaths, metadata)
    
    
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    