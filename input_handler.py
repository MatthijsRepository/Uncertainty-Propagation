import csv
from datetime import datetime
from my_dataclasses import CSVData, Uncertainty, Variable
import numpy as np



class InputHandler:
    def __init__(self):
        self.csv_data = None    #CSVData dataclass: Holds CSV datasets
        self.variables = {}     #dict: Holds variables
        
    def parse(self, filepath):
        """ Main parser that dispatches to sub-parsers """
        with open(filepath, "r") as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith("#"):  # skip empty/comment lines
                i += 1
                continue
            
            # Dispatch logic
            if line.lower().startswith("csv"):
                handler = CSVHandler()
                i = handler.parse(lines, i)
                self.csv_data = handler.compileCSVData()
                del handler
            elif line.lower().startswith("var") or line.lower().startswith("eq"):
                handler = VariableHandler(csv_data = self.csv_data)
                i = handler.parse(lines, i)
                new_variable = handler.createVariable()
                self.variables.update({new_variable.name: new_variable})
                del handler
            i += 1  # Move to next after handler finishes
                

class CSVHandler:
    def __init__(self):
        self.csv_name = None            #str: CSV name
        self.has_header = None          #bool: defines if csv has header
        self.has_timedata = None        #bool: defines whether CSV contains timedata
        self.data = None                #dict: data paired with variable names
        self.structure_list = None      #list: lists csv structure interpretation
        self.data = None                #dict: data paired with variable names
        self.timeformat = None          #str: time format
        self.file_path = None           #str: path to csv file
        
    def _postParseCheck(self):
        if self.delimiter is None:
            raise ValueError("Did not receive delimiter in input.")
        if self.structure_list is None:
            raise ValueError("Did not receive file structure in input.")
        if self.file_path is None:
            raise ValueError("Did not receive file location in input.")
    
    def _readCSV(self):
        """ Reads a CSV and populates the object data block """
        with open(self.file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=self.delimiter)
            if self.has_header:
                next(reader)
            for row in reader:
                for idx, name in enumerate(self.structure_list):
                    #Skip element if it should be excluded, flagged by '-'
                    if name != '-':
                        #Time section handling
                        if name.lower() == 'time':
                            value = datetime.strptime(row[idx], self.timeformat)
                        #All other data is attempted to be stored as a float
                        else:
                            try:
                                value = float(row[idx])
                            except ValueError:
                                #In case of value errors, paste data as string
                                value = row[idx]
                        self.data[name].append(value)
            
    def parse(self, lines, i):
        """ Parses single CSV input section and obtains variables """ 
        self.csv_name = lines[i].strip()[4:]
        i += 1
        
        while True:
            line = lines[i].strip()
            
            if not line or line.startswith("#"):
                i += 1 ; continue
            elif line.startswith("header"):
                if line.split(":")[1].strip().lower() == "true":
                    self.has_header = True
                i += 1 ; continue
            elif line.startswith("delimiter"):
                self.delimiter = line.split(':', maxsplit=1)[1].strip()[1]
                i += 1 ; continue
            elif line.startswith("structure"):
                line = line.split(':')[1]
                self.structure_list = [s.strip() for s in line.split(',')]
                self.data = {name: [] for name in self.structure_list if name != '-'}
                i += 1 ; continue
            elif line.startswith("timeformat"):
                self.timeformat = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("file"):
                self.file_path = line.split(":", maxsplit=1)[1].strip()
                #If filename is inside quotes: remove quotes
                if self.file_path.startswith('"') or self.file_path.startswith("'"):  
                    self.file_path = self.file_path[1:-1] 
                i += 1 ; continue
            elif line.lower().startswith("end csv"):
                break
            elif i>len(lines):
                raise ValueError("CSVHandler failed to quit. Please inspect.")
            else:
                print("Input line in section CSV not recognized, skipping. Input line:")
                print(f"{line}")
                i += 1 ; continue
            
        self._postParseCheck()
        return i  #If finished correctly, then the index now rests at an 'end csv' line
    
    def compileCSVData(self):
        """ Read CSV and populate self.data block """
        self._readCSV()
        
        #Check whether time was given as an input
        for name in self.structure_list:
            if name=="Time":                                        ###!!! Not robust: Time, time
                self.has_timedata = True
        
        #If csv contains timedata, get range and timestep
        if self.has_timedata:
            time_range = [self.data["Time"][0], self.data["Time"][-1]]      ###!!! Not robust: Time, time
            timestep = self.data["Time"][1] - self.data["Time"][0]          ###!!! Not robust: Time, time
            timestep = timestep.total_seconds()
        else:
            time_range = None ; timestep = None
        
        return CSVData(self.csv_name, self.data, timestep, time_range)
        
    
class VariableHandler:
    def __init__(self, csv_data=None):
        self.csv_data = csv_data        #list: csv data in inputhandler object
        self.value_string = None        #str: string containing input value definition
        self.var_values = None          #[int, float, array]: variable values
        self.timestep = None            #[float, str]: timestep of variable, or setting ###!!!
        self.description = None         #str: variable description
        self.is_basic_variable = None   #bool: flags if variable is basic or derived
        self.aggregation_rule = None    #str: aggregation rule of the quantity
        self.equation = None            #str: equation of the variable
        self.dependency_names = None    #list: lists all variables used in the equation
        self.uncertainties = []         #list: contains all uncertainty sources
        return  
        
    def parse(self, lines, i):
        """ Parses a single variable block in the input """
        #Check whether we are dealing with a basic or derived variable
        if lines[i].strip().startswith("var"):
            self.var_name = lines[i].strip()[4:]
            self.is_basic_variable = True
        else: #eq: "
            self.var_name = lines[i].strip()[3:]
            self.is_basic_variable = False

        i += 1
        while True:
            line = lines[i].strip()
            
            if not line or line.startswith("#"):
                i += 1 ; continue
            
            elif line.startswith("value"):
                self.value_string = line.split(":")[1].strip()
                i += 1 ; continue
            elif line.startswith("timestep"):
                self.timestep = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("description"):
                self.description = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("aggregation"):
                self.aggregation_rule = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("equation"):
                self.equation = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("variables"):
                self.dependency_names = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.lower().startswith("uncertainties"):
                i = self._uncertaintyParser(lines, i+1)
                i += 1 ; continue
            
            elif line.lower().startswith("end "):
                break
            elif i>len(lines):
                raise ValueError("VariableHandler failed to quit. Please inspect.")
            else:
                print("Input line in section Variable not recognized, skipping. Input line:")
                print(f"{line}")
                i += 1 ; continue
        return i
    
    def createVariable(self):       
        """ Dispatches the internal functions that handle the parsed input strings, creates a new variable and populates it """
        self._handleValueData()
        self._handleTimestep()
        
        new_variable = Variable(name=self.var_name, description=self.description, values=self.var_values, \
                                is_basic=self.is_basic_variable, aggregation_rule=self.aggregation_rule, equation=self.equation)
        if not self.timestep is None:
            new_variable.addTimeStep(self.time_range)
        if len(self.uncertainties)>0:
            new_variable.addUncertaintySource(self.uncertainties)
        return new_variable

    def _handleValueData(self):
        """ Handles value data string """
        if self.value_string is None:
            return None
        if self.value_string.lower().startswith("csv"):
            if self.csv_data is None:
                raise ValueError(f"CSV data referenced for variable {self.var_name}, but no CSV dataset was passed. Please input CSV data.")
            data_name = self.value_string.split('.')[1]
            try:
                self.var_values = np.asarray(self.csv_data.data[data_name])
            except:
                raise ValueError(f"CSV data does not contain data named {data_name}.")
        elif self.value_string.startswith("DETERMINE"):
            self.var_values = None
        elif self.value_string.startswith("NO"):
            self.var_values = None
        else:
            self.var_values = float(self.value_string)
    
    def _handleTimestep(self):
        """ Handles timestep settings, sets up timestep settings in object variables """
        if not self.timestep is None:
            if self.timestep.lower().startswith("none"):
                self.timestep = None
            elif self.timestep.lower().startswith("auto"):
                #self.timestep = (self.timestep.split("(")[1])[:-1] #Get values inside brackets
                #self.timestep = self.timestep.split(",") #List of start and end time strings
                #self.time_range = [self.timestep[0].strip()]
                raise ValueError("Automatic timestep definition currently not supported.")
            elif self.timestep.lower().startswith("csv"):
                if self.csv_data is None:
                    raise ValueError(f"CSV data referenced for variable {self.var_name}, but no CSV dataset was passed. Please input CSV data.")
                try:
                    self.time_range = self.csv_data.time_range #[self.csv_data.data["Time"][0], self.csv_data.data["Time"][-1]]      ###!!! Not robust: Time, time
                    #self.timestep = self.csv_data.timestep #"auto"
                except:
                    raise ValueError(f"CSV dataset provided for variable {self.var_name} does not contain time data.")
            else:
                raise ValueError(f"Timestep string {self.timestep} not recognized, please check input.")
    
    
    def _uncertaintyParser(self, lines, i):
        """ Parses uncertainties block and dispatches uncertainty handler """
        while True:
            line = lines[i].strip()
            if line.startswith("-"):
                self._handleUncertainty(line)
                i += 1 ; continue
            elif line.lower().startswith("end un"): #End uncertainties
                break
            elif i>len(lines):
                raise ValueError(f"Uncertainty handler failed to quit for variable {self.var_name}, please inspect.")
            else:
                print("Input line in section Variable not recognized, skipping. Input line:")
                print(f"{line}")
                i += 1 ; continue
        return i
    
    def _handleUncertainty(self, line):
        """ Handles a single line in an uncertainty block """
        #Check inclusion
        if line.lower().endswith("false"):
            return None

        line = line[1:].split(":") #Remove dash in front, split name and data
        name = line[0].strip()
        data = line[1].split(",") #Split data over commas
        #Absolute or relative
        if data[0].strip().lower().startswith("rel"):
            is_relative=True
        else:
            is_relative=False
        #Shape is passed directly
        shape = data[1].strip().lower()
        #Sidedness
        if data[2].strip().lower().startswith("sym"):
            is_symmetric = True
        else:
            is_symmetric = False
        #Correlation
        correlation = float(data[3].split("=")[1].strip())
        #Bound or sigma
        if data[4].strip().lower().startswith("sig"):
            bound = None
            sigma = float(data[4].split("=")[1].strip())
        else:
            bound = float(data[4].split("=")[1].strip())
            sigma = None
        
        new_uncertainty =  Uncertainty(name, is_relative, is_symmetric, shape, correlation, sigma=sigma, bound=bound)
        self.uncertainties.append(new_uncertainty)
        
    
    
                
        
        
        





