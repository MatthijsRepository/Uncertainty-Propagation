from my_dataclasses import Variable, UncertaintySource, ParsedVariableData
import numpy as np




class EquationTreeReader:
    def __init__(self):
        self.variables = {}
        self.csv_pointers = {}
        return
    
    def parse(self, filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith("#"):  # skip empty/comment lines
                i += 1
                continue
            if line.lower().startswith("var") or line.lower().startswith("eq"):
                i, variable_data = self._parseVariable(lines, i)
                new_var = self.createVariable(variable_data)
                self.variables[new_var.name] = new_var
                if variable_data.csv_pointer is not None:
                    self.csv_pointers[new_var.name] = variable_data.csv_pointer
            i += 1  # Move to next after handler finishes
        return self.variables, self.csv_pointers

    
    def _parseVariable(self, lines, i):
        """ Parses a single variable block in the input """
        variable_data = ParsedVariableData()
        
        #Check whether we are dealing with a basic or derived variable
        if lines[i].strip().startswith("var"):
            variable_data.var_name = lines[i].strip()[4:]
            variable_data.is_basic_variable = True
        else: #eq: "
            variable_data.var_name = lines[i].strip()[3:]
            variable_data.is_basic_variable = False

        i += 1
        while True:
            line = lines[i].strip()
            
            if not line or line.startswith("#"):
                i += 1 ; continue
            
            elif line.startswith("value"):
                variable_data.value_string = line.split(":")[1].strip()
                i += 1 ; continue
            elif line.startswith("timestep"):
                variable_data.timestep = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("description"):
                variable_data.description = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("rate"):
                variable_data.is_rate = line.split(":", maxsplit=1)[1].strip()
                if variable_data.is_rate.lower() == "true":
                    variable_data.is_rate = True
                i += 1 ; continue
            elif line.startswith("aggregation"):
                variable_data.aggregation_rule = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.startswith("equation"):
                variable_data.equation = line.split(":", maxsplit=1)[1].strip()
                i += 1 ; continue
            elif line.lower().startswith("uncertainties"):
                i = self._uncertaintyParser(variable_data, lines, i+1)
                i += 1 ; continue
            
            elif line.lower().startswith("end "):
                break
            elif i>len(lines):
                raise ValueError("VariableHandler failed to quit. Please inspect.")
            else:
                print("Input line in section Variable not recognized, skipping. Input line:")
                print(f"{line}")
                i += 1 ; continue
        return i, variable_data
    
    def createVariable(self, variable_data):       
        """ Dispatches the internal functions that handle the parsed input strings, creates a new variable and populates it """
        self._handleValueData(variable_data)
        
        new_variable = Variable(name             = variable_data.var_name, 
                                description      = variable_data.description, 
                                values           = variable_data.var_values, 
                                is_basic         = variable_data.is_basic_variable,
                                is_hardcoded     = variable_data.is_hardcoded,
                                is_rate          = variable_data.is_rate, 
                                aggregation_rule = variable_data.aggregation_rule,
                                equation         = variable_data.equation)
        if len(variable_data.uncertainties)>0:
            new_variable.addUncertaintySource(variable_data.uncertainties)
        return new_variable
    
    def _handleValueData(self, variable_data):
        """ Handles value data string """
        if variable_data.value_string is None:
            return
        if variable_data.value_string.lower().startswith("csv"):
            variable_data.var_values = None
            variable_data.csv_pointer = variable_data.value_string        
            return
        elif variable_data.value_string.startswith("DETERMINE") or variable_data.value_string.startswith("NO"):
            variable_data.var_values = None 
            return
        else:
            variable_data.is_hardcoded = True
            variable_data.var_values = float(variable_data.value_string) 
            return
    
    def _uncertaintyParser(self, variable_data, lines, i):
        """ Parses uncertainties block and dispatches uncertainty handler """
        while True:
            line = lines[i].strip()
            if not line or line.startswith("#"):
                i += 1 ; continue
            elif line.startswith("-"):
                self._handleUncertainty(variable_data, line)
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
    
    def _handleUncertainty(self, variable_data, line):
        """ Handles a single line in an uncertainty block """
        #Check inclusion
        if line.lower().endswith("include=false"):
            return None

        line = line[1:].split(":") #Remove dash in front, split name and data
        name = line[0].strip()
        data = line[1].split(",") #Split data over commas
        
        is_relative = shape = is_symmetric = correlation = bound = sigma = multiplier = None
        for i in range(len(data)):
            temp_data = data[i].strip().lower()
            #Relative or absolute
            if temp_data.startswith("relative"):
                is_relative=True ; continue
            elif temp_data.startswith("absolute"):
                is_relative=False ; continue
            #Shape
            elif temp_data.startswith("shape"):
                shape = temp_data.split("=")[1].strip().lower() ; continue
            #Symmetry
            elif temp_data.startswith("symme"):
                is_symmetric=True ; continue
            elif temp_data.startswith("one-"):
                is_symmetric=False ; continue
            #Correlation
            elif temp_data.startswith("corre"):
                correlation = float(temp_data.split("=")[1].strip()) ; continue
            #Bound and sigma
            elif temp_data.startswith("bound"):
                bound = float(temp_data.split("=")[1].strip())
                continue
            elif temp_data.startswith("sigma"):
                sigma = float(temp_data.split("=")[1].strip())
                continue
            #multiplier
            elif temp_data.startswith("multiplier"):
                multiplier = data[i].split("=")[1].strip() ; continue
            elif temp_data.startswith("include=false"):
                return None
            elif temp_data.startswith("include=true"):
                continue
            else:
                print(f"Uncertainty input field not recognized: {temp_data}")
                continue
        new_uncertainty =  UncertaintySource(name, is_relative, is_symmetric, shape, correlation, sigma=sigma, bound=bound, multiplier=multiplier)
        variable_data.uncertainties.append(new_uncertainty)
    
    
  
    
    
                
        
        
        
