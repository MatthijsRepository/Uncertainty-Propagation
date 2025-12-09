from my_dataclasses import Variable, UncertaintySource, ParsedVariableData, CSVData
import numpy as np
import csv
from datetime import datetime


import pandas as pd
class PandasCSVHandler:
    def __init__(self):
        return
    
    def readCSVData(self, filepath, delimiter, skip_rows=None, structure_list=None, structure_dict=None, timeformat=None, select_days=None):
        """ Reads CSV data into a pandas dataframe """
        
        #Skip header rows
        if skip_rows is None:
            skip_rows = 0
        
        #At least a structure list or dict must be provided
        if structure_list is None and structure_dict is None:
            raise ValueError("Please provide either a structure list or structure dictionary when reading out a CSV!")
        
        #Read out the columns and match to the desired names
        if structure_dict is not None:
            #Read CSV
            df = pd.read_csv(filepath, sep=delimiter, header=None, skiprows=skip_rows)
            #Extract the correct columns from the dataframe
            columns_to_read = list(structure_dict.keys())
            df = df.iloc[:, columns_to_read]
            #Rename columns
            df = df.rename(columns=structure_dict)
        else:
            #Usage of structure list is deprecated
            df = pd.read_csv(filepath, sep=delimiter)
            rename_map = dict(zip(list(df.columns), structure_list))
            df = df.rename(columns=rename_map)
            df.pop("-")
        
        # Parse datetime
        if timeformat is None:
            # Let pandas infer format
            df["Time"] = pd.to_datetime(df["Time"], errors='coerce')  ###!!! Not robust
        else:
            df["Time"] = pd.to_datetime(df["Time"], format=timeformat) ###!!!
        
        #If a specific "Date" column is present, convert its values to datetime objects and inform the "Time" column of its dates
        if "Date" in df.columns:
            df["Date"] = pd.to_datetime(df["Date"]).dt.date
            df["Time"] = pd.Series([
                pd.Timestamp.combine(d, t.time()) for d, t in zip(df["Date"], df["Time"])
            ])

        df.set_index("Time") ###!!!
        # Filter for selected days
        if select_days is not None:
            days = pd.to_datetime(select_days).normalize()
            df = df[df.index.normalize().isin(days)]
        return df
    
    def addDateColumn(self, df):
        """ Adds a date column, extracted from the datetime column. For easier subsetting by date. """
        df["Date"] = pd.to_datetime(df["Time"]).dt.date ###!!!
    
    def addDateToTime(self, df):
        """ If a csv has a date and a time column, uses the date column to inform the times of their date """
        df["Time"] = pd.Series([
            pd.Timestamp.combine(d.date(), t.time()) for d, t in zip(df["Date"], df["Time"])
        ])
        return df
    
    def addZenithColumn(self, df, coordinates, time_zone=None, UTC_offset=None):
        """ Adds the solar zenith angle for each timestep with the given coordinates. 
        User needs to provide the UTC offset of the local time, since PVLib always works in UTC """
        from solar_module import calculateZenithAngles
        solar_data = calculateZenithAngles(coordinates, df["Time"], time_zone=time_zone, UTC_offset=UTC_offset)
        df["zenith"] = np.array(solar_data["zenith"])
        
        
    def compileOneDayCSVData(self, df, date):
        """ Calls the compileCSVData function on the subset of the dataframe corresponding to the given date """
        subset = df[df["Date"] == date]
        return self.compileCSVData(subset)
        
    
    def compileCSVData(self, df):
        """ Compile a CSVdata object from the given (subset of a) dataframe """
        time_range = [df["Time"].iloc[0].to_pydatetime(), df["Time"].iloc[-1].to_pydatetime()] ###!!!
        timestep = (df["Time"].iloc[1].to_pydatetime() - df["Time"].iloc[0].to_pydatetime()) ###!!!
        timestep = timestep.total_seconds()
        
        data = {name: (df[name] if name in ["Time", "Date"] else df[name].to_numpy() ) for name in df.columns} ###!!!
        return CSVData(data, timestep, time_range)
        


class CSVHandler:
    """ Deprecated CSV handler, recommended to use the PandasCSVHandler """
    def _readCSV(self, filepath, delimiter, has_header, structure_list, timeformat=None):
        """ Reads a CSV and populates the object data block """
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=delimiter)
            if has_header:
                next(reader)
            for row in reader:
                for idx, name in enumerate(structure_list):
                    #Skip element if it should be excluded, flagged by '-'
                    if name != '-':
                        #Time section handling
                        if name.lower() == 'time':
                            value = datetime.strptime(row[idx], timeformat)
                        #All other data is attempted to be stored as a float
                        else:
                            try:
                                value = float(row[idx])
                            except ValueError:
                                #In case of value errors, paste data as string
                                value = row[idx]
                        self.data[name].append(value)
    
    def compileCSVData(self, filepath, delimiter, has_header, structure_list, timeformat=None):
        """ Read CSV and populate self.data block """
        self.data = {name: [] for name in structure_list if name != '-'}
        
        self._readCSV(filepath, delimiter, has_header, structure_list, timeformat)
        
        #Check whether time was given as an input
        has_timedata = False
        if "Time" in structure_list:
            has_timedata = True
        
        #If csv contains timedata, get range and timestep
        if has_timedata:
            time_range = [self.data["Time"][0], self.data["Time"][-1]]      ###!!! Not robust: Time, time
            timestep = self.data["Time"][1] - self.data["Time"][0]          ###!!! Not robust: Time, time
            timestep = timestep.total_seconds()
        else:
            time_range = None ; timestep = None
        
        return CSVData(self.data, timestep, time_range)


class EquationTreeReader:
    def __init__(self):
        self.variables = {}
        self.csv_pointers = {}
        return
    
    def parse(self, filepath):
        """ Main parses of the file """
        with open(filepath, "r") as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith("#"):  # skip empty/comment lines
                i += 1
                continue
            if line.lower().startswith("var") or line.lower().startswith("eq"):
                #Parse the variable section and populate the ParsedVariableData object
                i, variable_data = self._parseVariable(lines, i)
                #Create a new variable using the parsed data
                new_var = self.createVariable(variable_data)
                #Add new variable to variable registry
                self.variables[new_var.name] = new_var
                #Add pointer to csv column to the csv pointer registry
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
            elif line.startswith("maskable"):
                variable_data.is_maskable = line.split(":", maxsplit=1)[1].strip()
                if variable_data.is_maskable.lower() == "false":
                    variable_data.is_maskable = False
                else: variable_data.is_maskable = True
                i += 1 ; continue
            elif line.startswith("rate"):
                variable_data.is_rate = line.split(":", maxsplit=1)[1].strip()
                if variable_data.is_rate.lower() == "true":
                    variable_data.is_rate = True
                else: variable_data.is_rate = False
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
                                is_maskable      = variable_data.is_maskable,
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
    
    
  
    
    
                
        
        
        
