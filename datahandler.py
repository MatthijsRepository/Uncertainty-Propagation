import numpy as np
import pandas as pd


    
class GroupInfo:
    def __init__(self, group, start_time, end_time, timestep):
        self.group      = group
        self.start_time = start_time
        self.end_time   = end_time
        self.timestep   = timestep

class DataHandler:
    def __init__(self):
        self.dataframes     = {}  #Stores the dataframes and potentially their groupings as DataSet objects
        self.groups         = {}  #Stores grouping information (typically by date) of the dataframes for easy by-day access, keys identical to those of dataframes
        self.lookup_dict    = {}  #Dictionary coupling variable names to a specific dataframe
        self.blacklist      = {}  #Dictionary populated by preprocessing functions, blacklisting days for certain reasons
        
    def addDataFrame(self, df, index_by_date=True):
        """ Adds dataframe to internal registry """
        index = len(self.dataframes)
        df.df_index = index
        
        #Check if the column names do not already exist in the dataset
        for name in df.columns:
            if name.lower() in ["time", "date", "zenith"]:
                continue
            if name in self.lookup_dict.keys():
                raise ValueError(f"Data column of name {name} is doubly defined: in dataframe {self.lookup_dict[name]} and dataframe {index}.")
            self.lookup_dict[name] = index
            
        self.dataframes[index] = df
        
    def getDataFrame(self, name=None, coupled_name=None, df_index=None, return_index=False):
        """ Get dataframe from registry
            df_index > coupled_name > name  """
        #Get the index of the data in self.dataframes
        if df_index is None:
            if coupled_name is not None:
                df_index = self.lookup_dict.get(coupled_name)
            else:
                df_index = self.lookup_dict.get(name)
        if df_index is None:
            raise ValueError(f"Failed to retrieve dataframe, no dataframe contains column {name} or {coupled_name}. Be aware that you cannot retrieve dataframes on columns named 'time' or 'date' or 'zenith'")
        
        #If return index: return the index of the dataframe in self.dataframes, Else: return only the df
        if return_index:
            return df_index
        else:
            return self.dataframes[df_index]
           
    def setColumn(self, name, values, base_values=None, day=None, coupled_name=None, df_index=None, df=None):
        """ Set values of a column """
        #Get dataframe
        if df is None:
            df = self.getDataFrame(name=name, coupled_name=coupled_name, df_index=df_index)
        
        #Create column if it does not exist yet
        if name not in df.columns:
            self.createColumn(name, base_values=base_values, df=df)
        
        #Set values
        if day is None:
            df[name] = values
        else:
            self.ensureDateColumn(df=df)
            df.loc[df["Date"] == day, name] = values
            
    def getColumn(self, name, day=None, coupled_name=None, df_index=None, df=None, as_array=False, as_copy=False, blacklist=[]):
        """ Get data of name, potentially for selected days. Columns like 'time' or 'zenith' may be degenerate, so these can be retrieved using
            'coupled_name', then the function returns the times of zenith angles corresponding to this column """
        #Get correct dataframe
        df_index = self.getDataFrame(name=name, coupled_name=coupled_name, df_index=df_index, return_index=True)
        df = self.dataframes[df_index]

        #Check if column exists
        if not name in df.columns:
            raise ValueError(f"Data retrieval failed: column {name} not present in dataframe: {df.columns}")
        
        if day is not None:
            self.ensureGroupedByDate(df_index)
            
            group_data = self.groups[df_index][day]
            data       = df.loc[group_data.group, name].to_numpy()
            if day in blacklist:
                data   = np.zeros(len(data))
            
            start_time = group_data.start_time
            end_time   = group_data.end_time
            timestep   = group_data.timestep
        else:
            #Create mask for the blacklisted days
            mask = df["Date"].isin(blacklist)
            #Ensure no data is destroyed, mask values of a copy
            data = df[name].copy()
            data[mask] = 0
            data = data.to_numpy()
            
            times      = df["Time"]
            start_time = times.iloc[0]
            end_time   = times.iloc[-1]
            timestep   = pd.Timedelta.total_seconds(times.iloc[1] - times.iloc[0])
        return data, start_time, end_time, timestep
        
    def setColumn(self, name, values, base_values=None, day=None, coupled_name=None, df_index=None, df=None):
        """ Set values of a column """
        #Get dataframe
        if df is None:
            df = self.getDataFrame(name=name, coupled_name=coupled_name, df_index=df_index)
        
        #Create column if it does not exist yet
        if name not in df.columns:
            self.createColumn(name, base_values=base_values, df=df)
        
        #Set values
        if day is None:
            df[name] = values
        else:
            self.ensureDateColumn(df=df)
            df.loc[df["Date"] == day, name] = values
            
    def createColumn(self, name, base_values=None, coupled_name=None, df_index=None, df=None):
        """ Creates a column 'name' and sets values to 'base_values' """
        #Check if the name is already in the lookup dictionary
        if self.lookup_dict.get(name) is not None:
            raise ValueError(f"Column of name {name} already present in dataframe: {self.lookup_dict.get(name)}.")
        #Retrieve dataframe
        if df is None:
            df = self.getDataFrame(name=None, coupled_name=coupled_name, df_index=df_index)
        if base_values is None:
            base_values = np.nan
        if name in df.columns:
            raise ValueError(f"Column of name {name} is already present in this dataframe: {df.columns}.")
        
        self.lookup_dict[name] = df.df_index
        df[name] = base_values
        
        
    ############################
    ### Routines for adding, ensuring or getting metadata
    ############################
    
    def compileBlacklist(self): ###!!!
        temp = {}
        for code, days in self.blacklist.items():
            for day in days:
                entry = temp.get(day)
                if entry is None:
                    temp[day] = []
                temp[day].append(code)
        return temp
    
    def getTimeRange(self, name=None, day=None, df_index=None, df=None):
        """ Retrieve the timestamps of the first and last datapoints of the selected window """
        if df is None:
            df = self.getDataFrame(name=name, df_index=df_index)
        if day is None:
            times = df["Time"]
        else:
            times = df.loc[df["Date"] == day, "Time"]
        return [times.iloc[0], times.iloc[-1]]
        
    
    def getDays(self, name=None, df_index=None, df=None):
        """ Get all unique days in a timeseries """
        if df is None:
            df = self.getDataFrame(name=name, df_index=df_index)
        self.ensureDateColumn(df)
        return df["Date"].unique()
    
    
    def ensureGroupedByDate(self, df_index):
        """ Ensures a dataframe is grouped by date """
        if self.groups.get(df_index) is not None:
            return
        else:
            self.groupByDate(df_index)
    
    def groupByDate(self, df_index): ###!!!
        df = self.getDataFrame(df_index=df_index)
        self.ensureDateColumn(df)
        
        results = {}
        
        for date, g in df.groupby("Date").groups.items():
            times = df.loc[g, "Time"]
            if len(times) < 2:
                continue
            start_time = times.iloc[0]
            end_time   = times.iloc[-1]
            timestep   = pd.Timedelta.total_seconds(times.iloc[1] - times.iloc[0])
            results[date] = GroupInfo(g, start_time, end_time, timestep)
        
        self.groups[df_index] = results
        
    def ensureDateColumn(self, df=None, df_index=None):
        """ Ensures a dataframe has a date column for indexing by days """
        if df is None:
            df = self.getDataFrame(df_index=df_index)
        if "Date" not in df.columns:
            self.addDateColumn(df=df)
        
    def addDateColumn(self, df=None, df_index=None):
        """ Adds a date column, extracted from the datetime column. For easier subsetting by date. """
        if df is None:
            df = self.getDataFrame(df_index=df_index)
        df["Date"] = pd.to_datetime(df["Time"]).dt.date ###!!!
    
    def addDateToTime(self, df=None, df_index=None):
        """ If a csv has a date and a time column, uses the date column to inform the times of their date """
        if df is None:
            df = self.getDataFrame(df_index=df_index)
        df["Time"] = pd.Series([
            pd.Timestamp.combine(d.date(), t.time()) for d, t in zip(df["Date"], df["Time"])
        ])
        return df
    
    def addZenithColumn(self, coordinates, df=None, df_index=None, time_zone=None, UTC_offset=None):
        """ Adds the solar zenith angle for each timestep with the given coordinates. 
        User needs to provide the UTC offset of the local time, since PVLib always works in UTC """
        from solar_module import calculateZenithAngles
        if df is None:
            df = self.getDataFrame(df_index=df_index)
        solar_data = calculateZenithAngles(coordinates, df["Time"], time_zone=time_zone, UTC_offset=UTC_offset)
        df["zenith"] = np.array(solar_data["zenith"])
        
    
    
    ############################
    ### Routines for filtering and cleaning data
    ############################   

    def deleteNaTAtEnds(self, df, col_name = "Time"):
        """ Deletes any Not a Time rows from the start and end of the datasets, if present.
            Use this function with caution: deleting data will break any precomputed groupings by date, these must be recomputed """
        valid = df[col_name].notna()
        
        if not valid.any():
            return df
        
        first_valid = valid.idxmax()
        last_valid = valid[::-1].idxmax()
        
        return df.loc[first_valid:last_valid]

    def checkForValidValues(self, column_name, day=None): ###!!! needs revision
        """ Checks if the column contains any defined value except for 0 
            Returns True if there are any values inside the dataset"""
        df, date_mask = self.getColumnView(column_name, day=day)
        mask = np.isnan(df.loc[date_mask, column_name]) | (df.loc[date_mask, column_name]==0)
        return not np.all(mask)
        
    def checkForExtremeValues(self, column_name, value_limit, day=None):
        """ Checks if there is any case where the data assumes a value (in absolute terms) greater than the value limit, returns true if the data does not contain extreme values """        
        df = self.getDataFrame(column_name)
        self.ensureDateColumn(df=df)
        
        mask = df[column_name].abs() > value_limit
        self.blacklist[f"Extreme_{column_name}"] = list(df.loc[mask, "Date"].unique())
    
    
    def checkForNaNInBody(self, column_name, day, body_start_after=0): ###!!!
        """ Checks if there are any NaN's or zeros in the body of a dataset. The body is defined to be the interval between the first non-zero value and the last non-zero value.
            Allows user to let the interval bounds to start after a certain amount of nonzero values have passed at both ends, to allow for alternating NaN and nonzero values at dawn and dusk. 
            Returns True if there are not any nan's or zeroes inside the body of the data. """
        column = self.getColumn(column_name, day=day)
        #True means element is nan or 0
        mask = np.isnan(column) | (column==0)
        indices = np.nonzero(~mask)[0]
        if len(indices) <= (2*body_start_after)+1:
            return True
        first_index, last_index = indices[body_start_after], indices[::-1][body_start_after]
        
        return not np.any(mask[first_index:last_index])
    
    def compareNaNToZenith(self, column_name, day=None, zenith_limit=85): 
        """ Checks if there are any nan's or zeroes in a column after the solar zenith angle is above a certain height
            Used to check whether a dataset contains undefined values during the day 
            Returns True if there are not any nan's or zeroes when the zenith is under the zenith limit (i.e. sun is high) """
        df = self.getDataFrame(name=column_name)
        if not "zenith" in df.columns:
            raise ValueError(f"Tried to compare {column_name} to solar zenith angle, while no zenith angle is defined for this csv.")
        
        self.ensureDateColumn(df=df)
        
        mask = (df[column_name].isna() | (df[column_name]==0)) & (df["zenith"] < zenith_limit)
        self.blacklist[f"Nan_zenith_{column_name}"] = list(df.loc[mask, "Date"].unique())
    
    def cleanNonZeroToZenith(self, column_name, day=None, max_value=1.01, zenith_limit=100):
        """ Checks if the column is nonzero after the solar zenith angle is below a certain height
            If this is the case and the nonzero values are below max_value, the values are replaced by zeroes
        """
        df = self.getDataFrame(name=column_name)
        if not "zenith" in df.columns:
            raise ValueError(f"Tried to compare {column_name} to solar zenith angle, while no zenith angle is defined for this csv.")
        
        mask = (~(df[column_name].isna() | df[column_name]==0) 
                & (df[column_name] < max_value) 
                & (df["zenith"] > zenith_limit) )
        df.loc[mask, column_name] = 0
    
    def compareNonZeroToZenith(self, column_name, day=None, zenith_limit=100):
        """ Checks if the column is nonzero after the solar zenith angle is below a certain height
            Used to check whether a column is nonzero during the night. 
            Returns true if there are no nonzero values when the zenith is over the zenith limit (i.e. sun is low). """
        df = self.getDataFrame(name=column_name)
        if not "zenith" in df.columns:
            raise ValueError(f"Tried to compare {column_name} to solar zenith angle, while no zenith angle is defined for this csv.")
        
        self.ensureDateColumn(df=df)
        
        mask = ~(df[column_name].isna() | (df[column_name]==0)) & (df["zenith"] > zenith_limit)
        self.blacklist[f"Nonzero_zenith_{column_name}"] = list(df.loc[mask, "Date"].unique())

        
    def interpolateNaN(self, column_name, day=None, min_value=0, else_value=np.nan):
        """ Interpolates isolated nan's, that are neighboured by two defined values, as the average of their neighbours.
            In case the average is less than min_value, the nan is instead replaced by else_value.
            This is to allow nan's near the beginning and end of the day to be set to 0 instead of to the interpolated value. """
        #Obtain nan's that is neighboured by two defined values
        df = self.getDataFrame(name=column_name)
        
        isnan       = df[column_name].isna()
        left_valid  = ~isnan.shift(1, fill_value=True)
        right_valid = ~isnan.shift(-1, fill_value=True)
        isolated    = isnan & left_valid & right_valid

        interpolated_values = ( df[column_name].shift(1) + df[column_name].shift(-1) ) / 2
        
        df.loc[isolated & (interpolated_values > min_value), column_name] = interpolated_values
        df.loc[isolated & (interpolated_values < min_value), column_name] = else_value
    
    def interpolateExtremeValues(self, column_name, value_limit):
        """ Interpolates isolated values that exceed (in absolute terms) the extreme_limit. Used to remove unphysical values from dataset """
        df = self.getDataFrame(column_name)
        
        extremes    = df[column_name].abs() > value_limit
        left_valid  = ~extremes.shift(1, fill_value=True)
        right_valid = ~extremes.shift(-1, fill_value=True)
        isolated    = extremes & left_valid & right_valid
        
        interpolated_values = ( df[column_name].shift(1) + df[column_name].shift(-1) ) / 2
        
        df.loc[isolated, column_name] = interpolated_values
        
    def cleanExtremeValues(self, column_name, value_limit=10000, new_value=np.nan):
        """ Replaces extreme values by new value """
        df = self.getDataFrame(column_name)
        df.loc[df[column_name].abs() > value_limit, column_name] = new_value
        
    
    def cleanNegatives(self, column_name, new_value=0):
        """ Replaces negative values by specified value """
        df = self.getDataFrame(column_name)
        df.loc[df[column_name] < 0, column_name] = new_value

        
    def cleanNaN(self, column_name, new_value):
        """ Cleans NaN instances and replaces them by a new value. Replacement value can be manually specified """
        df = self.getDataFrame(column_name)
        df.loc[df[column_name].isna(), column_name] = new_value
    
    def cleanNaNAtNight(self, column_name, day): ###!!!
        """ Cleans NaN's before and after the main day data to assume the first and last value of the data that is not NaN, respectively. """
        column = self.getColumn(column_name, day=day, as_array=True)
        #Identifiying start and end indices
        indices = ~np.isnan(column)
        indices = np.where(indices)[0]
        first, last = indices[0], indices[-1]

        #Replacing values
        column[:first] = column[first]
        column[last:] = column[last]
        self.setColumn(column_name, column, day=day)
        
    def cleanAllNaN(self, new_value=0):
        """ Replaces all NaN's for the given time window by the new value """
        for df_index, df in self.dataframes.items():
            df.fillna(new_value, inplace=True)
    
    
































        