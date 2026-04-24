from job_handler import JobHandler
from input_handler_modules import PandasCSVHandler

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



##############################################

def preprocessing(handler):
    handler.data.interpolateNaN("Pout")
    handler.data.interpolateNaN("G")
    handler.data.interpolateExtremeValues("T", value_limit=100)
    
    handler.data.cleanNonZeroToZenith("G", max_value=5.01, zenith_limit=100)
    
    handler.data.compareNaNToZenith("Pout", zenith_limit=80)
    handler.data.compareNaNToZenith("G", zenith_limit=80)
    handler.data.compareNonZeroToZenith("Pout", zenith_limit=100)
    handler.data.compareNonZeroToZenith("G", zenith_limit=100)
    handler.data.checkForExtremeValues("T", value_limit=100)
    
    handler.data.cleanAllNaN(new_value=0)
    handler.data.cleanNegatives("Pout")
    handler.data.cleanNegatives("G")
    return True, None


def main(handler, identifier=None):
    #Evaluate values and uncertainties of PR and temperature-corrected PR
    handler.evaluateVariable("PR")
    handler.evaluateVariable("PR_temp_corr")
    
    handler.calculateTotalUncertainty("PR", mask=True)
    handler.calculateTotalUncertainty("PR_temp_corr", mask=True)
    
    #You can define your own post-calculation validity checks if desired. 
    if handler.variables["PR"].values < 0.4:
        return False, "Unreliable_PR"
  
    if handler.variables["PR_temp_corr"].values < 0.4:
        return False, "Unreliable_PR_T"
 
    #Retrieve uncertainty contribution splits
    PR_u_split = handler.uncertainty_engine.calculateRootContributions(handler.variables["PR"])
    PR_source_names = handler.variables["PR"].uncertainty.getSourceNames()
    
    PR_T_u_split = handler.uncertainty_engine.calculateRootContributions(handler.variables["PR_temp_corr"])
    PR_T_source_names = handler.variables["PR_temp_corr"].uncertainty.getSourceNames()
    
    #Choose which results to store
    handler.store("PR values", "var.PR.values")
    handler.store("PR uncertainty", "var.PR.uncertainty.total_uncertainty")
    handler.store("PR u split", PR_u_split)
    handler.storeUniqueResult("PR u sources", PR_source_names)
    
    handler.store("PR T values", "var.PR_temp_corr.values")
    handler.store("PR T uncertainty", "var.PR_temp_corr.uncertainty.total_uncertainty")
    handler.store("PR T u split", PR_T_u_split)
    handler.storeUniqueResult("PR T u sources", PR_T_source_names)

    return True, None
   
##############################################
    
#Creating a dataframe from CSV data

CSV_filepath = "C:\\Users\\mate\\Desktop\\local_work\\code\\Uncertainty-Propagation\\Dataset-SolarTechLab.csv"
structure_list = ["Time", "Pout", "T", "-", "G", "W", "-"]
timeformat = None

#Create Pandas datahander instance, read a CSV, add a date column, add a solar zenith column
data_handler = PandasCSVHandler()
df = data_handler.readCSVData(CSV_filepath, ";", structure_list=structure_list, timeformat=timeformat, select_days=None)
del data_handler

##############################################

#Equation tree input file
equation_tree_filepath = "C:\\Users\\mate\\Desktop\\local_work\\code\\Uncertainty-Propagation\\test_tree.txt"

#Define coordinates and UTC offset of the location
coordinates = (45.30103, 9.092366)
UTC_offset = pd.Timedelta(1, 'hour')


#Create JobHandler instance, load equation tree, populate preprocessing and main functions
job = JobHandler()
job.loadEquationTree(equation_tree_filepath)

job.preprocessing = preprocessing
job.main          = main

#Add data, add a Date column, add a solar zenith column
job.addDataFrame(df)
job.data.addDateColumn(df_index=0)
job.data.addZenithColumn(coordinates, UTC_offset=UTC_offset, df_index=0)





#Execution over all data
job.execute(identifier="year", blacklist_mode="mask", results_group="year")
res = job.results.getResult(identifier="year", group="year")
print()
print("Dataset results:")
print(f"PR  : {np.round( res.data['PR values'], decimals=5)} +/- {np.round( res.data['PR uncertainty']*2, decimals=5)} (k=2)")
print(f"PR_T: {np.round( res.data['PR T values'], decimals=5)} +/- {np.round( res.data['PR T uncertainty']*2, decimals=5)} (k=2)")
print()





#Execution over daily data
days = job.data.getDays(name="Pout")
unique_days = days[:-1]
#Loop through the data day-by-day and execute the job
for i, day in enumerate(unique_days):
    #print(day)
    job.execute(day=day, identifier=day, blacklist_mode="fail")






print()
print("Daily results:")
print()

job.results.summariseFails()

#Retrieve results
PR_array, identifiers   = job.results.getResultArray("PR values", give_identifiers=True)
PR_u_array              = job.results.getResultArray("PR uncertainty")
PR_avg                  = job.results.getAverageResult("PR values")
PR_u_avg                = job.results.getAverageResult("PR uncertainty")
PR_u_split_avg          = job.results.getAverageResult("PR u split") * 100
PR_u_sources            = job.results.getUniqueResult("PR u sources")

PR_T_array       = job.results.getResultArray("PR T values")
PR_T_u_array     = job.results.getResultArray("PR T uncertainty")
PR_T_avg         = job.results.getAverageResult("PR T values")
PR_T_u_avg       = job.results.getAverageResult("PR T uncertainty")
PR_T_u_split_avg = job.results.getAverageResult("PR T u split") * 100
PR_T_u_sources   = job.results.getUniqueResult("PR T u sources")


job.results.summariseFails()

print(f"Avg daily PR                        : {PR_avg} +/- {PR_u_avg*2} (k=2)")
print(f"Avg daily PR (temperature corrected): {PR_T_avg} +/- {PR_T_u_avg*2} (k=2)")

print()
print("Source contribution splits [%]:")
print("PR")
for i, s in enumerate(PR_u_sources):
    print(f"{s}  {PR_u_split_avg[i]}")
print()
print("PR temperature corrected")
for i, s in enumerate(PR_T_u_sources):
    print(f"{s}  {PR_T_u_split_avg[i]}")


success_booleans, identifiers = job.results.getSuccessBooleans(as_array=True)
success_identifiers = identifiers[success_booleans]

plt.errorbar(x=success_identifiers, y=PR_array*100, yerr=2*PR_u_array*100, linestyle="", marker=".", label="PR")
plt.errorbar(x=success_identifiers, y=PR_T_array*100, yerr=2*PR_T_u_array*100, linestyle="", marker=".", label="PR (T25)")
plt.ylim(60,120)
plt.ylabel("PR [%]")
plt.xlabel("Date")
plt.legend()
plt.grid()
plt.show()





