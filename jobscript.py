from job_handler import JobHandler
from input_handler_modules import PandasCSVHandler

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#Defining filepaths and CSV structure, and coordinates + UTC offset of the dataset
#When defining multiple CSVs to read the data from, make sure to use unique identifiers for all columns

#Equation tree input file
equation_tree_filepath = "C:\\Users\\mate\\Desktop\\local_work\\code\\Uncertainty-Propagation\\test_tree.txt"


CSV_filepath = "C:\\Users\\mate\\Desktop\\local_work\\code\\Uncertainty-Propagation\\Dataset-SolarTechLab.csv"
structure_list = ["Time", "Pout", "T", "-", "G", "W", "-"]
timeformat = None

coordinates = (45.30103, 9.092366)
UTC_offset = pd.Timedelta(1, 'hour')

##############################################


#Define the the data preprocessing steps to be carried out
def preprocessing(handler, identifier=None):
    day=identifier
    #Initial missing value interpolation
    handler.data.interpolateNaN("Pout", day=day)
    handler.data.interpolateNaN("G", day=day)
    handler.data.interpolateExtremeValues("T", value_limit=100, day=day)
    
    handler.data.cleanNonZeroToZenith("G", day=day, max_value=5.01, zenith_limit=100)


    #Consistency checks    
    #success, error_code = handler.data.compareNaNToZenith("Pout", day, zenith_limit=80)
    success = handler.data.compareNaNToZenith("Pout", day, zenith_limit=80)
    if not success:
        
        """
        Pout = handler.data.getColumn("Pout", day=day, as_array=True)
        zenith = handler.data.getColumn("zenith", day=day, df_index=0, as_array=True)
        fig, ax = plt.subplots()
        ax.plot(Pout/245)
        ax.plot((90-zenith)/90)
        ax.hlines((90-80)/90, xmin=0, xmax=len(zenith))
        plt.show()
        """
        
        error_code = "NaNtoZenith_Pout"
        #handler.cleanNaN("Pout", new_value=245)
        #data = handler.getCSVColumn("Pout")
        #zenith = handler.getCSVColumn("zenith") ###!!!
        
        #plt.hlines((90-zen_lim) / 90, 0, 1440)
        #plt.plot(data/245)
        #plt.plot( (90 - zenith) / 90)
        #plt.grid()
        #plt.ylim(-0.8, 1.2)
        #plt.show()
        return success, error_code
    
    #success, error_code = handler.data.compareNaNToZenith("G", day=day, zenith_limit=80)
    success = handler.data.compareNaNToZenith("G", day=day, zenith_limit=80)
    if not success:
        error_code = "NanToZenith_G"
        return success, error_code
    
    #success, error_code = handler.data.compareNonZeroToZenith("Pout", day=day, zenith_limit=100)
    success = handler.data.compareNonZeroToZenith("Pout", day=day, zenith_limit=100)
    if not success:
        error_code = "NonzeroToZenith_Pout"
        return success, error_code
    
    #success, error_code = handler.data.compareNonZeroToZenith("G", day=day, zenith_limit=100)
    success = handler.data.compareNonZeroToZenith("G", day=day, zenith_limit=100)
    if not success:
        error_code = "NonzeroToZenith_G"
        return success, error_code
    
    #success, error_code = handler.data.checkForExtremeValues("T", day=day, value_limit=100)
    success = handler.data.checkForExtremeValues("T", day=day, value_limit=100)
    if not success:
        error_code = "Extreme_T"
        return success, error_code    
        
    #Final data cleaning    
    handler.data.cleanAllNaN(new_value=0, day=day)
    handler.data.cleanNegatives("Pout", day=day)
    handler.data.cleanNegatives("G", day=day)
    return True, None

#Define the main calculation procedure
def main(handler, identifier=None):
    #Evaluate values and uncertainties of PR and temperature-corrected PR
    handler.evaluateVariable("PR")
    handler.evaluateVariable("PR_temp_corr")
    
    handler.calculateTotalUncertainty("PR", mask=True)
    handler.calculateTotalUncertainty("PR_temp_corr", mask=True)
    
    #You can define your own post-calculation validity checks if desired. 
    if handler.variables["PR"].values < 0.4:
        return False, "Unreliable_PR"
        #plt.plot(handler.variables["G"].values / 1000)
        #plt.plot(handler.variables["Pout"].values / 245)
        #plt.grid()
        #plt.show()
    """
    if handler.variables["PR"].uncertainty.total_uncertainty>0.20:
        return False, "Unreliable_PR_uncertainty"
        #handler.calculateTotalUncertainty("G", mask=True)
        #plt.plot(handler.variables["G"].values)
        #plt.show()
        #handler.uncertainty_engine.plotAbsoluteRootContributions(handler.variables["G"])
    """
        
    if handler.variables["PR_temp_corr"].values < 0.4:
        #plt.plot(handler.variables["G"].values / 1000)
        #plt.plot(handler.variables["Pout"].values / 245)
        #plt.grid()
        #plt.show()
        #plt.plot(handler.variables["T"].values)
        #plt.title("T")
        #plt.show()
        #plt.plot(handler.variables["T_mod"].values)
        #plt.title("T_mod")
        #plt.show()
        return False, "Unreliable_PR_T"
    """
    if handler.variables["PR"].uncertainty.total_uncertainty>0.20:
        #handler.calculateTotalUncertainty("G", mask=True)
        #plt.plot(handler.variables["G"].values)
        #plt.title("G")
        #plt.show()
        #plt.plot(handler.variables["T"].values)
        #plt.title("T")
        #plt.show()
        handler.uncertainty_engine.plotAbsoluteRootContributions(handler.variables["G"])
        return False, "Unreliable_PR_T_uncertainty"
    """
    
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
    


import time
t0 = time.time()

#Create JobHandler instance, load equation tree, populate preprocessing and main functions
job = JobHandler()
job.loadEquationTree(equation_tree_filepath)

job.preprocessing = preprocessing
job.main          = main

#Create Pandas datahander instance, read a CSV, add a date column, add a solar zenith column
data_handler = PandasCSVHandler()
df = data_handler.readCSVData(CSV_filepath, ";", structure_list=structure_list, timeformat=timeformat, select_days=None)
del data_handler
#data_handler.addDateColumn(df)
#data_handler.addZenithColumn(df, coordinates, UTC_offset=UTC_offset)




job.addDataFrame(df)
job.data.addDateColumn(df_index=0)
job.data.addZenithColumn(coordinates, UTC_offset=UTC_offset, df_index=0)


days = job.data.getDays(name="Pout")

prep_time = time.time()-t0


#import sys; sys.exit()
t1 = time.time()

#Loop through the data day-by-day and execute the job
unique_days = days[:-1]
for i, day in enumerate(unique_days):
    if i>80: break
    print(day)
    #job.populateVariablesFromCSV(day=day)
    job.execute(day=day, identifier=day)
print()

main_time = np.round(time.time()-t1, decimals=2)

execution_time = np.round(time.time()-t0, decimals=2)
print(f"Total execution time: {execution_time}")
print(f"Total preparation time: {np.round(prep_time, decimals=2)} s - {np.round(prep_time/execution_time*100, decimals=2)} % ")
print(f"Total preprocessing time: {np.round(job.t_prepro, decimals=2)} s - {np.round(job.t_prepro/execution_time*100, decimals=2)} % ")
print(f"Total main time: {np.round(job.t_main, decimals=2)} s - {np.round(job.t_main/execution_time*100, decimals=2)} % ")
print(f"Total population time: {np.round(job.t_treepop, decimals=2)} s - {np.round(job.t_treepop/execution_time*100, decimals=2)} % ")
print(f"Total result time: {np.round(job.t_runresult, decimals=2)} s - {np.round(job.t_runresult/execution_time*100, decimals=2)} % ")
print()
print(f"Main time external: {main_time} s")
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


#print(PR_array)
#print(PR_u_array)

print(f"Avg PR                        : {PR_avg} +/- {PR_u_avg*2} (k=2)")
print(f"Avg PR (temperature corrected): {PR_T_avg} +/- {PR_T_u_avg*2} (k=2)")


success_booleans, identifiers = job.results.getSuccessBooleans(as_array=True)
success_identifiers = identifiers[success_booleans]

plt.errorbar(x=success_identifiers, y=PR_array*100, yerr=2*PR_u_array*100, linestyle="", marker=".", label="PR")
plt.errorbar(x=success_identifiers, y=PR_T_array*100, yerr=2*PR_T_u_array*100, linestyle="", marker=".", label="PR (T25)")

#plt.errorbar(x=np.arange(len(PR_array)), y=PR_array*100, yerr=2*PR_u_array*100, linestyle="", marker=".", label="PR")
#plt.errorbar(x=np.arange(len(PR_T_array))+1/3, y=PR_T_array*100, yerr=2*PR_T_u_array*100, linestyle="", marker=".", label="PR (T25)")

plt.ylim(60,120)
plt.ylabel("PR [%]")
plt.xlabel("Succesful run no.")
plt.legend()
plt.grid()
plt.show()

print()
print("Source contribution splits [%]:")
print("PR")
for i, s in enumerate(PR_u_sources):
    print(f"{s}  {PR_u_split_avg[i]}")
print()
print("PR temperature corrected")
for i, s in enumerate(PR_T_u_sources):
    print(f"{s}  {PR_T_u_split_avg[i]}")








