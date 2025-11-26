from job_handler import JobHandler
from input_handler_modules import PandasCSVHandler

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



#filepath = "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"
#structure_list = ["Time", "zenith", "G", "-", "Pout", "-"]
#timeformat = "%H:%M:%S"

filepath = "\\\\Office\\RedirectedFolders\\mate\\My Documents\\local files Matthijs\\Dataset-SolarTechLab.csv"
structure_list = ["Time", "Pout", "T", "-", "G", "W", "-"]
timeformat = None

inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\test_tree.txt"


def preprocessing(handler):
    handler.interpolateNaN("Pout")
    handler.interpolateNaN("G")
    
    if not handler.compareNaNToZenith("Pout", zenith_limit=80):
        
        #Convert nan of Pout to a specific value
        #handler.cleanNaN("Pout", new_value=245)
        #data = handler.getCSVColumn("Pout")
        #zenith = handler.getCSVColumn("zenith") ###!!!
        
        #plt.hlines((90-zen_lim) / 90, 0, 1440)
        #plt.plot(data/245)
        #plt.plot( (90 - zenith) / 90)
        #plt.grid()
        #plt.ylim(-0.8, 1.2)
        #plt.show()
        
        return False, "Pout_NaNToZenith"
    if not handler.compareNaNToZenith("G"):
        return False, "G_NaNToZenith"
    if not handler.compareNonZeroToZenith("Pout"):
        return False, "Pout_NonZeroToZenith"
    if not handler.compareNonZeroToZenith("G"):
        return False, "G_NonZeroToZenith"
    
    handler.cleanAllNaN()
    handler.cleanNegatives("Pout")
    handler.cleanNegatives("G")
        
    return True, None

def main(handler):
    handler.evaluateVariable("PR")
    handler.evaluateVariable("PR_temp_corr")
    
    handler.calculateTotalUncertainty("PR")
    handler.calculateTotalUncertainty("PR_temp_corr")
    
    #if handler.variables["PR"].uncertainty.total_uncertainty>0.15:
    #    print(f"Detected PR uncertainty of {handler.variables['PR'].uncertainty.total_uncertainty} (k=2)")
    #    uncertainty_split = handler.uncertainty_engine.calculateRootContributions(handler.variables["PR"])
    #    print("Uncertainty split: ")
    #    for i,source in enumerate(handler.variables["PR"].uncertainty.root_sources):
    #        print(f"{source.name}  :  {uncertainty_split[i]}")
    
    
    handler.store("PR values", "var.PR.values")
    handler.store("PR T values", "var.PR_temp_corr.values")
    handler.store("PR uncertainty", "var.PR.uncertainty.total_uncertainty")
    handler.store("PR T uncertainty", "var.PR_temp_corr.uncertainty.total_uncertainty")
    


job = JobHandler()
job.loadEquationTree(inputfile)

job.preprocessing = preprocessing
job.main          = main


data_handler = PandasCSVHandler()
df = data_handler.readCSVData(filepath, ";", structure_list, timeformat=timeformat, select_days=None)
data_handler.addDateColumn(df)


coordinates = (45.30103, 9.092366)
UTC_offset = pd.Timedelta(1, 'hour')
data_handler.addZenithColumn(df, coordinates, UTC_offset)


#Loop through the data day-by-day
unique_days = df["Date"].unique()[:-1]
i=0
for day in unique_days:
    #i+=1
    #if i>80: break
    print(day)
    job.addCSVData(data_handler.compileOneDayCSVData(df, day))
    job.execute(identifier=day)
print()
    

    




PR_array, PR_identifiers = job.results.getResultArray("PR values")
PR_u_array, PR_u_identifiers = job.results.getResultArray("PR uncertainty")


PR_T_array, PR_T_identifiers = job.results.getResultArray("PR T values")
PR_T_u_array, PR_T_u_identifiers = job.results.getResultArray("PR T uncertainty")


job.results.summariseFails()

print(PR_array)
print(PR_u_array)

plt.errorbar(x=np.arange(len(PR_array))*3, y=PR_array, yerr=2*PR_u_array, linestyle="", marker=".")
plt.errorbar(x=np.arange(len(PR_T_array))*3+1, y=PR_T_array, yerr=2*PR_T_u_array, linestyle="", marker=".")
plt.ylim(0.7,1.15)
plt.grid()
plt.show()





