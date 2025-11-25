from job_handler import JobHandler
from input_handler_modules import PandasCSVHandler




#filepath = "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"
#structure_list = ["Time", "zenith", "G", "-", "Pout", "-"]
#timeformat = "%H:%M:%S"

filepath = "\\\\Office\\RedirectedFolders\\mate\\My Documents\\local files Matthijs\\Dataset-SolarTechLab.csv"
structure_list = ["Time", "Pout", "T", "-", "G", "W", "-"]
timeformat = None

inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\test_tree.txt"

job = JobHandler()
job.loadEquationTree(inputfile)




handler = PandasCSVHandler()
df = handler.readCSVData(filepath, ";", structure_list, timeformat=timeformat, select_days=None)
handler.addDateColumn(df)


import pandas as pd
coordinates = (45.30103, 9.092366)
UTC_offset = pd.Timedelta(1, 'hour')
handler.addZenithColumn(df, coordinates, UTC_offset)



job.addPreprocessingTask(job.interpolateNaN, "Pout")
job.addPreprocessingTask(job.interpolateNaN, "G")

job.addPreprocessingTask(job.compareNaNToZenith, "Pout", zenith_limit=82)
job.addPreprocessingTask(job.compareNaNToZenith, "G")
job.addPreprocessingTask(job.compareNonZeroToZenith, "Pout")
job.addPreprocessingTask(job.compareNonZeroToZenith, "G")

job.addPreprocessingTask(job.cleanAllNaN)
job.addPreprocessingTask(job.cleanNegatives, "Pout")
job.addPreprocessingTask(job.cleanNegatives, "G")




job.addTask(job.evaluateVariable, "var.PR")
job.addTask(job.evaluateVariable, "var.PR_temp_corr")
job.addTask(job.calculateTotalUncertainty, "var.PR")
job.addTask(job.calculateTotalUncertainty, "var.PR_temp_corr")
job.addTask(job.store, "PR values", "var.PR.values")
job.addTask(job.store, "PR T values", "var.PR_temp_corr.values")
job.addTask(job.store, "PR uncertainty", "var.PR.uncertainty.total_uncertainty")
job.addTask(job.store, "PR T uncertainty", "var.PR_temp_corr.uncertainty.total_uncertainty")



unique_days = df["Date"].unique()[:-1]

for day in unique_days:
    print(day)
    job.addCSVData(handler.compileOneDayCSVData(df, day))
    job.executeJob()
    



import sys
sys.exit()

import matplotlib.pyplot as plt
import time
t0 = time.time()
print("Started running the code")

#handler.addZenithColumn(df, coordinates, UTC_offset)
#print(f"Time elapsed calculating zenith angles: {time.time()-t0}")


unique_days = df["Date"].unique()[:-1]

Pout_skips = G_skips = T_skips = 0

for day in unique_days:
    print(day)
    data = handler.compileOneDayCSVData(df, day)
    
    if (not data.compareNaNToZenith("Pout", zenith_limit=82)): #or (not data.compareNonZeroToZenith("Pout")):
        
        import numpy as np
        
        data.cleanAllNaN(new_value=245)
        plt.plot((90-data.data["zenith"])/90)
        plt.plot(data.data["Pout"] / 245)
        plt.grid()
        plt.ylim(-0.8, 1.2)
        plt.title(day)
        plt.show()
        if Pout_skips>20:
            break
        
        Pout_skips += 1
        continue
    if (not data.compareNaNToZenith("G")): #or (not data.compareNonZeroToZenith("G")):
        G_skips += 1
        continue
    if not data.compareNaNToZenith("T"):
        T_skips += 1
        continue
    continue
    
    data.cleanAllNaN()
    data.cleanNegatives("Pout")
    data.cleanNegatives("G")
    
    job.populateVariablesFromCSV([data])
    job.executeJob()
    
    
    plt.plot((90-data.data["zenith"]) / 90)
    #plt.plot(data.data["Pout"] / max(data.data["Pout"]))
    #plt.plot(data.data["G"] / max(data.data["G"]))
    plt.plot(data.data["Pout"] / 245)
    plt.plot(data.data["G"] / 1000)
    
    plt.ylim(-0.8,1.2)
    plt.title(day)
    plt.grid()
    plt.show()
    

print(f"Pout skips: {Pout_skips}, G skips: {G_skips}, T skips: {T_skips}")

import matplotlib.pyplot as plt

plt.plot(job.getResult("PR values"))
plt.plot(job.getResult("PR T values"))
plt.grid()
plt.show()


plt.plot(job.getResult("PR uncertainty"))
plt.plot(job.getResult("PR T uncertainty"))
plt.grid()
plt.show()


job.printResult("PR values")
job.printResult("PR uncertainty")
job.printResult("PR T values")
job.printResult("PR T uncertainty")



print("Finished running, total time:")
print(time.time()-t0)











import sys
sys.exit()


inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\equation_tree.txt"

CSV_metadata_1 = (";", True, ["Time", "zenith", "G", "-", "Pout", "-"], "%H:%M:%S")
CSV_metadata_2 = (";", True, ["Time", "-", "-", "T", "-", "-"], "%H:%M:%S")
metadata = [CSV_metadata_1, CSV_metadata_2]
filepaths = ["C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv", "C:\\Users\\mate\\Desktop\\python\\Experimental\\testdata.csv"]


job = JobHandler()
job.loadEquationTree(inputfile)



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


