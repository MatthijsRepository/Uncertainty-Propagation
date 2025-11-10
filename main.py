from input_handler import InputHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine


if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Experimental\\uncertainty_testinput.txt"
    
    #Handling input
    print("Started parsing input")
    input_handler = InputHandler()
    input_handler.parse(inputfile)
    variables = input_handler.variables
    del input_handler
    print("Ended parsing input \n")
    
    
    #""" 
    print("Hardcoding testing variables for time series")
    from my_dataclasses import Variable
    import numpy as np
    from datetime import datetime
    
    timeformat = "%H:%M:%S"
    first_time_a = datetime.strptime("12:00:00", timeformat)        #10-minute
    last_time_a = datetime.strptime("16:10:00", timeformat)
    first_time_b = datetime.strptime("13:00:00", timeformat)        #Hourly
    last_time_b = datetime.strptime("16:00:00", timeformat)
    
    variables["time_a"] = Variable("time_a", values=np.ones(26)*5, is_basic=True, aggregation_rule="sum", first_time=first_time_a, last_time=last_time_a)
    #variables["time_a"] = Variable("time_a", values=np.arange(26), is_basic=True, aggregation_rule="sum", first_time=first_time_a, last_time=last_time_a)
    variables["time_b"] = Variable("time_b", values=np.array([2,4,3,2]), is_basic=True, aggregation_rule="sum", first_time=first_time_b, last_time=last_time_b)
    variables["time_c"] = Variable("time_c", values=None, is_basic=False, aggregation_rule="sum", equation="'time_a' * 'time_b'")
    
    
    import copy
    variables['time_a'].addUncertaintySource(copy.deepcopy(variables['G'].uncertainty.direct_uncertainty_sources[2]))
    variables['time_b'].addUncertaintySource(copy.deepcopy(variables['G'].uncertainty.direct_uncertainty_sources[1]))

    print("End hardcoding testing variables")
    print()
    #""" 
    
    
    #Verifying equation tree consistency, building equation tree in SymPy
    print("Verifying root-consistency of equation tree")
    equation_engine = EquationEngine(variables)
    equation_engine.checkEquationTreeConsistency(silent=True)
    print("Equation tree is root-consistent \n")
    print("Populating equation tree dependencies \n")
    equation_engine.populateEquationTreeDependencies()
    
    equation_engine.populateEquationTreeTimeSumSettings()
    
    print("Creating variable equation executables \n")
    equation_engine.buildEquationTreeExecutables()
    
    time_engine = TimeEngine()
    calculation_engine = CalculationEngine(variables, time_engine=time_engine)
    calculation_engine.validateBasicVariables(equation_engine)
    
    #calculation_engine.harmonizeTimeSeries(variables['time_c'].dependencies)

    equation_engine.buildPartialDerivativeExecutables(variables['time_c'])
    calculation_engine.evaluateVariable(variables['time_c'], update_var=True, silent=True)
    calculation_engine.executeAllPartials(variables['time_c'], absolute_values=False, store_results=True, force_recalculation=False)


    #calculation_engine.decreaseTemporalResolution(variables['G'], new_timestep=60, update_var=True, smuggle_limit=0)
    
    
    derived_variables = equation_engine.splitBasicDerived(variables)[1]    
    for var_name in derived_variables:
        calculation_engine.evaluateVariable(variables[var_name], update_var=True, silent=True)
    
    print(f"Calculated performance ratio: {variables['PR'].values}")
    print(f"Calculated temperature corrected performance ratio: {variables['PR_temp_corr'].values}")
    
    print("\n \n \n")
    print("TESTING AREA")
    print()
    print("Calculating uncertainties... - TESTING FROM HERE")
    
    
   

    #second_harm_dat = time_engine.decreaseVariableTemporalResolution(variables['G'], new_timestep=3600, benchmark_time="12:59:00", update_var=False, smuggle_limit=0)
    #print("Hourly totals")
    #print(second_harm_dat.new_values * 60)
    
    uncertainty_engine = UncertaintyEngine(variables, equation_engine = equation_engine, calculation_engine = calculation_engine)
    
    
    

    
    uncertainty_engine.prepareAllDirectUncertainties(variables.values())
    
    
    uncertainty_engine.calculateTotalUncertainty(variables['time_c'], recurse=True)
    
    print(variables['time_c'].uncertainty.total_uncertainty)
    
    
    
    uncertainty_engine.calculateTotalUncertainty(variables['G'], recurse=True)
    
    harm_dat = time_engine.calculateTimeHarmonizationData(variables['G'], new_timestep=120, benchmark_time = "12:59:00")
    
    
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('G')"], recurse=True)
    print(variables["TS_('G')"].uncertainty.total_uncertainty)
    
    
    uncertainty_engine.calculateTotalUncertainty(variables['my_test'], recurse=True)
    
    
    #hourly_aggregation = uncertainty_engine.partialAggregation(variables['G'], harm_dat)
    
    

    

    

    print()


    
    print(variables['G'].uncertainty.total_uncertainty)
    
    print(variables['G'].uncertainty.total_uncertainty[800])
    print(variables['G'].uncertainty.total_uncertainty[750])
    print(variables['G'].uncertainty.total_uncertainty[700])
    print(variables['G'].uncertainty.total_uncertainty[650])
    print(variables['G'].uncertainty.total_uncertainty[600])


    print(np.shape(variables['G'].uncertainty.total_uncertainty))

    
    

    
    

    
    
    
    
    