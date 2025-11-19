from input_handler import InputHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine
from time_engine import TimeEngine


if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Uncertainty-Propagation\\uncertainty_testinput.txt"
    
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
    variables["time_b"] = Variable("time_b", values=np.array([3,2,3,4]), is_basic=True, aggregation_rule="sum", first_time=first_time_b, last_time=last_time_b)
    variables["time_c"] = Variable("time_c", values=None, is_basic=False, aggregation_rule="sum", equation="'time_a' * 'time_b'")
    
    #variables['time_a'].addUncertaintySource(variables['G'].uncertainty.direct_uncertainty_sources[2])
    #variables['time_b'].addUncertaintySource(variables['G'].uncertainty.direct_uncertainty_sources[1])

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
    print(variables['time_c'].partial_values['time_b'])

    
    #calculation_engine.decreaseTemporalResolution(variables['G'], new_timestep=120, update_var=True, smuggle_limit=30)
    
    
    derived_variables = equation_engine.splitBasicDerived(variables)[1]    
    for var_name in derived_variables:
        calculation_engine.evaluateVariable(variables[var_name], update_var=True, silent=True)
    
    print(f"Calculated performance ratio: {variables['PR'].values}")
    print(f"Calculated temperature corrected performance ratio: {variables['PR_temp_corr'].values}")
    
    print("\n \n \n")
    print("TESTING AREA")
    print()
    print("Calculating uncertainties... - TESTING FROM HERE")
    
    
    uncertainty_engine = UncertaintyEngine(variables, equation_engine = equation_engine, calculation_engine = calculation_engine)
    
    
    

    
    #uncertainty_engine.calculateUncertainty(variables['time_c'], recurse=True)
    #uncertainty_engine.calculateCorrelation(variables['time_c'], auto_calculate=True, recurse=True, force_recalculation=False)
    
    
    
    uncertainty_engine.calculateUncertainty(variables["G"], recurse=True)
    uncertainty_engine.splitDirectUncertaintyContributions(variables['G'])
    uncertainty_engine.splitTotalUncertaintyContributions(variables['G'])
    uncertainty_engine.splitToSourceContributions(variables['G'])
    uncertainty_engine.calculateCorrelation(variables['G'], auto_calculate=True, recurse=True, force_recalculation=False)
    
    
    two_minute_harm_dat = time_engine.calculateTimeHarmonizationData(variables['G'], new_timestep=120, benchmark_time = "12:59:00")
    u_new, c_new = uncertainty_engine.partialAggregation(variables['G'], two_minute_harm_dat)
    print(u_new[300])
    import sys
    sys.exit()
    
    
    
    second_harm_dat = time_engine.decreaseVariableTemporalResolution(variables['G'], new_timestep=3600, benchmark_time="12:00:00", update_var=False, smuggle_limit=0)
    print(second_harm_dat.new_values * 60)
    
    #Aggregation test
    harm_dat = time_engine.calculateTimeHarmonizationData(variables['G'], new_timestep=3600, benchmark_time = "12:59:00")
    u_new, c_new = uncertainty_engine.partialAggregation(variables['G'], harm_dat)
    print("Hourly aggregation G")
    print(u_new*60 / 1000 * 2)
    print(c_new)
    print()
    c_new = np.repeat(c_new[:,np.newaxis],axis=1, repeats=24)
    c_new = (c_new + c_new.T) / 2
    np.fill_diagonal(c_new, 1)
    
    #u_new[:4] = 0
    #u_new[-4:] = 0
    print(u_new*60 / 1000 * 2)
    
    tot = np.sqrt(u_new.T @ c_new @ u_new)
    tot = tot * 60 * 2 / 1000 / 1000
    print(tot)
    print()
    tot_2 = uncertainty_engine.timeAggregateTotalUncertainty(variables['G'])
    print(tot_2 * 60 * 2 / 1000 / 1000)
    print("END section")
    

    
    #variables['G'].uncertainty.plotRelativeRootSplit(k=2)
    #variables['G'].uncertainty.plotAbsoluteRootSplit(k=2)
        

    
    from datetime import datetime
    timeformat = "%H:%M:%S"
    
    
    test_time = datetime.strptime("12:59:00", timeformat)
    test_harmonization_data = time_engine.decreaseVariableTemporalResolution(variables['G'], new_timestep=3600, benchmark_time=test_time, update_var=False, smuggle_limit=0)
    
    
    
    uncertainty_engine.partialAggregation(variables['G'], test_harmonization_data)
    
    #"""

    start_time_a = datetime.strptime("12:04:15", timeformat)
    #start_time_a = datetime.strptime("12:59:01", timeformat)
    
    #variables['G'].values = variables['G'].values[1:]
    #variables['G'].start_time = datetime.strptime("00:00:30", timeformat)
    #variables['G'].first_time = datetime.strptime("00:01:00", timeformat)
    #import numpy as np
    #variables['G'].values = np.ones(len(variables['G'].values))
    #calculation_engine.decreaseVariableTemporalResolution(variables['G'], 3600, benchmark_time=start_time_a, smuggle_limit=60)
    
    #"""
    print(variables['PR'].values)
    
    uncertainty_engine.calculateUncertainty(variables['PR'], recurse=True)
    uncertainty_engine.calculateCorrelation(variables['PR'], auto_calculate=True)
    #u = uncertainty_engine.timeAggregateTotalUncertainty(variables['PR'])
    #print(u)
    #print(u / variables['PR'].values)
    print(variables['G'].uncertainty.total_uncertainty[800])
    print(variables['G'].uncertainty.correlation[800])
    print(variables["TS_('G')"].uncertainty.total_uncertainty * 60 * 2 / 1000)
    print(variables["TS_('Pout')"].uncertainty.total_uncertainty / 1000)
    print(variables['PR'].uncertainty.total_uncertainty)
    print()
    
    
    uncertainty_engine.calculateUncertainty(variables['PR_temp_corr'], recurse=True)
    uncertainty_engine.calculateCorrelation(variables['PR_temp_corr'], auto_calculate=True)
    print(variables['PR_temp_corr'].uncertainty.total_uncertainty)
    
    
    


    
    

    
    
    
    
    