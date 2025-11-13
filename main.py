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
    #variables['time_a'].addUncertaintySource(copy.deepcopy(variables['G'].uncertainty.direct_uncertainty_sources[2]))
    #variables['time_b'].addUncertaintySource(copy.deepcopy(variables['G'].uncertainty.direct_uncertainty_sources[1]))

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

    
    uncertainty_engine = UncertaintyEngine(variables, equation_engine = equation_engine, calculation_engine = calculation_engine)
    
    
    

    
    uncertainty_engine.prepareAllDirectUncertainties(variables.values())
    
    """ 
    uncertainty_engine.calculateTotalUncertainty(variables['P_error_test'], recurse=True)
    rel_error = np.divide(variables["P_error_test"].uncertainty.total_uncertainty,
                          variables["P_error_test"].values,
                          out = np.zeros_like(variables["P_error_test"].values),
                          where=(variables["P_error_test"].values) != 0)
    print(variables['Pout'].values[500])
    print(variables['Pout'].uncertainty.direct_uncertainty_sources[0].values[500])
    print(variables["P_error_test"].uncertainty.total_uncertainty[500])
    print(variables["P_error_test"].values[500])
    print(rel_error[500])
    
    
    uncertainty_engine.calculateTotalUncertainty(variables['P_error_agg_test'], recurse=True)
    print(variables['P_error_agg_test'].uncertainty.total_uncertainty*60 / variables['P_error_agg_test'].values)
    print(variables['P_error_agg_test'].uncertainty.total_uncertainty)
    print(variables['P_error_agg_test'].values * 60)
    
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('Pout' / 'P0')"], recurse=True)
    print(variables["TS_('Pout' / 'P0')"].values)
    print(variables["TS_('Pout' / 'P0')"].uncertainty.total_uncertainty)
    print(variables["TS_('Pout' / 'P0')"].uncertainty.total_uncertainty*60 / variables["TS_('Pout' / 'P0')"].values)
    import sys
    sys.exit()
    """ 
    
    uncertainty_engine.calculateTotalUncertainty(variables['G'], recurse=True, mask=True)
    uncertainty_engine.calculateTotalUncertainty(variables['PR'], recurse=True, mask=True)
    uncertainty_engine.calculateTotalUncertainty(variables['PR_temp_corr'], recurse=True, mask=True)
    
    
    uncertainty_engine.plotRootContributions(variables['PR'])
    uncertainty_engine.plotRootContributions(variables['G'])
    uncertainty_engine.plotAbsoluteRootContributions(variables['G'], k=2)
    uncertainty_engine.plotRelativeRootContributions(variables['G'], k=2)
    
    print("PRINTING OUTPUTS")
    
    print(f"Calculated performance ratio: {variables['PR'].values} +/- {2*variables['PR'].uncertainty.total_uncertainty[0]} (k=2)")
    print(f"Calculated temperature corrected performance ratio: {variables['PR_temp_corr'].values}  +/- {2*variables['PR_temp_corr'].uncertainty.total_uncertainty[0]} (k=2)")
    
    
    variables['G'].giveReport()
    variables['PR'].giveReport(decimals=5)
    variables['PR_temp_corr'].giveReport(decimals=5)
    
    
    
    uncertainty_engine.calculateTotalUncertainty(variables['Pout'], recurse=True)
    rel_error = np.divide(variables['Pout'].uncertainty.total_uncertainty,
                      variables['Pout'].values,
                      out = np.zeros_like(variables['Pout'].values),
                      where=(variables['Pout'].values != 0))
    print(f"Relative Pout error: {rel_error[500]}")
    
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('Pout' / 'P0')"], recurse=True)
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('G' / 'G_STC')"], recurse=True)
    
    
    variables["TS_('Pout' / 'P0')"].giveReport(short_report=True)
    variables["TS_('G' / 'G_STC')"].giveReport(short_report=True)
    #print(variables["TS_('Pout' / 'P0')"].values)
    #print(variables["TS_('G' / 'G_STC')"].values)
    
    
    
    
    print(variables['PR'].uncertainty.getSource("directional response"))
    
    
    