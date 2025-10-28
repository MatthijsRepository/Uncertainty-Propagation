from input_handler import InputHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine
from uncertainty_engine import UncertaintyEngine


if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Uncertainty-Propagation\\myinput.txt"
    
    #Handling input
    print("Started parsing input")
    input_handler = InputHandler()
    input_handler.parse(inputfile)
    variables = input_handler.variables
    del input_handler
    print("Ended parsing input \n")
    
    
    """ 
    print("Hardcoding testing variables for time series")
    from my_dataclasses import Variable
    import numpy as np
    from datetime import datetime
    
    timeformat = "%H:%M:%S"
    start_time_a = datetime.strptime("12:00:00", timeformat)
    end_time_a = datetime.strptime("16:10:00", timeformat)
    start_time_b = datetime.strptime("12:30:00", timeformat)
    end_time_b = datetime.strptime("15:30:00", timeformat)
    
    variables["time_a"] = Variable("time_a", values=np.ones(26), is_basic=True, aggregation_rule="sum", first_time=start_time_a, last_time=end_time_a)
    variables["time_b"] = Variable("time_b", values=np.array([1,2,3,4]), is_basic=True, aggregation_rule="sum", first_time=start_time_b, last_time=end_time_b)
    variables["time_c"] = Variable("time_c", values=None, is_basic=False, aggregation_rule="sum", equation="'time_a' + 'time_b'")
    print("End hardcoding testing variables")
    print()
    #""" 
    
    
    #Verifying equation tree consistency, building equation tree in SimPy
    print("Verifying root-consistency of equation tree")
    equation_engine = EquationEngine(variables)
    equation_engine.checkEquationTreeConsistency(silent=True)
    print("Equation tree is root-consistent \n")
    print("Populating equation tree dependencies \n")
    equation_engine.populateEquationTreeDependencies()
    
    
    print("Creating variable equation executables \n")
    equation_engine.buildEquationTreeExecutables()
    
    
    calculation_engine = CalculationEngine(variables)
    calculation_engine.validateBasicVariables(equation_engine)
    
    
    
    derived_variables = equation_engine.splitBasicDerived(variables)[1]
    for var_name in derived_variables:
        print(var_name)
        calculation_engine.calculateValues(variables[var_name], update_var=True, silent=True)
    
    
    uncertainty_engine = UncertaintyEngine(variables, equation_engine = equation_engine)
    
    uncertainty_engine.calculateUncertainty(variables["G"], recurse=True)
    uncertainty_engine.splitDirectUncertaintyContributions(variables['G'])
    uncertainty_engine.splitTotalUncertaintyContributions(variables['G'])
    uncertainty_engine.splitToSourceContributions(variables['G'])
    
    uncertainty_engine.calculateCorrelation(variables['G'], auto_calculate=True, recurse=True, recalculate=False)
    
    
    print(variables['G'].uncertainty.correlation)
    print(variables['G'].uncertainty.correlation[400])
    print(variables['G'])
    
    
    print("TESTING")
        
    
    #variables['G'].uncertainty.plotRelativeRootSplit(k=2)
    #variables['G'].uncertainty.plotAbsoluteRootSplit(k=2)
    
    from datetime import datetime
    timeformat = "%H:%M:%S"
    start_time_a = datetime.strptime("12:04:15", timeformat)
    #start_time_a = datetime.strptime("12:59:01", timeformat)
    
    #variables['G'].values = variables['G'].values[1:]
    #variables['G'].start_time = datetime.strptime("00:00:30", timeformat)
    variables['G'].first_time = datetime.strptime("00:01:00", timeformat)
    import numpy as np
    variables['G'].values = np.ones(len(variables['G'].values))
    calculation_engine.decreaseTemporalResolution(variables['G'], 3600, start_time_a, smuggle_limit=60)
    
    
    
    mytest = np.arange(10)
    start=-1
    factor=4
    print(mytest[start+1:start+factor-1])
    print(mytest[start+factor-1])


    """ 
    for i, name in enumerate(variables['G'].uncertainty.root_uncertainty_sources):
        #plt.plot(variables['G'].uncertainty.root_uncertainty_contribution_split[i], label=name)
        plt.plot(absolute_uncertainty_split[i], label=name)
        #plt.plot(relative_uncertainty_split[i], label=name)

    plt.legend()
    plt.show()
    """ 
    
    
    
    
    
    """
    #####
    #print(variables['time_c'].equation)
    #print(variables['time_c'].dependency_names)
    #print(variables['time_c'].dependencies.keys())
    
    #calculation_engine.harmonizeTimeSeries(variables['time_c'].dependencies)
    
    print("TEST TRIVIAL")
    calculation_engine.calculateValues(variables['test_trivial'])
    print(variables['test_trivial'].values)
    
    
    calculation_engine.calculateValues(variables['timesumtest'])
    
    print('HA')
    #var = variables["timesumtest"]
    
    var = variables["TS_('G' * TS_('C_25'))"]
    
    equation_engine.buildPartialDerivativeExecutables(var)
    var.executeAllPartials()
    

    #"""
    
    #calculation_engine.calculateValues(variables["PR"], update_var=True)
    
    #calculation_engine.calculateValues(variables["PR_temp_corr"], update_var=True)
    
    #print("Derived variables:")
    #for name in equation_engine.derived_variables:
    #    print(name)
    #print()
    
    

    
    
    
    
    