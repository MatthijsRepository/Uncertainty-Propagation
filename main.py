from input_handler import InputHandler
from equation_engine import EquationEngine
from calculation_engine import CalculationEngine


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
    """ 
    
    
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
    
    #####
    #print(variables['time_c'].equation)
    #print(variables['time_c'].dependency_names)
    #print(variables['time_c'].dependencies.keys())
    
    #calculation_engine.harmonizeTimeSeries(variables['time_c'].dependencies)
    
    
    
    
    
    calculation_engine.calculateValues(variables["PR"], update_var=True)
    
    calculation_engine.calculateValues(variables["PR_temp_corr"], update_var=True)
    
    #print("Derived variables:")
    #for name in equation_engine.derived_variables:
    #    print(name)
    #print()
    
    

    
    
    
    
    