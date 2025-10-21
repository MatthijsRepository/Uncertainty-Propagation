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
    
    calculation_engine.calculateValues(variables["PR"], update_var=True)
    
    
    #print("Derived variables:")
    #for name in equation_engine.derived_variables:
    #    print(name)
    #print()
    
    #print(variables["testC"].dependency_names)
    #equation_engine.buildEquation(variables["testC"])
    
    #variables['testC'].executeEquation()
    #print(variables['testC'].values)
    
    #variables['C_25'].executeEquation()
    #print(variables['T'].values[0])
    #print(variables['C_25'].values[0])
    
    #print()
    #print(f"{variables['testD'].name} = {variables['testD'].equation} with dependencies: ")
    #print(f"{variables['testD'].dependency_names}")
    #print(f"Is it a timesum: {variables['testD'].is_timesum}")
    #variables["testD"].executeEquation()
    #print()
    
    #print(variables['PR'].sympy_equation)
    #print(variables['testD'].sympy_equation)
    
    
    #variables['testD'].executeEquation()
    #print(variables['testD'].values)
    #print(variables['T'].values)
    
    #variables['PR'].executeEquation()
    #print(variables['PR'].values)
    #variables['PR_temp_corr'].executeEquation()
    #print(variables['PR_temp_corr'].values)
    

    
    
    
    
    