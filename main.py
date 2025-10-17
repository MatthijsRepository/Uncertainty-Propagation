from input_handler import InputHandler
from equation_engine import EquationEngine


if __name__=="__main__":
    inputfile = "C:\\Users\\mate\\Desktop\\python\\Uncertainty-Propagation\\myinput.txt"
    
    #Handling input
    print("Started parsing input")
    input_handler = InputHandler()
    input_handler.parse(inputfile)
    variables = input_handler.variables
    del input_handler
    print("Ended parsing input \n")
    
    #for var in variables.values():
    #    print(var)
    
    #Verifying equation tree consistency, building equation tree in SimPy
    print("Verifying root-consistency of equation tree")
    equation_engine = EquationEngine(variables, update_dependencies=True)
    #equation_engine.checkEquationTreeConsistency(variables="T")
    #equation_engine.checkEquationTreeConsistency(variables_to_check='PR')
    equation_engine.checkEquationTreeConsistency(silent=False)
    print("Equation tree is root-consistent \n")
    
    equation_engine.populateEquationTreeDependencies()
    
    print(variables.keys())
    
    print("Building sympy equation tree")
    equation_engine.buildSymPyEquationTree()
    
    print(variables["testC"].dependency_names)
    #equation_engine.buildEquation(variables["testC"])
    
    variables['testC'].executeEquation()
    print(variables['testC'].values)
    
    variables['C_25'].executeEquation()
    print(variables['T'].values[0])
    print(variables['C_25'].values[0])
    
    print(variables['PR'].sympy_equation)
    #print(variables['testD'].sympy_equation)
    
    
    #variables['testD'].executeEquation()
    #print(variables['testD'].values)
    #print(variables['T'].values)
    
    variables['PR'].executeEquation()
    print(variables['PR'].values)
    #variables['PR_temp_corr'].executeEquation()
    #print(variables['PR_temp_corr'].values)
    
    
    
    
    