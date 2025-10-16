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
    equation_engine.checkEquationTreeConsistency()
    print("Equation tree is root-consistent \n")
    
    equation_engine.populateEquationTreeDependencies()
    
    print("Building sympy equation tree")
    equation_engine.buildSymPyEquationTree()
    
    print(variables["testC"].dependency_names)
    #equation_engine.buildEquation(variables["testC"])
    #print(variables['testC'].executable_equation())
    
    
    
    
    