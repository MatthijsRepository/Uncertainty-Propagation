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
    

    derived_variables = equation_engine.splitBasicDerived(variables)[1]    
    for var_name in derived_variables:
        calculation_engine.evaluateVariable(variables[var_name], update_var=True, silent=True)  

    


    uncertainty_engine = UncertaintyEngine(variables, equation_engine = equation_engine, calculation_engine = calculation_engine)
    uncertainty_engine.prepareAllDirectUncertainties(variables.values())
    
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
    
    
    
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('Pout' / 'P0')"], recurse=True)
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('G' / 'G_STC')"], recurse=True)
    
    
    variables["TS_('Pout' / 'P0')"].giveReport(short_report=True)
    variables["TS_('G' / 'G_STC')"].giveReport(short_report=True)
    
    
    
    
    print(variables['PR'].uncertainty.getSource("directional response"))
    
    
    