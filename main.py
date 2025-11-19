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
    
    #"""
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
    
    
    #variables['G'].giveReport()
    variables['PR'].giveReport(decimals=5)
    variables['PR_temp_corr'].giveReport(decimals=5)
    
    
    
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('Pout' / 'P0')"], recurse=True)
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('G' / 'G_STC')"], recurse=True)
    
    
    #variables["TS_('Pout' / 'P0')"].giveReport(short_report=True)
    variables["TS_('G' / 'G_STC')"].giveReport(short_report=True)
    
    
    variables['G'].plotValuesAndUncertainty(k=2)
    
    
    
    
    
    
    
    #print(variables['PR'].uncertainty.getSource("directional response"))
    import matplotlib.pyplot as plt
    plt.plot(variables['PR'].uncertainty.getSource("directional response").values)
    plt.show()
    
    
    
    
    
    """
    print()
    print([source.name for source in variables["TS_('G' / 'G_STC')"].uncertainty.root_sources])
    
    import numpy as np
    
    
    
            #          sum of V                dt    S^2        G_STC
    factor =  np.sum(variables['V'].values) * 60 / ( 15**2 * 1000) 
    print(variables["TS_('G' / 'G_STC')"].uncertainty.aggregated_weighted_uncertainties[0] / factor /15*3) #Should be 0.015 - 1.5% bound
    
    
    
    print()
    uncertainty_engine.calculateTotalUncertainty(variables['TSGC'])
    
    print([source.name for source in variables["TSGC"].uncertainty.root_sources])
    
    
    factor = -1*variables['gamma'].values * np.sum(variables['G'].values) * 60
    print(variables['TSGC'].uncertainty.aggregated_weighted_uncertainties[0] / factor)
    
    print()
    uncertainty_engine.calculateTotalUncertainty(variables['timesumtest'])
    
    print([source.name for source in variables["timesumtest"].uncertainty.root_sources])
    
    factor = np.sum(variables['G'].values) * 60 #* -1 * variables['gamma'].values
    
    print(variables['timesumtest'].uncertainty.aggregated_weighted_uncertainties[0] / factor)
    
    
    
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('C_25')"])
    print(variables["TS_('C_25')"].values)
    print(variables["TS_('C_25')"].uncertainty.total_uncertainty)
    """
    
    
    """
    
    uncertainty_engine.calculateTotalUncertainty(variables['timesumtest'], recurse=True)
    uncertainty_engine.calculateTotalUncertainty(variables['C_25'], recurse=True)
    uncertainty_engine.calculateTotalUncertainty(variables["TS_('C_25')"], recurse=True)
    print(variables["TS_('C_25')"].uncertainty.total_uncertainty)
    import numpy as np
    
    print([source.name for source in variables['timesumtest'].uncertainty.root_sources])
    
    print(np.average(variables['C_25'].values))
    
    print(np.average(variables['C_25'].uncertainty.total_uncertainty))



    print(variables['timesumtest'].uncertainty.aggregated_weighted_uncertainties[0] / variables["TSG"].values /60)
    
    
    variables["TS_('C_25')"].giveReport()
    
    
    print()
    uncertainty_engine.calculateTotalUncertainty(variables['TStest_agg'], recurse=True)
    print([source.name for source in variables['TStest_agg'].uncertainty.root_sources])
    
    print(variables['TStest_agg'].uncertainty.aggregated_weighted_uncertainties[0])
    print(variables['TStest_agg'].uncertainty.aggregated_weighted_uncertainties[0] / variables['PR'].values / abs(variables['gamma'].values))
    
    #Check if calibration is as expected
    print(variables['TStest_agg'].uncertainty.aggregated_weighted_uncertainties[2] / np.average(variables['C_25'].values) / abs(variables['gamma'].values))
    #variables['timesumtest'].giveReport()
    
    
    
    uncertainty_engine.calculateTotalUncertainty(variables['PR'])
    print([source.name for source in variables['PR'].uncertainty.root_sources])
    print(variables['PR'].uncertainty.aggregated_weighted_uncertainties[1])
    #print(variables['PR'].uncertainty.getSource("directional response"))
    """
    
    
    