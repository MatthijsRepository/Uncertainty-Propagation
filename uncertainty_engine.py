from my_dataclasses import Variable, UncertaintySource
import numpy as np


class UncertaintyEngine:
    def __init__(self, variables, equation_engine=None, calculation_engine=None):
        self.variables = variables
        self.equation_engine = equation_engine
        self.calculation_engine = calculation_engine
        
     
        
    def _calculateDirectUncertainty(self, var):
        """ Calculates the magnitudes of the direct uncertainties for each source and returns these in an array; or 0 if there is no uncertainty """
        if len(var.uncertainty.direct_uncertainty_sources)==0:
            return 0
        else:
            #Prepare uncertainty array
            if isinstance(var.values, np.ndarray):
                u_array = np.zeros((len(var.uncertainty.direct_uncertainty_sources), len(var.values)))
            else: #values is float or int
                u_array = np.zeros(len(var.uncertainty.direct_uncertainty_sources))
            #Fill the uncertainty array
            for index, u in enumerate(var.uncertainty.direct_uncertainty_sources):
                if u.is_relative is False:
                    u_array[index] = u.sigma #for absolute uncertainty: cast value directly into array
                else:
                    u_array[index] = u.sigma * var.values
            return u_array
                    
    
    def _retrieveDependencyUncertainties(self, var, recurse):
        if var.is_basic:
            return 0            
        #Initialize dictionary of dependencies and associated uncertainties
        total_dep_uncertainties = {}
        #Loop through dependencies and populate uncertainties dictionary
        for dep in var.dependencies.values():
            #Check if dependency has uncertainty sources and whether these are calculated already
            if not dep.uncertainty.is_calculated:
                if recurse:
                    self.calculateUncertainty(dep, recurse)
                else:
                    raise ValueError(f"Uncertainty of dependency {dep.name} of variable {var.name} not yet calculated, please ensure this or enable recursive calculation")
            total_dep_uncertainties[dep] = dep.uncertainty.total_uncertainty
        return total_dep_uncertainties
            
            
    def _getDependencyPartialsValues(self, var, equation_engine=None):
        if equation_engine is None:
            equation_engine = self.equation_engine
        
        #Populate variable partial executables
        equation_engine.buildPartialDerivativeExecutables(var)
        return var.executeAllPartials(store_results=False, force_recalculation=False)
    
    def calculateUncertainty(self, var, recurse=False):
        #Populate direct uncertainty array in variable uncertainty data
        var.uncertainty.direct_uncertainty_contributions = self._calculateDirectUncertainty(var)
        #Calculate uncertainty due to dependencies
        if not var.is_basic:
            total_dep_uncertainties = self._retrieveDependencyUncertainties(var, recurse)
            dep_partial_values = self._getDependencyPartialsValues(var)
            
                
            
        """
        for i, source in enumerate(var.uncertainty.direct_uncertainty_sources):
            print(source.name)
            print(f"Is onesided: {source.is_symmetric}")
            print(var.uncertainty.direct_uncertainty_contributions[i,0])
        """
        
        #Flag variable uncertainty to be calculated
        var.uncertainty.is_calculated = True
            



