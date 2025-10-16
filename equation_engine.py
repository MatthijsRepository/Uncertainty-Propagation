from my_dataclasses import Uncertainty, Variable
import re


""" This engine handles the initialization of the equation tree, 
    verifies its internal consistency and populates missing basic values """
class EquationEngine:
    def __init__(self, variables):
        self.variables = variables          #dict: dictionary of variable names and Variable objects
        self.basic_variables, self.derived_variables = self.splitBasicDerived()
        for name in self.derived_variables:
            self.checkVariableEquationConsistency(self.variables[name])
        
    def splitBasicDerived(self, variables=None):
        """ Splits a tree into basic and derived variables """
        #If no set provided, act on own registry
        if variables is None:
            variables = self.variables
        
        basic_variables = []
        derived_variables = []
        #Separate variables into basic and derived variables
        for name, var in variables.items():
            if var.is_basic:
                basic_variables.append(name)
            else:
                derived_variables.append(name)
        return basic_variables, derived_variables

    def equationReader(self, variable):
        """ This function interprets written equations to python """
        equation_variables = re.findall(r"'(.*?)'", variable.equation)                   
        return equation_variables
    
    def checkVariableEquationConsistency(self, variable):
        """ This function verifies whether specified variable dependency matches a given equation """
        equation_variables = self.equationReader(variable)
        if not set(equation_variables).issubset(variable.equation_variables):
            raise ValueError(f"Listed equation variables of variable {variable.name} do not match equation. Please check inputs. \nEquation: {variable.equation} \nDetected variables: {equation_variables} \nListed variables: {variable.equation_variables}")


    def checkEquationTreeConsistency(self, variables=None, variables_to_check=None, silent=True):
        """ This function verifies the equation tree is consistent """
        #If no variables provided: act on own registry
        if variables is None:
            variables = self.variables
        #If no variables to check are provided: split provided variables
        if variables_to_check is None:
            if variables is self.variables:
                variables_to_check = self.derived_variables
            else:
                variables_to_check = self.splitBasicDerived(variables)
  
        #If variables_to_check is a single variable, we convert it to a list
        if isinstance(variables_to_check, str):
            variables_to_check = [variables_to_check]

        
        #Recursively navigates down the tree and checks if the equation tree is defined in a consistent way
        for name in variables_to_check:
            #Check if variable in variables list
            if name not in variables:
                raise ValueError(f"Equation tree consistency check failed: variable {name} not included in the variable set. \nVariable set: {variables.keys()}.")
            #If variable already checked, pass this variable
            if variables[name].is_basic:
                continue
            if variables[name].is_root_consistent: 
                continue
            #Else: loop through all variables in dependencies
            for var in variables[name].equation_variables:
                #If listed variable not included in variables list then the equation tree is not consistent
                if var not in variables:
                    raise ValueError(f"Equation tree of variable {name} not consistent: variable {var} not recognized.")
                #If variable is basic we can continue
                if variables[var].is_basic:
                    continue
                #If tree section is verified we can continue
                elif variables[var].is_root_consistent:
                    continue
                #If variable is not basic and not yet verified: check tree consistency
                #If tree is consistent we can continue
                elif self.checkEquationTreeConsistency(variables, variables[var].equation_variables, silent=silent):
                    continue
                else:
                    raise ValueError(f"Equation of {name} is not root-consistent: check for {var} failed")
            #If we succesfully looped through all dependencies of this variable, we can flag it as safe, and potentially log this to output
            if not silent:
                print(f"{name} is root-consistent")
            variables[name].is_root_consistent = True
        #If no errors were encountered we can safely return True
        return True
                

    
                
            
            
                            







