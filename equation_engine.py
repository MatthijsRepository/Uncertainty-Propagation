from my_dataclasses import Variable
import re
import sympy as sp


""" This engine handles the initialization of the equation tree, 
    verifies its internal consistency and populates missing basic values """
class EquationEngine:
    def __init__(self, variables, update_dependencies=False):
        self.variables = variables          #dict: dictionary of variable names and Variable objects
        self.basic_variables, self.derived_variables = self.splitBasicDerived()
        for name in self.derived_variables:
            self.checkVariableEquationConsistency(self.variables[name], update_dependencies)
        
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
    
    def checkVariableEquationConsistency(self, variable, update_dependencies=False):
        """ This function verifies whether specified variable dependency matches a given equation, optionally updates dependencies of variable to detected variables """
        equation_variables = self.equationReader(variable)
        if not set(equation_variables).issubset(variable.equation_variables):
            raise ValueError(f"Listed equation variables of variable {variable.name} do not match equation. Please check inputs. \nEquation: {variable.equation} \nDetected variables: {equation_variables} \nListed variables: {variable.equation_variables}")
        #update equation variables to detected variables
        if update_dependencies:
            variable.equation_variables = equation_variables
            
    def updateEquationVariables(self, variables=None):
        """ Updates listed equation variables for a variable set to only include those detected in the listed equation """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        else:
            derived_variables = self.splitBasicDerived(variables)[1]
            
        for name in derived_variables:
            var = variables[name]
            equation_variables = self.equationReader(var)
            var.equation_variables = equation_variables
            
    def _checkEquationTreeRecursive(self, variables, variables_to_check, silent):
        """ Internal function that recursively checks equation tree consistency """
        #Recursively navigates down the tree and checks if the equation tree is defined in a consistent way
        for name in variables_to_check:
            var = variables.get(name)
            #Check if variable in variables list
            if var is None:
                raise ValueError(f"Equation tree consistency check failed: variable {name} not included in the variable set. \nVariable set: {variables.keys()}.")
            
            #If variable already checked, pass this variable
            if var.is_basic:
                continue
            if var.is_root_consistent: 
                continue
            #Else: loop through all variables in dependencies
            for dep_name in var.equation_variables:
                dep = variables.get(dep_name)
                #If listed variable not included in variables list then the equation tree is not consistent
                if dep is None:
                    raise ValueError(f"Equation tree of variable {name} not consistent: variable {dep_name} not recognized.")
                #If variable is basic we can continue
                if dep.is_basic:
                    continue
                #If tree section is verified we can continue
                elif dep.is_root_consistent:
                    continue
                #If variable is not basic and not yet verified: check tree consistency
                elif self._checkEquationTreeRecursive(variables, dep.equation_variables, silent):
                    continue
                else:
                    raise ValueError(f"Equation of {name} is not root-consistent: check for {dep_name} failed")
            #If we succesfully looped through all dependencies of this variable, we can flag it as safe, and potentially log this to output
            if not silent:
                print(f"{name} is root-consistent")
            var.is_root_consistent = True
        #If no errors were encountered we can safely return True
        return True


    def checkEquationTreeConsistency(self, variables=None, variables_to_check=None, silent=True):
        """ This function that handles input and executes the equation tree consistency checker """
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
        
        return self._checkEquationTreeRecursive(variables, variables_to_check, silent)
    
    
    def buildEquation(self, var, symbol_map):
        #Get symbol list for this variable
        symbols = [symbol_map[name] for name in var.equation_variables]
        #Remove quotes from equation
        expr_str = var.equation
        for dep in var.equation_variables:
            expr_str = expr_str.replace(f"'{dep}'", dep)
        
        #Create sympy equation - this is later also used for differentiation
        var.sympy_equation = sp.sympify(expr_str, locals=symbol_map)
        
        #Create callable function
        executable_equation = sp.lambdify(var.equation_variables, var.sympy_equation)
        def wrapper(variable_map):
            return
            
        
    
    def buildSymPyEquationTree(self, variables=None):
        if variables is None:
            variables = self.variables
        
        symbol_map = {name: sp.Symbol(name) for name in variables.keys()}
        for name, var in variables.items():
            if not var.is_basic:
                expr_str = var.equation
                for dep in var.equation_variables:
                    expr_str = expr_str.replace(f"'{dep}'", dep)
                
                var.sympy_equation = sp.sympify(expr_str, locals=symbol_map)
                #var.sympy_equation = sp.lambdify(symbol_map)

                print(var.sympy_equation)
                var.executable_equation = sp.lambdify(var.equation_variables, var.sympy_equation)

                

    
                
            
            
                            







