import sympy as sp

class timesum(sp.Function):
    """ Symbolic SymPy representation of a timesum routine, required for equation interpretation and differentiation. For evaluation, see the calculation engine """
    
    @classmethod
    def eval(cls, arg):
        #Never automatically simplify, always keep timesum as a node in the tree.
        return None
    
    def _eval_derivative(self, sym):
        # d/dsym of timesum(f) = timesum(df/dsym)
        return timesum(sp.diff(self.args[0], sym), *self.args[1:])
    
    
    
    
    
    
    
    
    
    
    