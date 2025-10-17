from my_dataclasses import Variable


def TIMESUM_TEMP(var, options=None):
    print("TIMESUM_TEMP test")
    print(var)
    if isinstance(var, (float,int)):
        return var
    elif isinstance(var.values, (float, int)):
        raise ValueError(f"Timesums of constants are not supported yet (variable {var.name})")
    else:
        return var.TimeSum()
    
    
    
class CalculationEngine:
    def __init__(self):
        return
    
    def add(self, var1, var2):
        return var1 + var2
    
    def timeSum(var):
        if isinstance(var, (float,int)):
            return var
        elif isinstance(var.values, (float, int)):
            raise ValueError(f"Timesums of constants are not supported yet (variable {var.name})")
        else:
            return var.TimeSum()


