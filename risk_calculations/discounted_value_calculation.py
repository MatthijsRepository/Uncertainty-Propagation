import numpy as np

#### Small Explainer ####

#We define the start of the operational lifetime as the date at which the final acceptance certificate is issued: 25 - 2 years
#Moreover, we define the date at which the final acceptance certificate is issued as the start of the 'fiscal year'
#We assume PPA payments are given at the end of each fiscal year
#At the end of each fiscal year, the degradation of the power plant is increased

#So, the first payment will already have to be discounted by the discount rate to the power 1
#But, the first year of energy yield is not yet affected by the degradation

#NPV discount is exponential in time, just as in regular NPV calculations
#Degradation is linear over time and is an absolute percentage subtraction from the PR at the end of each year

#########################



rating = 100 * 10**6         #Rating in W
c_cost = 1.61                  #Construction cost in $/W
acceptance_frac = 0.10      #Fraction of contract value related to FAC

contract_value = rating * c_cost #Total contract value in $


tolerance_level = 0.97      #Tolerance level in %/100
PR_G = 0.90                 #Guaranteed PR
PR_T = 0.90                 #True PR
PR_M = 0.89                #Measured PR
U_PR = 0.022                #PR uncertainty at k=1


deg_rate = 0.016            #Plant degradation rate in %/100 / year
dis_rate = 0.02             #Discount rate in %/100 / year
op_lifetime = 23            #Total operational lifetime in years


#STC_hours = 1642.5            #Number of kWh/m^2/year solar irradiance available

STC_hours_day = 6
STC_hours = 365*STC_hours_day


revenue_per_MWh = 29        #Revenue in $/MWh delivered to the grid


def calculateYearYield(PR, STC_hours, rating, deg_rate=None, years_old=None):
    if deg_rate is not None:
        #PR *= 1 - deg_rate * years_old
        rating *= 1 - deg_rate * years_old
    return PR * rating * STC_hours

def calculateYearMissedYield(PR_G, PR_M, STC_hours, rating, deg_rate=None, years_old=None):
    if deg_rate is not None:
        return (PR_G - PR_M) * STC_hours * rating * (1 - deg_rate * years_old)
    else:
        return (PR_G - PR_M) * STC_hours * rating


def calculateLifeTimeYield(PR, STC_hours, rating, deg_rate, op_lifetime, total_sum=False):
    years_old = np.arange(op_lifetime)
    #yearly_PR = PR - deg_rate*years_old
    yearly_rating = rating * (1 - deg_rate*years_old)
    if total_sum:
        #return np.sum(yearly_PR * STC_hours * yearly_rating)
        return np.sum(PR * STC_hours * yearly_rating)
    else:
        #return yearly_PR *  STC_hours * yearly_rating
        return PR * STC_hours * yearly_rating
    

def calculate_NPV_LifeTimeYield(lifetime_yield, revenue_per_MWh, dis_rate, total_sum=False):
    payments = lifetime_yield * revenue_per_MWh / 10**6
    years = np.arange(1,len(payments)+1)
    discount_rates = (1+dis_rate)**(-1*years)
    if total_sum:
        return np.sum(payments * discount_rates)
    else:
        return payments * discount_rates





meas_lifetime_yield = calculateLifeTimeYield(PR_M, STC_hours, rating, deg_rate, op_lifetime, total_sum=True)
guar_lifetime_yield = calculateLifeTimeYield(PR_G, STC_hours, rating, deg_rate, op_lifetime, total_sum=True)

lifetime_missed_yield = guar_lifetime_yield - meas_lifetime_yield

lifetime_guar_income = guar_lifetime_yield / 10**6 * revenue_per_MWh
lifetime_meas_income = meas_lifetime_yield / 10**6 * revenue_per_MWh
lifetime_missed_income = lifetime_missed_yield / 10**6 * revenue_per_MWh

print()
print(f"Total Plant Capex (EPC Contract Value): {contract_value / 10**6} [$million]")
print()
print(f"Theoretical lifetime revenue (guaranteed): {lifetime_guar_income/10**6} [$million]")
print(f"Theoretical lifetime revenue (meas): {lifetime_meas_income/10**6} [$million]")
print(f"Lifetime revenue (meas) - EPC_capex: {(lifetime_meas_income - contract_value)/10**6} [$million]")
print()
print(f"Theoretical lifetime missed revenue: {lifetime_missed_income/10**6} [$million]")
print(f"Theoretical lifetime missed revenue as percentage of contract value: {lifetime_missed_income / contract_value * 100} %")








