import numpy as np
from scipy.integrate import quad, dblquad
import matplotlib.pyplot as plt


#####################################################
############### Basic probability functions 
#####################################################

#Calculation of PR_R
def upperLimit(case, PR_G, delta=None, u_PR=None, k=None):
    #Calculates PR_R for different acceptance testing types
    
    #Basic methods
    if case==1:
        return PR_G - 100*delta
    if case==2:
        return PR_G + k*u_PR
    
    #Optimal error overall
    if case==3:
        return PR_G + -0.318703*max(u_PR-1.349325, 0)
    if case==4: 
        return PR_G + -0.346437*max(u_PR-0.999948, 0)
    

#Probability density function of PR_M
def pdf(PR_M, PR_T, PR_R, u_PR):
    return (1/np.sqrt(2 * np.pi * u_PR**2)) * np.exp(-1* (PR_M - PR_T)**2/(2*u_PR**2))

#Probability density of PR_M multiplied by payable damages at the measured value
def expectedDamages(PR_M, PR_T, PR_R, u_PR):
    return (PR_R - PR_M) * (1/np.sqrt(2 * np.pi * u_PR**2)) * np.exp(-1* (PR_M - PR_T)**2/(2*u_PR**2))




#Calculate the heatmap values of pdf or expected loss using this function
def calculateHeatmap(function, PR_T_series, PR_G, u_PR_series, case, delta=None, k=None):
    results = np.zeros((len(u_PR_series), len(PR_T_series)))
    
    for i in range(len(u_PR_series)):
        PR_R = upperLimit(case=case, PR_G=PR_G, delta=delta, u_PR=u_PR_series[i], k=k)
        
        for j in range(len(PR_T_series)):
            a, b = quad(function, 0, PR_R, args=(PR_T_series[j], PR_R, u_PR_series[i])) 
            results[i,j] = a
    return results


#Calculate the heatmap of error in expected payable damages
def calculateErrormap(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None):
    results = np.zeros((len(u_PR_series), len(PR_T_series)))
    
    PR_R = None
    if case==1:
        PR_R = upperLimit(case=1, PR_G=PR_G, delta=delta)
    
    for i in range(len(u_PR_series)):
        if case==2:
            PR_R = upperLimit(case=2, PR_G=PR_G, u_PR=u_PR_series[i], k=k)
        
        for j in range(len(PR_T_series)):
            a, b = quad(expectedDamages, 0, PR_R, args=(PR_T_series[j], PR_R, u_PR_series[i])) 
            
            fair_loss = PR_T_series[j] - PR_G
            fair_loss = min(fair_loss, 0)
            
            results[i,j] = a +fair_loss
    return results


#Calculate the cross-section of a heatmap at a given uncertainty level
def calculate_line_uPR(function, PR_T_series, PR_G, u_PR, case, delta=None, k=None):
    results = np.zeros(len(PR_T_series))
    PR_R = upperLimit(case=case, PR_G=PR_G, u_PR=u_PR, delta=delta, k=k)
    
    for j in range(len(PR_T_series)):
        a, b = quad(function, 0, PR_R, args=(PR_T_series[j], PR_R, u_PR)) 
        results[j] = a
    return results


#####################################################
############### Value-at-Risk
#####################################################


#Value at risk calculations
def calculateValueAtRisk(PR_T, PR_R, u_PR, risk_level=5):
    #Relevant k-value for given percentage
    k_dict = {1 : 2.33, 5: 1.64, 10: 1.28}
    
    k = k_dict.get(risk_level)
    if k is None:
        raise ValueError("given percentage value-at-risk not implemented!")
        
    #Get PR_M at the 5% level
    level = PR_T - k*u_PR
    #Calculate value using cost function
    value_at_risk = -1*(level - PR_R)
    return max(value_at_risk, 0)

#Calculate a heat map for the value-at-risk
def calculateValueAtRiskHeatmap(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None, risk_level=5):
    #Initialize results array
    results = np.zeros((len(u_PR_series), len(PR_T_series)))
    
    for i in range(len(u_PR_series)):
        #Calcualte PR_R
        PR_R = upperLimit(case=case, PR_G=PR_G, delta=delta, u_PR=u_PR_series[i], k=k)
        for j in range(len(PR_T_series)):
            #Populate results
            results[i,j] =  calculateValueAtRisk(PR_T_series[j], PR_R, u_PR_series[i], risk_level=risk_level) 
    return results


def plotValueAtRiskImpact(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None, risk_level=5, colors=None, labels=None, limits=None):
    #Building the results
    lines = []
    for u_PR in u_PR_series:
        results = np.zeros(len(PR_T_series)) 
        PR_R = upperLimit(case=case, PR_G=PR_G, delta=delta, u_PR=u_PR, k=k)
        for j in range(len(PR_T_series)):
            results[j] = calculateValueAtRisk(PR_T_series[j], PR_R, u_PR, risk_level)
        lines.append(results)
    
    #Plotting the results
    fig, ax = plt.subplots(dpi=200)
    ax.set_aspect('equal', adjustable='box')
    for i in range(len(lines)):
        
        if labels is not None:
            label = r"$u_{PR}$="+f"{u_PR_series[i]*2}\n" + labels[i]
        else:
            label=r"$u_{PR}$="+f"{u_PR_series[i]*2}"
            
        if colors is None:
            ax.plot(PR_T_series-PR_G , lines[i] , label=label)
        else:
            ax.plot(PR_T_series-PR_G , lines[i] , label=label, color=colors[i])
        if i%2==0:
            diff = lines[i+1]-lines[i]
            max_index = np.argmax(diff)
            print(f"Value at risk absolute difference {i//2+1}: {diff[max_index]}")
    
    
    ax.set_xlim(PR_T_series[0]-PR_G, PR_T_series[-1]-PR_G)
    ax.set_ylim(-0.3001772878604755, 10.303725807641181)
    
    ax.set_xlabel(r"$PR_T - PR_R$")
    ax.set_ylabel(str(risk_level) + r"% value-at-risk / $\xi$")
    ax.grid()
    
    ax.legend(bbox_to_anchor=(1,0.5), loc="center left")
    #ax.legend(bbox_to_anchor=(0.5,1.02), loc="lower center", ncols=2)
    
    plt.show()
    return

def plotValueAtRiskHeatmap(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None, risk_level=5, contour_levels=None, color_range=None):
    results = calculateValueAtRiskHeatmap(PR_T_series, PR_G, u_PR_series, case, delta=delta, k=k, risk_level=risk_level)
    
    #Plotting the results
    fig, ax = plt.subplots(dpi=200)
    u_PR_series *= 2   
    domain = (PR_T[0]-PR_G, PR_T[-1]-PR_G, u_PR_series[0], u_PR_series[-1])
    aspect_ratio = (domain[1]-domain[0]) / (domain[3]-domain[2]) * 0.7
    
    if color_range is None:
        im = ax.imshow(results, cmap='magma', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain)
    else:
        im = ax.imshow(results, cmap='magma', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain,
                       vmin=color_range[0], vmax=color_range[1])
    
    #"""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "top",           # position
        size="6%",       # height relative to axes
        pad=0.15         # gap between axes and colorbar
    )
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.invert_xaxis()
    cbar.set_label(fr"$\frac{{{risk_level}\text{{\% Value-at-risk}}}}{{\xi}}$", rotation=0, labelpad=10, fontsize=14)
    #cbar.set_label(r"$\frac{\text{5% value-at-risk}}{\xi}$", rotation=0, labelpad=10, fontsize=14)
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")
    #"""
    #cbar = fig.colorbar(im, ax=ax)
    #cbar.set_label(r"$\frac{\text{expected loss}}{\xi}$", rotation=0, labelpad=40, fontsize=14)
    
    ax.set_yticks([2,4,6,8])
    
    ax.set_xlabel(r"$PR_T - PR_G$ [%]", fontsize=12)
    ax.set_ylabel(r"2$u_{PR}$"+"\n [%]", rotation=0, labelpad=15, fontsize=12)
    
    
    if contour_levels is not None:
        X, U = np.meshgrid(PR_T-PR_G, u_PR)
        contour_levels = cbar.get_ticks()
        cs = ax.contour(X, U, results,
                        levels=contour_levels,
                        colors="white",
                        linewidths=0.9)
    
    plt.show()
    u_PR_series /= 2
    
    return


#####################################################
############### Failure probabilty plots
#####################################################


def plotFailProbability(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None):
    results = calculateHeatmap(pdf, PR_T_series, PR_G, u_PR_series, case, delta, k)
    
    fig, ax = plt.subplots(dpi=200)
    #fig = plt.figure(figsize=(15,5), dpi=200)
    #ax = fig.add_axes([0.12, 0.15, 0.76, 0.65])
    
    u_PR_series *= 2
    domain = (PR_T[0]-PR_G, PR_T[-1]-PR_G, u_PR_series[0], u_PR_series[-1])
    aspect_ratio = (domain[1]-domain[0]) / (domain[3]-domain[2]) * 0.7
    

    
    im = ax.imshow(results, cmap='magma', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain,
                   vmin=0, vmax=1)
    
    #"""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "top",           # position
        size="6%",       # height relative to axes
        pad=0.15         # gap between axes and colorbar
    )
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal", fraction=0.046)
    cbar.ax.invert_xaxis()
    cbar.set_label(r"$\text{P}(\text{fail})$", rotation=0, labelpad=10, fontsize=14)
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")
    #"""
    
    #cbar = fig.colorbar(im, ax=ax)
    #cbar.set_label(r"$\text{P}(\text{fail})$", rotation=0, labelpad=20, fontsize=12)

    ax.set_yticks([2,4,6,8])
    
    ax.set_xlabel(r"$PR_T - PR_G$ [%]", fontsize=12)
    ax.set_ylabel(r"2$u_{PR}$"+"\n [%]", rotation=0, labelpad=15, fontsize=12)
    
    X, U = np.meshgrid(PR_T-PR_G, u_PR_series)
    cs = ax.contour(X, U, results,
                    levels=[0.2, 0.4, 0.6, 0.8],
                    colors="white",
                    linewidths=0.9)
        
    plt.show()
    u_PR_series /= 2 
    return


#####################################################
############### Expected damages plots
#####################################################


def plotExpectedDamages(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None, contour_levels=None, color_range=None):
    results = calculateHeatmap(expectedDamages, PR_T_series, PR_G, u_PR_series, case, delta=delta, k=k)
    
    fig, ax = plt.subplots(dpi=200)
    
    
    #define color range:
    if case==2 and k>0:
        loss_color_range = [0,12.5]
    else:
        loss_color_range = [0, 5.5]
    
    
    u_PR_series *= 2   
    domain = (PR_T_series[0]-PR_G, PR_T_series[-1]-PR_G, u_PR_series[0], u_PR_series[-1])
    aspect_ratio = (domain[1]-domain[0]) / (domain[3]-domain[2]) * 0.7
    
    if color_range is None:
        im = ax.imshow(results, cmap='magma', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain)
    else:
        im = ax.imshow(results, cmap='magma', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain,
                       vmin=color_range[0], vmax=color_range[1])
    
    #"""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "top",           # position
        size="6%",       # height relative to axes
        pad=0.15         # gap between axes and colorbar
    )
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.invert_xaxis()
    cbar.set_label(r"$\frac{\text{expected loss}}{\xi}$", rotation=0, labelpad=10, fontsize=14)
    #cbar.set_label(r"$\frac{\text{5% value-at-risk}}{\xi}$", rotation=0, labelpad=10, fontsize=14)
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")
    #"""
    #cbar = fig.colorbar(im, ax=ax)
    #cbar.set_label(r"$\frac{\text{expected loss}}{\xi}$", rotation=0, labelpad=40, fontsize=14)
    
    ax.set_yticks([2,4,6,8])
    
    ax.set_xlabel(r"$PR_T - PR_G$ [%]", fontsize=12)
    ax.set_ylabel(r"2$u_{PR}$"+"\n [%]", rotation=0, labelpad=15, fontsize=12)
    
    
    X, U = np.meshgrid(PR_T_series-PR_G, u_PR_series)
    contour_levels = cbar.get_ticks()
    cs = ax.contour(X, U, results,
                    levels=contour_levels,
                    colors="white",
                    linewidths=0.9)
    
    plt.show()
    u_PR_series /= 2


#####################################################
############### Expected error plots
#####################################################

def plotExpectedError(PR_T_series, PR_G, u_PR_series, case, delta=None, k=None):
    results = calculateErrormap(PR_T_series, PR_G, u_PR_series, case, delta=delta, k=k)
    
    fig, ax = plt.subplots(dpi=200)
    
    
    #define color range:
    color_range = np.array([-1, 1]) * 4
    
    
    u_PR_series *= 2   
    domain = (PR_T_series[0]-PR_G, PR_T_series[-1]-PR_G, u_PR_series[0], u_PR_series[-1])
    aspect_ratio = (domain[1]-domain[0]) / (domain[3]-domain[2]) * 0.7
    
    if color_range is None:
        im = ax.imshow(results, cmap='seismic', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain)
    else:
        im = ax.imshow(results, cmap='seismic', interpolation='nearest', origin='lower', aspect=aspect_ratio, extent=domain,
                       vmin=color_range[0], vmax=color_range[1])
    
    #"""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "top",           # position
        size="6%",       # height relative to axes
        pad=0.15         # gap between axes and colorbar
    )
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.invert_xaxis()
    cbar.set_label(r"$\frac{\text{expected error}}{\xi}$", rotation=0, labelpad=10, fontsize=14)
    #cbar.set_label(r"$\frac{\text{5% value-at-risk}}{\xi}$", rotation=0, labelpad=10, fontsize=14)
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")
    #"""
    #cbar = fig.colorbar(im, ax=ax)
    #cbar.set_label(r"$\frac{\text{expected loss}}{\xi}$", rotation=0, labelpad=40, fontsize=14)
    
    ax.set_yticks([2,4,6,8])
    
    ax.set_xlabel(r"$PR_T - PR_G$ [%]", fontsize=12)
    ax.set_ylabel(r"2$u_{PR}$"+"\n [%]", rotation=0, labelpad=15, fontsize=12)
    
    
    X, U = np.meshgrid(PR_T_series-PR_G, u_PR_series)
    contour_levels = cbar.get_ticks()
    cs = ax.contour(X, U, results,
                    levels=contour_levels,
                    colors="white",
                    linewidths=0.9)
    
    plt.show()
    u_PR_series /= 2



#####################################################
############### Expected damages impact plots
#####################################################

def plotHeatmapHorizontalCrossSection(function, PR_T_series, PR_G, u_PR_series, case, delta=None, k=None, colors=None, labels=None):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
    
    lines = []
    for u_PR in u_PR_series:
        data = calculate_line_uPR(function, PR_T_series, PR_G, u_PR, case, delta=delta, k=k)
        
        #data += np.minimum(PR_T_series-PR_G, 0)
        
        lines.append(data)
    
    fig, ax = plt.subplots(dpi=200)
    ax.set_aspect('equal', adjustable='box')
    
    for i in range(len(lines)):
        
        if labels is not None:
            label = r"$u_{PR}$="+f"{u_PR_series[i]*2}\n" + labels[i]
        else:
            label=r"$u_{PR}$="+f"{u_PR_series[i]*2}"

        if colors is None:
            ax.plot(PR_T_series-PR_G , lines[i] , label=label)
        else:
            ax.plot(PR_T_series-PR_G , lines[i] , label=label, color=colors[i])
            #if i%2==0:
                #ax.plot(PR_T_series-PR_G, (lines[i+1]-lines[i]), label=label, color=colors[i])
                #ax.plot(PR_T_series-PR_G, (lines[i+1]-lines[i])/lines[i+1], label=label, color=colors[i])
    
    
    ax.set_xlim(PR_T_series[0]-PR_G, PR_T_series[-1]-PR_G)
    ax.set_ylim(-0.3001772878604755, 6.303725807641181)
    
    ax.set_xlabel(r"$PR_T - PR_G$")
    ax.set_ylabel(r"expected loss / $\xi$")
    ax.grid()
    
    #ax.legend(bbox_to_anchor=(1,0.5), loc="center left")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    
    axins = inset_axes(ax, width="55%", height="55%", loc="upper right", borderpad=0)
    for i in range(len(lines)):
        if colors is None:
            axins.plot(PR_T_series-PR_G , lines[i] , label=label)
        else:
            axins.plot(PR_T_series-PR_G , lines[i] , label=label, color=colors[i])
            if i%2==0:
                #axins.plot(PR_T_series-PR_G, (lines[i+1]-lines[i]), label=label, color=colors[i])
                diff = lines[i+1]-lines[i]
                max_index = np.argmax(diff)
                print()
                print(f"Maximum difference in expected loss {i//2+1}: {diff[max_index]}")
                print(f"Relative difference in expected loss {i//2+1}: {diff[max_index] / lines[i+1][max_index]}")
                #axins.plot(PR_T_series-PR_G, (lines[i+1]-lines[i])/lines[i+1], label=label, color=colors[i])
    
    #axins.set_xlim(-6, -2)
    axins.set_xlim(-2,2)
    
    axins.set_ylim(0, 2)
    
    axins.set_xticks([])
    axins.set_yticks([])
    
    axins.axvline(-4, color="0.7", linestyle="-", linewidth=0.7)
    axins.axvline(0, color="0.7", linestyle="-", linewidth=0.7)
    axins.axhline(1, color="0.7", linestyle="-", linewidth=0.7)
    
    
    #axins.axhline(0.4, linewidth=0.5)
    
    patch, connector1, connector2 = mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="dimgray", lw=1)
    patch.set_edgecolor("black")
    patch.set_linewidth(1.2)
    
    plt.show()
    return



#####################################################
############### Calculations
#####################################################


PR_G = 90 #-10
PR_R = None

PR_T = np.linspace(84, 96, num=121, endpoint=True) #-10
u_PR = np.linspace(0, 4, num=41, endpoint=True)[2:]
#PR_T = np.linspace(84, 96, num=361, endpoint=True) #-10
#u_PR = np.linspace(0, 4, num=121, endpoint=True)[6:]

#PR_T = 90
#u_PR = 2.5

#case 1 and 2 correspond to acceptance types 1 and 2; case 3 and 4 are two optimized metrics
case = 1
delta = 0
k = 1   #1.96 for P97.5,   1.28 for P90

risk_level = 5



expected_loss_contours = [0.1, 1, 2, 3, 4, 5]
pdf_contours = [0.2, 0.4, 0.6, 0.8]#, 0.975]

loss_color_range_1 = [0, 5.5]
loss_color_range_2 = [0,12.5]
if case == 1:
    loss_color_range = loss_color_range_1
elif case==2:
    loss_color_range = loss_color_range_2
else:
    loss_color_range = None

loss_color_range = loss_color_range_1
#loss_color_range = None

error_color_range = np.array([-1,1]) * 1.5




#plot_u_PR_series = np.array([1,2,3,4,5,6,7,8]) / 2

#basic directional response
plot_u_PR_series = np.array([3.74, 4.64, 3.11, 3.67]) / 2
labels = ["DS1 - halved DiR", "DS1 - standard", "DS2 - halved DiR", "DS2 - standard"]

#sr300
plot_u_PR_series = np.array([3.89, 4.64, 3.03, 3.67]) / 2
labels = ["DS1 - SR300", "DS1 - standard", "DS2 - SR300", "DS2 - standard"]

#sr300 vs sr300 directional response
plot_u_PR_series = np.array([2.84, 3.891, 2.40, 3.03]) / 2
labels = ["DS1 - SR300, halved DiR", "DS1 - SR300", "DS2 - SR300 halved DiR", "DS2 - SR300"]

#sr300 directional response
#plot_u_PR_series = np.array([2.84, 4.64, 2.40, 3.67]) / 2
#labels = ["DS1 - SR300, halved DiR", "DS1 - standard", "DS2 - SR300 halved DiR", "DS2 - standard"]

#class A vs class B
#plot_u_PR_series= np.array([8.91, 4.64, 6.96, 3.67]) / 2
#labels = ["DS1 - class B", "DS1 - standard", "DS2 - class B", "DS2 - standard"]


year = True
if year:
    #basic directional response
    plot_u_PR_series = np.array([3.15, 3.91, 2.97, 3.35]) / 2
    labels = ["DS1 - halved DiR", "DS1 - standard", "DS2 - halved DiR", "DS2 - standard"]

    #sr300
    #plot_u_PR_series = np.array([3.39, 3.91, 2.95, 3.35]) / 2
    #abels = ["DS1 - SR300", "DS1 - standard", "DS2 - SR300", "DS2 - standard"]

    #sr300 vs sr300 directional response
    plot_u_PR_series = np.array([2.48, 3.39, 2.31, 2.95]) / 2
    labels = ["DS1 - SR300, halved DiR", "DS1 - SR300", "DS2 - SR300 halved DiR", "DS2 - SR300"]

    #sr300 directional response
    #plot_u_PR_series = np.array([2.48, 3.91, 2.31, 3.35]) / 2
    #abels = ["DS1 - SR300, halved DiR", "DS1 - standard", "DS2 - SR300 halved DiR", "DS2 - standard"]

    #class A vs class B
    #plot_u_PR_series= np.array([]) / 2
    #labels = ["DS1 - class B", "DS1 - standard", "DS2 - class B", "DS2 - standard"]



colors = ["darkgray", "dimgray", "darkkhaki",  "darkgoldenrod"]
colors = ["tomato", "maroon", "lightskyblue",  "teal"]


#labels=None


#####################################################
############### Value-at-risk plots
#####################################################

print()
plotValueAtRiskImpact(PR_T, PR_G, u_PR_series=plot_u_PR_series, case=case, delta=delta, k=k, colors=colors, labels=labels, risk_level=risk_level)
plotValueAtRiskHeatmap(PR_T, PR_G, u_PR_series=u_PR, case=case, delta=delta, k=k, risk_level=risk_level, contour_levels=[0,2,4,6,8,10,12])

#####################################################
############### Failure probability and expected damages plots
#####################################################

print()
plotFailProbability(PR_T, PR_G, u_PR, case, delta=delta, k=k)
plotExpectedDamages(PR_T, PR_G, u_PR, case, delta=delta, k=k)

plotExpectedError(PR_T, PR_G, u_PR, case, delta=delta, k=k)


plotHeatmapHorizontalCrossSection(expectedDamages, PR_T, PR_G, plot_u_PR_series, case, delta=delta, k=k, colors=colors, labels=labels)



#####################################################
############### Criteria optimization
#####################################################

import sys; sys.exit()

print()
print()
print()
print("Acceptance criterion optimization running, this takes a few minutes...")

def total_error_func(deps, PR_G, PR_T_series, u_PR_series):
    #k, m = deps[0], deps[1]
    delta = deps
    
    results = np.zeros((len(u_PR_series), len(PR_T_series)))
    for i in range(len(u_PR_series)):
        #R_R = PR_G + k*max(u_PR[i]-m, 0)
        PR_R = PR_G - delta
        
        for j in range(len(PR_T_series)):
            a, b = quad(expectedDamages, 0, PR_R, args=(PR_T_series[j], PR_R, u_PR_series[i]))
            
            fair_loss = PR_T_series[j] - PR_G
            fair_loss = min(fair_loss, 0)
            
            results[i,j] = a + fair_loss
            
    total_error = np.sqrt(np.sum(results**2))
    
    mask = (u_PR < 2.25) & (u_PR > 1.25)
    total_error = np.sqrt(np.sum(results[mask,:])**2)
    
    return total_error


from scipy.optimize import minimize



args = (PR_G, PR_T, u_PR)

#x = k, m
x0 = [-0.5, 1]
bounds = [(-1.1, 0), (0,2)]

x0 = [0]
bounds = [(-2, 1)]

result = minimize(total_error_func, x0, args=args, method="Powell", bounds=bounds,
                  options={
    "maxiter": 200,     # hard stop
    "xtol": 1e-4,       # parameter change tolerance
    "ftol": 1e-4,       # function change tolerance
    "disp": True
})


delta_opt = result.x
print(delta_opt)

k_opt, m_opt = result.x
min_loss = result.fun

print(f"k: {k_opt}, m: {m_opt}")
print(f"Loss: {min_loss}")

print(f"Succes: {result.success}")
print(f"Message: {result.message}")
print(f"Iterations: {result.nit}")










