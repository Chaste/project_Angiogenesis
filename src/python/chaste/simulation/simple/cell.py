"""
A collection of simple cell population growth models.
"""

import numpy as np

def exponential_growth_model(times, n_cells_initial, growth_rate):
    
    """
    Return the population 'n' for the input times using an exponential growth model.
    @param times np.array of times to evaluate the growth rate
    
    """
    
    return n_cells_initial * np.exp(growth_rate * times)

def logistic_growth_model(times, n_cells_initial, growth_rate, carrying_capacity):
    
    """
    Return the population 'n' for the input times using a logistic growth model.
    @param times np.array of times to evaluate the growth rate
    
    """
    numerator = carrying_capacity * n_cells_initial * np.exp(growth_rate * times)
    denominator = carrying_capacity + n_cells_initial * (np.exp(growth_rate * times) - 1)
    
    return numerator / denominator

def gompertz_growth_model(times, n_cells_initial, growth_rate, carrying_capacity):
    
    """
    Return the population 'n' for the input times using a gompertz growth model.
    @param times np.array of times to evaluate the growth rate
    
    """
    
    return carrying_capacity * np.exp(np.log(n_cells_initial / carrying_capacity) * np.exp(-growth_rate*times))