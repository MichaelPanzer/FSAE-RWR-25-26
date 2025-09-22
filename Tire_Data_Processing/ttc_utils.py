"""
Contains functions for processing ttc data
"""

import numpy as np

def process_stepwise_data(stepped_data, extr_data = None, threshold=0.2, min_step_length=100):
    """
    Process data by identifying steps, removing transients, and returning step averages and data.
    
    Parameters:
    - stepped_data: array-like, input stepwise data
    - extr_data: array-like, input extra data to be processed using step data
    - threshold: float, relative change threshold to detect a new step
    - min_step_length: int, minimum number of points to consider a stable step
    
    Returns:
    - list of tuples, each containing (step_average, step_data_array)
    """
    if len(stepped_data) == 0:
        return []

    
    diffs = [np.abs(stepped_data[i+1]-stepped_data[i]) for i in range(len(stepped_data)-1)]

    #find all indices of transients
    absolute_threshold = threshold*np.max(diffs)
    transient_indices = np.where(diffs > absolute_threshold)[0] + 1
    transient_indices = np.concatenate(([0], transient_indices, [len(stepped_data)]))

    #find step edges from transient indices 
    step_indices = []
    for i in range(len(transient_indices)-1):
        start = transient_indices[i] 
        end = transient_indices[i + 1]

        # Only include steps with sufficient length
        if (end-start) >= min_step_length:
            step_indices.append((start,end))
    

    #process step indices into list of tuples with average step value and step arrays
    if extr_data != None:
        return [
            (np.mean(step := stepped_data[start:end]), 
            [step, *[d[start:end] for d in extr_data]])
            for start, end in step_indices
        ]    
    
    return [(np.mean(step := stepped_data[start:end]), step) for start, end in step_indices ]
    

def regroup_similar_steps(tuples_array, rel_tol=0.05, absolute_tol=0.05):
    """
    Groups tuples with similar averages and concatenates their data arrays.
    
    Args:
        tuples_array: List of tuples in format [(avg, data_array), ...]
        rel_tol: Relative tolerance for considering averages similar (default: 5%)
        absolute_tol: Absolute tolerance to be used when values are < 1 (default: 0.5)
    
    Returns:
        List of regrouped tuples in same format
    """
    if not tuples_array:
        return []
    
    # Sort by average to group similar values together
    sorted_tuples = sorted(tuples_array, key=lambda x: x[0])
    
    regrouped = []
    current_avg, current_data = sorted_tuples[0]
    
    for avg, data in sorted_tuples[1:]:
        # Check if averages are similar within relative tolerance
        if np.isclose(avg, current_avg, rtol=rel_tol) or np.isclose(current_avg, avg, atol=absolute_tol):
            # Concatenate along first axis (vertically)
            current_data = np.concatenate((current_data, data), axis=1)
        else:
            # Add to results and start new group
            regrouped.append((current_avg, current_data))
            current_avg, current_data = avg, data
    
    # Add the last group
    regrouped.append((current_avg, current_data))
    
    return regrouped

