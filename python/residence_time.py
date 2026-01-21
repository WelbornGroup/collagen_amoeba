########################################################################
###### Functions to calculate residence times
########################################################################
import numpy as np

def frame_waters(linesf,skips,wats):
    # ats = int(float(linesf[0].rstrip().split()[0]))
    watt = np.zeros([wats],dtype=int)
    # watt_coords = np.zeros([wats,3])
    for ii,linef in enumerate(linesf):
        wordsf = linef.rstrip().split()
        if float(wordsf[4]) == 349:
            wintemp = int((int(float(wordsf[0])) - skips)/3)
            # wc_t = np.array([float(wordsf[1]),float(wordsf[2]),wordsf[3]])
            watt[wintemp] = 1
            # watt_coords[wintemp] = wc_t
    return np.array(watt)
    # return np.array(watt),np.array(watt_coords)
    
def calc_residence_times(wt,tol):
    residence_times = []
    for mol_idx in range(wt.shape[1]):
        molecule_data = wt[:, mol_idx]  # Extract time series for molecule
        diff = np.diff(molecule_data.astype(int), prepend=0, append=0)  # Detect transitions
        start_indices = np.where(diff == 1)[0]  # Start of residence events (0→1)
        end_indices = np.where(diff == -1)[0]  # End of residence events (1→0)
        # Adjust for tolerance: merge events with short false sequences
        merged_start = [start_indices[0]] if len(start_indices) > 0 else []
        merged_end = []
        for i in range(1, len(start_indices)):
            prev_end = end_indices[i - 1]
            curr_start = start_indices[i]
            if (curr_start - prev_end) <= tol:  # Merge if gap ≤ tolerance
                continue  # Skip adding a new event
            else:
                merged_end.append(prev_end)  # Close previous event
                merged_start.append(curr_start)  # Start a new event
        if len(start_indices) > 0:
            merged_end.append(end_indices[-1])  # Add last end
        # Compute residence times
        for s, e in zip(merged_start, merged_end):
            residence_times.append(e - s + 1)  # Add 1 to include last frame
    residence_times = np.array(residence_times)
    max_res_time = max(residence_times) if len(residence_times) > 0 else 10  # Handle empty case
    return max_res_time,residence_times

def calculate_exchange_rate(water_table):
    """
    Calculates the exchange rate of water molecules in the hydration shell.
    Parameters:
    - water_table (numpy array): A boolean array of shape (num_frames, num_molecules)
                                 indicating whether a molecule is in the hydration shell.

    Returns:
    - exchange_rate (float): The number of exchange events per frame.
    """
    num_frames, num_molecules = water_table.shape
    total_exchanges = 0  # Counter for exchange events
    # Compute exchanges for each molecule
    for mol_idx in range(num_molecules):
        molecule_data = water_table[:, mol_idx].astype(int)  # Convert boolean to int
        diff = np.diff(molecule_data, prepend=0)  # Detect changes (0->1 or 1->0)
        num_exchanges = np.count_nonzero(diff)  # Count all transition events
        total_exchanges += num_exchanges
    # Compute exchange rate as total exchanges per frame
    exchange_rate = total_exchanges / num_frames
    return exchange_rate

def percentilex(rt, x):
    percentile_x = np.percentile(rt, x)
    filtered_rt = rt[rt <= percentile_x]
    mean_filtered_rt = np.mean(filtered_rt)
    std_filtered_rt = np.std(filtered_rt)
    median_filtered_rt = np.median(filtered_rt)
    return mean_filtered_rt, std_filtered_rt, median_filtered_rt

def glc_indices(fname,type):
    glcs = []
    nglc = 0
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    for i,line in enumerate(lines):
        words = line.rstrip().split()
        if (i > 1):
            if int(words[5]) == type:
                glcs.append(int(words[0])-1)
    return len(glcs), np.array(glcs)

def parse_glucose_presence(filename, c1_indices, total_glucoses):
    """
    Parse a file to create a 2D array indicating glucose presence across frames.
    
    Parameters:
    filename (str): Path to the input file
    c1_indices (list): List of C1 indices for all glucoses
    total_glucoses (int): Total number of glucoses in the system
    
    Returns:
    numpy.ndarray: 2D array of shape [frames, number_of_glucoses] with 1s and 0s
                   indicating presence (1) or absence (0) of each glucose
    """
       
    # Create a mapping from C1 index to glucose position for faster lookup
    c1_to_glucose_idx = {c1_idx: i for i, c1_idx in enumerate(c1_indices)}
    
    frames_data = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    line_idx = 0
    frame_num = 0
    
    while line_idx < len(lines):
        # Initialize current frame array (all zeros)
        current_frame = np.zeros(total_glucoses, dtype=int)
        
        # Read the number of glucoses that fulfill criteria in this frame
        num_glucoses_present = int(lines[line_idx].strip())
        line_idx += 1
        
        # If there are glucoses present, read their details
        if num_glucoses_present > 0:
            for i in range(num_glucoses_present):
                if line_idx < len(lines):
                    # Extract the C1 index (first word) from the detail line
                    detail_line = lines[line_idx].strip()
                    c1_index = int(detail_line.split()[0])
                    
                    # Find which glucose this C1 belongs to and mark it as present
                    if c1_index in c1_to_glucose_idx:
                        glucose_position = c1_to_glucose_idx[c1_index]
                        current_frame[glucose_position] = 1
                    
                    line_idx += 1
        
        frames_data.append(current_frame)
        frame_num += 1
    
    # Convert to 2D numpy array
    result_array = np.array(frames_data)
    
    print(f"Processed {len(frames_data)} frames")
    print(f"Result array shape: {result_array.shape}")
    
    return result_array

def analyze_coordination_events(presence_array, threshold=0):
    """
    Analyze coordination events from glucose presence array.
    
    Parameters:
    presence_array (numpy.ndarray): 2D array of shape [frames, glucoses] with 1s and 0s
    threshold (int): Maximum number of consecutive 0s to merge events (default: 0)
    
    Returns:
    dict: Dictionary containing all analysis results
    """
    frames, num_glucoses = presence_array.shape
    
    # Arrays to store results
    events_per_glucose = np.zeros(num_glucoses, dtype=int)
    all_residence_times = []
    residence_times_per_glucose = [[] for _ in range(num_glucoses)]
    
    # Process each glucose separately
    for glucose_idx in range(num_glucoses):
        glucose_timeline = presence_array[:, glucose_idx]
        events, residence_times = find_coordination_events(glucose_timeline, threshold)
        
        events_per_glucose[glucose_idx] = events
        all_residence_times.extend(residence_times)
        residence_times_per_glucose[glucose_idx] = residence_times
    
    # Calculate statistics
    total_events = np.sum(events_per_glucose)
    
    if all_residence_times:
        mean_residence_time = np.mean(all_residence_times)
        median_residence_time = np.median(all_residence_times)
        residence_distribution = np.array(all_residence_times)
        
        # Calculate exchange rates (events per frame)
        mean_exchange_rate = total_events / frames if frames > 0 else 0
        exchange_rates_per_glucose = events_per_glucose / frames if frames > 0 else np.zeros(num_glucoses)
    else:
        mean_residence_time = 0
        median_residence_time = 0
        residence_distribution = np.array([])
        mean_exchange_rate = 0
        exchange_rates_per_glucose = np.zeros(num_glucoses)
    
    return {
        'events_per_glucose': events_per_glucose,
        'total_events': total_events,
        'mean_residence_time': mean_residence_time,
        'median_residence_time': median_residence_time,
        'residence_distribution': residence_distribution,
        'residence_times_per_glucose': residence_times_per_glucose,
        'mean_exchange_rate': mean_exchange_rate,
        'exchange_rates_per_glucose': exchange_rates_per_glucose,
        'threshold_used': threshold
    }

def find_coordination_events(timeline, threshold):
    """
    Find coordination events in a single glucose timeline with gap merging.
    
    Parameters:
    timeline (numpy.ndarray): 1D array of 1s and 0s for a single glucose
    threshold (int): Maximum gap size to merge events
    
    Returns:
    tuple: (number_of_events, list_of_residence_times)
    """
    if len(timeline) == 0:
        return 0, []
    
    # Find all sequences of 1s
    events = []
    in_event = False
    start_idx = 0
    
    for i, value in enumerate(timeline):
        if value == 1 and not in_event:
            # Start of new event
            in_event = True
            start_idx = i
        elif value == 0 and in_event:
            # End of event
            in_event = False
            events.append((start_idx, i - 1))
    
    # Handle case where timeline ends with 1s
    if in_event:
        events.append((start_idx, len(timeline) - 1))
    
    if not events:
        return 0, []
    
    # Merge events based on threshold
    merged_events = merge_events_with_threshold(events, threshold)
    
    # Calculate residence times
    residence_times = []
    for start, end in merged_events:
        residence_time = end - start + 1  # +1 because both endpoints are inclusive
        residence_times.append(residence_time)
    
    return len(merged_events), residence_times

def merge_events_with_threshold(events, threshold):
    """
    Merge events if the gap between them is <= threshold.
    
    Parameters:
    events (list): List of (start, end) tuples
    threshold (int): Maximum gap size to merge
    
    Returns:
    list: List of merged (start, end) tuples
    """
    if not events or threshold < 0:
        return events
    
    merged = []
    current_start, current_end = events[0]
    
    for i in range(1, len(events)):
        next_start, next_end = events[i]
        gap = next_start - current_end - 1  # Gap between events
        
        if gap <= threshold:
            # Merge events by extending current_end
            current_end = next_end
        else:
            # Gap too large, finalize current event and start new one
            merged.append((current_start, current_end))
            current_start, current_end = next_start, next_end
    
    # Add the final event
    merged.append((current_start, current_end))
    
    return merged

def print_coordination_analysis(results):
    """
    Print a comprehensive summary of coordination analysis results.
    
    Parameters:
    results (dict): Results dictionary from analyze_coordination_events
    """
    print(f"=== Coordination Events Analysis (Threshold: {results['threshold_used']}) ===")
    print(f"\nTotal Events: {results['total_events']}")
    print(f"Mean Exchange Rate: {results['mean_exchange_rate']:.4f} events/frame")
    
    if results['total_events'] > 0:
        print(f"\nResidence Time Statistics:")
        print(f"  Mean: {results['mean_residence_time']:.2f} frames")
        print(f"  Median: {results['median_residence_time']:.2f} frames")
        print(f"  Min: {np.min(results['residence_distribution']):.2f} frames")
        print(f"  Max: {np.max(results['residence_distribution']):.2f} frames")
        print(f"  Std Dev: {np.std(results['residence_distribution']):.2f} frames")
    
    print(f"\nEvents per Glucose:")
    for i, events in enumerate(results['events_per_glucose']):
        rate = results['exchange_rates_per_glucose'][i]
        # print(f"  Glucose {i}: {events} events (rate: {rate:.4f} events/frame)")
    
    print(f"\nResidence Time Distribution:")
    print(f"  Total residence times recorded: {len(results['residence_distribution'])}")
    if len(results['residence_distribution']) > 0:
        unique_times, counts = np.unique(results['residence_distribution'], return_counts=True)
        print(f"  Unique residence times: {len(unique_times)}")
        print(f"  Most common residence time: {unique_times[np.argmax(counts)]} frames ({np.max(counts)} occurrences)")