####################################################################
# Functions to calculate Pyramidalization angles and distances
####################################################################
import numpy as np
import re

def parse_atomic_coordinates(filename):
    """
    Parse atomic coordinates file and return arrays organized by atom type and strand.
    
    Args:
        filename (str): Path to the input file
    
    Returns:
        dict: Dictionary containing 12 arrays (4 atom types × 3 strands)
              Keys: 'N_strand1', 'N_strand2', 'N_strand3', 
                    'CA_strand1', 'CA_strand2', 'CA_strand3',
                    'C_strand1', 'C_strand2', 'C_strand3',
                    'O_strand1', 'O_strand2', 'O_strand3'
              Values: numpy arrays of shape (num_frames, atoms_per_strand, 3)
    """
    
    frames = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    i = 0
    while i < len(lines):
        # Read number of lines in this frame
        num_lines = int(lines[i].strip())
        i += 1
        
        frame_coords = []
        
        # Read the coordinates for this frame
        for j in range(num_lines):
            line = lines[i + j].strip()
            
            # Parse coordinates using regex to find all [x,y,z] patterns
            coord_pattern = r'\[([^\]]+)\]'
            matches = re.findall(coord_pattern, line)
            
            line_coords = []
            for match in matches:
                # Split by comma and convert to float
                coords = [float(x.strip()) for x in match.split(',')]
                line_coords.append(coords)
            
            frame_coords.extend(line_coords)
        
        frames.append(frame_coords)
        i += num_lines
    
    # Convert to numpy array
    frames = np.array(frames)  # Shape: (num_frames, total_atoms, 3)
    
    num_frames = frames.shape[0]
    total_atoms = frames.shape[1]
    
    # Assuming atoms are organized as [N, CA, C, O, N, CA, C, O, ...]
    # We need to separate by atom type first
    atoms_per_residue = 4  # N, CA, C, O
    num_residues = total_atoms // atoms_per_residue
    atoms_per_strand = num_residues // 3
    
    # Initialize result dictionary
    result = {}
    atom_types = ['N', 'CA', 'C', 'O']
    
    for atom_idx, atom_type in enumerate(atom_types):
        # Extract all atoms of this type across all residues
        atom_coords = []
        for residue in range(num_residues):
            atom_position = residue * atoms_per_residue + atom_idx
            atom_coords.append(frames[:, atom_position, :])
        
        atom_coords = np.array(atom_coords)  # Shape: (num_residues, num_frames, 3)
        atom_coords = np.transpose(atom_coords, (1, 0, 2))  # Shape: (num_frames, num_residues, 3)
        
        # Split into three strands
        for strand in range(3):
            start_idx = strand * atoms_per_strand
            end_idx = (strand + 1) * atoms_per_strand
            strand_coords = atom_coords[:, start_idx:end_idx, :]
            
            key = f"{atom_type}_strand{strand + 1}"
            result[key] = strand_coords
    
    return result

def calculate_protein_metrics(N_coords, CA_coords, C_coords, O_coords):
    """
    Calculate protein structural metrics for a single strand.
    
    Args:
        N_coords: numpy array of shape (num_frames, num_atoms, 3)
        CA_coords: numpy array of shape (num_frames, num_atoms, 3)
        C_coords: numpy array of shape (num_frames, num_atoms, 3)
        O_coords: numpy array of shape (num_frames, num_atoms, 3)
    
    Returns:
        tuple: (oi_ci_distances, oci_angles, plane_vector_angles)
               Each array has shape (num_frames, num_atoms)
    """
    num_frames, num_atoms, _ = N_coords.shape
    
    # Initialize output arrays
    oi_ci_distances = np.zeros((num_frames, num_atoms))
    oci_angles = np.zeros((num_frames, num_atoms))
    plane_vector_angles = np.zeros((num_frames, num_atoms))
    
    # 1. O(i-1) distance to C(i) over time - set i=0 to 0
    for i in range(1, num_atoms):
        # Distance between O(i-1) and C(i)
        diff = O_coords[:, i-1, :] - C_coords[:, i, :]
        oi_ci_distances[:, i] = np.linalg.norm(diff, axis=1)
    
    # 2. O(i-1)-C(i)-O(i) angle over time - set i=0 to 0
    for i in range(1, num_atoms):
        # Vectors from C(i) to O(i-1) and from C(i) to O(i)
        vec1 = O_coords[:, i-1, :] - C_coords[:, i, :]  # C(i) to O(i-1)
        vec2 = O_coords[:, i, :] - C_coords[:, i, :]    # C(i) to O(i)
        
        # Normalize vectors
        vec1_norm = vec1 / (np.linalg.norm(vec1, axis=1, keepdims=True) + 1e-10)
        vec2_norm = vec2 / (np.linalg.norm(vec2, axis=1, keepdims=True) + 1e-10)
        
        # Calculate angle using dot product
        dot_product = np.sum(vec1_norm * vec2_norm, axis=1)
        dot_product = np.clip(dot_product, -1.0, 1.0)  # Avoid numerical errors
        oci_angles[:, i] = np.arccos(dot_product) * 180.0 / np.pi  # Convert to degrees
    
    # 3. Angle between CA(i)-C(i)-N(i+1) plane and C(i)-O(i) vector - set last atom to 0
    for i in range(num_atoms - 1):
        # Define vectors for the plane CA(i)-C(i)-N(i+1)
        vec_ca_c = C_coords[:, i, :] - CA_coords[:, i, :]     # CA(i) to C(i)
        vec_c_n = N_coords[:, i+1, :] - C_coords[:, i, :]    # C(i) to N(i+1)
        
        # Calculate plane normal using cross product
        plane_normal = np.cross(vec_ca_c, vec_c_n)
        
        # Normalize plane normal
        plane_normal_norm = plane_normal / (np.linalg.norm(plane_normal, axis=1, keepdims=True) + 1e-10)
        
        # C(i)-O(i) vector
        co_vector = O_coords[:, i, :] - C_coords[:, i, :]
        co_vector_norm = co_vector / (np.linalg.norm(co_vector, axis=1, keepdims=True) + 1e-10)
        
        # Calculate angle between plane normal and C-O vector
        # The angle between plane and vector is 90° - angle between normal and vector
        dot_product = np.sum(plane_normal_norm * co_vector_norm, axis=1)
        dot_product = np.clip(dot_product, -1.0, 1.0)
        angle_with_normal = np.arccos(np.abs(dot_product)) * 180.0 / np.pi
        plane_vector_angles[:, i] = 90.0 - angle_with_normal
    
    return oi_ci_distances, oci_angles, plane_vector_angles

def score_metrics(distances, angles, dihedrals, mean_distance, mean_angle, mean_dihedral, threshold=0.1):
    """
    Score protein metrics based on deviation from expected mean values.
    
    Args:
        distances: numpy array of shape (num_frames, num_atoms) - distance values
        angles: numpy array of shape (num_frames, num_atoms) - angle values
        dihedrals: numpy array of shape (num_frames, num_atoms) - dihedral angle values
        mean_distance: float - expected mean distance value
        mean_angle: float - expected mean angle value
        mean_dihedral: float - expected mean dihedral value
        threshold: float - threshold percentage (default 0.1 for 10%)
    
    Returns:
        tuple: (distance_scores, angle_scores, dihedral_scores)
               Each array has shape (num_frames, num_atoms) with scores 0-1
    """
    
    def calculate_scores(values, mean_value):
        """Calculate scores for a single metric type."""
        # Calculate absolute percentage difference from mean
        abs_diff = np.abs(values - mean_value)
        percent_diff = abs_diff / np.abs(mean_value)
        
        # Initialize scores array
        scores = np.zeros_like(values)
        
        # Values within threshold get scaled scores (1 at mean, 0 at threshold)
        within_threshold = percent_diff <= threshold
        scores[within_threshold] = 1.0 - (percent_diff[within_threshold] / threshold)
        
        # Values beyond threshold get score of 0 (already initialized to 0)
        
        return scores
    
    # Calculate scores for each metric type
    distance_scores = calculate_scores(distances, mean_distance)
    angle_scores = calculate_scores(angles, mean_angle)
    dihedral_scores = calculate_scores(dihedrals, mean_dihedral)
    
    return distance_scores, angle_scores, dihedral_scores

def score_metrics_new(distances, angles, dihedrals, mean_distance, mean_angle, mean_dihedral, threshold=0.1):
    """
    Score protein metrics based on deviation from expected mean values.
    Values beyond threshold get NaN instead of 0.
    
    Args:
        distances: numpy array of shape (num_frames, num_atoms) - distance values
        angles: numpy array of shape (num_frames, num_atoms) - angle values
        dihedrals: numpy array of shape (num_frames, num_atoms) - dihedral angle values
        mean_distance: float - expected mean distance value
        mean_angle: float - expected mean angle value
        mean_dihedral: float - expected mean dihedral value
        threshold: float - threshold percentage (default 0.1 for 10%)
    
    Returns:
        tuple: (distance_scores, angle_scores, dihedral_scores)
               Each array has shape (num_frames, num_atoms) with scores 0-1 or NaN
    """
    
    def calculate_scores(values, mean_value):
        """Calculate scores for a single metric type."""
        # Calculate absolute percentage difference from mean
        abs_diff = np.abs(values - mean_value)
        percent_diff = abs_diff / np.abs(mean_value)
        
        # Initialize scores array with NaN
        scores = np.full_like(values, np.nan)
        
        # Values within threshold get scaled scores (1 at mean, 0 at threshold)
        within_threshold = percent_diff <= threshold
        scores[within_threshold] = 1.0 - (percent_diff[within_threshold] / threshold)
        
        # Values beyond threshold remain as NaN (already initialized)
        
        return scores
    
    # Calculate scores for each metric type
    distance_scores = calculate_scores(distances, mean_distance)
    angle_scores = calculate_scores(angles, mean_angle)
    dihedral_scores = calculate_scores(dihedrals, mean_dihedral)
    
    return distance_scores, angle_scores, dihedral_scores

def count_zero_nonzero_scores(distance_scores, angle_scores, dihedral_scores):
    """
    Count zero and non-zero values in each score array.
    
    Args:
        distance_scores: numpy array of shape (num_frames, num_atoms) - distance scores
        angle_scores: numpy array of shape (num_frames, num_atoms) - angle scores
        dihedral_scores: numpy array of shape (num_frames, num_atoms) - dihedral scores
    
    Returns:
        dict: Dictionary containing counts for each metric type
    """
    
    # Count zeros and non-zeros for distance scores
    distance_zeros = np.sum(distance_scores == 0)
    distance_nonzeros = np.sum(distance_scores > 0)
    
    # Count zeros and non-zeros for angle scores
    angle_zeros = np.sum(angle_scores == 0)
    angle_nonzeros = np.sum(angle_scores > 0)
    
    # Count zeros and non-zeros for dihedral scores
    dihedral_zeros = np.sum(dihedral_scores == 0)
    dihedral_nonzeros = np.sum(dihedral_scores > 0)
    
    # Calculate total elements for verification
    total_elements = distance_scores.size
    
    results = {
        'distance': {
            'zeros': distance_zeros,
            'non_zeros': distance_nonzeros,
            'total': total_elements,
            'zero_percentage': (distance_zeros / total_elements) * 100,
            'non_zero_percentage': (distance_nonzeros / total_elements) * 100
        },
        'angle': {
            'zeros': angle_zeros,
            'non_zeros': angle_nonzeros,
            'total': total_elements,
            'zero_percentage': (angle_zeros / total_elements) * 100,
            'non_zero_percentage': (angle_nonzeros / total_elements) * 100
        },
        'dihedral': {
            'zeros': dihedral_zeros,
            'non_zeros': dihedral_nonzeros,
            'total': total_elements,
            'zero_percentage': (dihedral_zeros / total_elements) * 100,
            'non_zero_percentage': (dihedral_nonzeros / total_elements) * 100
        }
    }
    
    return results

def print_score_counts(results):
    """
    Print the score count results in a readable format.
    
    Args:
        results: Dictionary returned by count_zero_nonzero_scores
    """
    print("Score Count Summary:")
    print("=" * 50)
    
    for metric in ['distance', 'angle', 'dihedral']:
        data = results[metric]
        print(f"\n{metric.capitalize()} Scores:")
        print(f"  Zeros: {data['zeros']:,} ({data['zero_percentage']:.1f}%)")
        print(f"  Non-zeros: {data['non_zeros']:,} ({data['non_zero_percentage']:.1f}%)")
        print(f"  Total: {data['total']:,}")

def calculate_binary_conditions(distance_scores, angle_scores, dihedral_scores):
    """
    Calculate binary condition arrays based on non-zero score combinations.
    
    Args:
        distance_scores: numpy array of shape (num_frames, num_atoms) - distance scores
        angle_scores: numpy array of shape (num_frames, num_atoms) - angle scores
        dihedral_scores: numpy array of shape (num_frames, num_atoms) - dihedral scores
    
    Returns:
        tuple: (distance_angle_binary, distance_dihedral_binary, all_three_binary)
               Each array has shape (num_frames, num_atoms) with values 0 or 1
    """
    
    # Case 1: Distance and angle scores are both non-zero
    distance_angle_binary = ((distance_scores > 0) & (angle_scores > 0)).astype(int)
    
    # Case 2: Distance and dihedral scores are both non-zero
    distance_dihedral_binary = ((distance_scores > 0) & (dihedral_scores > 0)).astype(int)
    
    # Case 3: All three scores are non-zero
    all_three_binary = ((distance_scores > 0) & (angle_scores > 0) & (dihedral_scores > 0)).astype(int)
    
    return distance_angle_binary, distance_dihedral_binary, all_three_binary
