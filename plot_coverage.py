import sys
import math

def generate_coverage_pseudo_plot(mpileup_file_path, plot_width=80, plot_height=20):
    """
    Reads an mpileup file and generates a text-based pseudo-plot of coverage.

    Args:
        mpileup_file_path (str): Path to the input mpileup file.
        plot_width (int): The maximum width of the plot area (excluding Y-axis labels and separator).
        plot_height (int): The maximum height of the plot area (representing coverage levels).
    """
    # Ensure plot_width is at least 1
    if plot_width < 1:
        print("Error: plot_width must be at least 1.", file=sys.stderr)
        return

    positions = []
    coverages = []
    max_coverage = 0

    print(f"Reading mpileup file: {mpileup_file_path}")

    try:
        with open(mpileup_file_path, 'r') as f:
            for line in f:
                # Skip header lines if any (lines starting with @)
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    try:
                        # Position is the 2nd column (index 1), Coverage is the 4th (index 3)
                        pos = int(fields[1])
                        cov = int(fields[3])
                        positions.append(pos)
                        coverages.append(cov)
                        if cov > max_coverage:
                            max_coverage = cov
                    except ValueError:
                        # Skip lines where position or coverage are not valid integers
                        continue
    except FileNotFoundError:
        print(f"Error: File not found at {mpileup_file_path}", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return

    if not positions:
        print("No valid data found in the mpileup file.", file=sys.stderr)
        return

    print(f"Max coverage found: {max_coverage}")

    # --- Generate Pseudo Plot ---

    # Adjust plot height if max_coverage is less than desired height
    # Use max_coverage + 1 to include the 0 coverage level
    actual_plot_height = min(plot_height, max_coverage + 1)
    if actual_plot_height < 1:
         actual_plot_height = 1 # Ensure at least one row for 0 coverage

    # Create a grid for the plot characters
    # Grid is height x width, initialized with spaces
    plot_grid = [[' ' for _ in range(plot_width)] for _ in range(actual_plot_height)]

    # Determine which positions to sample if there are more positions than plot_width
    num_positions = len(positions)
    indices_to_plot = []
    if num_positions > 0:
        if num_positions > plot_width:
            # Simple sampling: take plot_width evenly spaced points
            indices_to_plot = [int(i * (num_positions - 1) / (plot_width - 1)) for i in range(plot_width)]
        else:
            # Plot all positions if they fit within the width
            indices_to_plot = list(range(num_positions))

    # Populate the plot grid
    for col_idx, data_idx in enumerate(indices_to_plot):
        pos = positions[data_idx]
        cov = coverages[data_idx]

        # Determine the height for this position based on scaled coverage
        # Scale coverage to the plot height (0 to actual_plot_height - 1)
        scaled_height = 0
        if max_coverage > 0 and actual_plot_height > 1:
             scaled_height = int((cov / max_coverage) * (actual_plot_height -1))
        elif max_coverage > 0 and actual_plot_height == 1:
             # If height is 1, plot '#' if coverage > 0
             scaled_height = 0 if cov == 0 else 0 # Always plot at the single row if coverage > 0


        # Fill the column from the bottom up to the scaled height
        # The bottom row of the grid is index actual_plot_height - 1
        for row_idx in range(scaled_height + 1):
             grid_row = actual_plot_height - 1 - row_idx
             if grid_row >= 0: # Ensure we don't go above the grid
                 plot_grid[grid_row][col_idx] = '#' # Use '#' for coverage


    # --- Add Axes (Text Representation) ---

    # Y-axis (Coverage levels)
    y_axis_labels = []
    # Calculate coverage levels for each row from top to bottom
    for i in range(actual_plot_height):
        # Calculate coverage value for this row based on scaling
        if actual_plot_height > 1:
            coverage_level = int(((actual_plot_height - 1 - i) / (actual_plot_height - 1)) * max_coverage)
        else:
            # If height is 1, the only level is max_coverage (or 0 if max_coverage is 0)
            coverage_level = max_coverage if max_coverage > 0 else 0

        y_axis_labels.append(str(coverage_level).ljust(len(str(max_coverage)) + 2)) # Left-align label, pad based on max_coverage length

    # Calculate the width needed for the Y-axis labels and separator
    y_axis_width = len(y_axis_labels[0]) + 3 # Label width + " | "

    # Print the plot with Y-axis labels
    print("\nCoverage Pseudo-Plot:")
    # Separator line width
    separator_width = y_axis_width + plot_width
    print("-" * separator_width)

    for i in range(actual_plot_height):
        # Print Y-axis label and then the plot row
        print(f"{y_axis_labels[i]}| {''.join(plot_grid[i])}") # Print separator directly

    print("-" * separator_width) # Separator line

    # Print X-axis markers below the plot (simplified)
    # This is a basic representation, not a precise scale
    total_line_width = y_axis_width + plot_width
    x_axis_marker_line = [' '] * total_line_width

    # Place the separator below the plot
    separator_start_index = len(y_axis_labels[0])
    # Fill the separator column with '-' or spaces below the plot
    # For the marker line, just put a space and '|'
    x_axis_marker_line[separator_start_index] = ' '
    x_axis_marker_line[separator_start_index + 1] = '|'
    x_axis_marker_line[separator_start_index + 2] = ' '


    # Calculate column indices for markers within the plot area (0 to plot_width - 1)
    # These indices are relative to the start of the plot area (after Y-axis and separator)
    start_col = 0
    mid_col = plot_width // 2
    end_col = plot_width - 1

    # Calculate the actual list indices for placing markers in the total line width
    start_marker_list_idx = y_axis_width + start_col
    mid_marker_list_idx = y_axis_width + mid_col
    end_marker_list_idx = y_axis_width + end_col

    # Place the markers, checking bounds
    if start_marker_list_idx < total_line_width:
        x_axis_marker_line[start_marker_list_idx] = '|'
    if mid_marker_list_idx < total_line_width:
        x_axis_marker_line[mid_marker_list_idx] = '|'
    if end_marker_list_idx < total_line_width:
        x_axis_marker_line[end_marker_list_idx] = '|'


    print("".join(x_axis_marker_line))

    # Print position labels below the markers
    pos_label_line = [' '] * total_line_width

    # Get the positions corresponding to the plotted columns
    first_pos = positions[indices_to_plot[0]] if indices_to_plot else 0
    last_pos = positions[indices_to_plot[-1]] if indices_to_plot else 0
    mid_pos = positions[indices_to_plot[len(indices_to_plot) // 2]] if indices_to_plot else 0


    start_label = str(first_pos)
    mid_label = str(mid_pos)
    end_label = str(last_pos)

    # Calculate the starting list index for placing each label
    # Attempt to center the label below the marker
    start_label_list_idx = start_marker_list_idx - len(start_label) // 2
    mid_label_list_idx = mid_marker_list_idx - len(mid_label) // 2
    end_label_list_idx = end_marker_list_idx - len(end_label) // 2


    # Place the labels, ensuring they don't go out of bounds or overwrite separator
    # Place start label
    for i in range(len(start_label)):
        current_idx = start_label_list_idx + i
        # Ensure index is within bounds and not in the separator column
        if current_idx >= 0 and current_idx < total_line_width and \
           not (current_idx >= separator_start_index and current_idx < separator_start_index + 3):
            pos_label_line[current_idx] = start_label[i]

    # Place mid label
    for i in range(len(mid_label)):
        current_idx = mid_label_list_idx + i
        if current_idx >= 0 and current_idx < total_line_width and \
           not (current_idx >= separator_start_index and current_idx < separator_start_index + 3):
             pos_label_line[current_idx] = mid_label[i]

    # Place end label
    for i in range(len(end_label)):
        current_idx = end_label_list_idx + i
        if current_idx >= 0 and current_idx < total_line_width and \
           not (current_idx >= separator_start_index and current_idx < separator_start_index + 3):
             pos_label_line[current_idx] = end_label[i]


    print("".join(pos_label_line))


# --- Script Entry Point ---
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py <mpileup_file>")
        sys.exit(1)

    mpileup_file = sys.argv[1]
    generate_coverage_pseudo_plot(mpileup_file)
